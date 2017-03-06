#include "streakline.h"

// Simulation parameters
double Mcli, Mclf, Rcl, dt;

int stream(double *x0, double *v0, double *xm1, double *xm2, double *xm3, double *xp1, double *xp2, double *xp3, double *vm1, double *vm2, double *vm3, double *vp1, double *vp2, double *vp3, double *par, double *offset, int potential, int integrator, int N, int M, double mcli, double mclf, double rcl, double dt_)
{
	int i,j, k=0, Napar, Ne, imin=0;
	double x[3], v[3], xs[3], vs[3], omega[3], om, sign=1., back=-1., r, rp, rm, vtot, vlead, vtrail, dM, Mcl, dR, dRRj;
	double *xc1, *xc2, *xc3, *Rj, *dvl, *dvt;
	long s1=560;
	
	// number of output particles
	Ne=ceil((float)N/(float)M);
	xc1 = (double*) malloc(N*sizeof(double));
	xc2 = (double*) malloc(N*sizeof(double));
	xc3 = (double*) malloc(N*sizeof(double));
	Rj = (double*) malloc(Ne*sizeof(double));
	dvl = (double*) malloc(Ne*sizeof(double));
	dvt = (double*) malloc(Ne*sizeof(double));
	
	// Initial position and velocities
	t2t(x0,x);
	t2t(v0,v);
	
	// Initial mass and mass loss rate
	Mcli = mcli;
	Mclf = mclf;
	Mcl = Mcli;
	dM = (Mcli-Mclf)/(double)N;
	
	// Cluster size
	Rcl =rcl;
	
	// Position offset
	dR = offset[0];
	
	// Initialize velocity offsets, drawn from a Maxwell distribution
	double r1,r2,r3;
	for(i=0;i<Ne;i++){
		// Leading tail velocity offsets
		r1=gasdev(&s1);
		r2=gasdev(&s1);
		r3=gasdev(&s1);
		dvl[i]=sqrt(r1*r1 + r2*r2 + r3*r3)*offset[1]/3.;
		
		// Trailing tail velocity offsets
		r1=gasdev(&s1);
		r2=gasdev(&s1);
		r3=gasdev(&s1);
		dvt[i]=sqrt(r1*r1 + r2*r2 + r3*r3)*offset[1]/3.;
	}
	
	// Set up actual potential parameters;
	if(potential==1){
		Napar=3;
	}else if(potential==2 || potential==3){
		Napar=6;
	}else if(potential==4){
		Napar=11;
	}else if(potential==5){
		Napar=4;
	}else if(potential==6){
        Napar=15;
    }else if(potential==7){
        Napar=13;
    }else{
		Napar=1;
	}
	double apar[Napar];
	initpar(potential, par, apar);
	
	// Integrator switch
	void (*pt2dostep)(double*, double*, double*, int, double, double) = NULL;
	
	if(integrator==0){
		// Leapfrog
		pt2dostep=&dostep;
	}
	else if(integrator==1){
		// Runge-Kutta
		pt2dostep=&dostep_rk;
	}
	
	// Time step
	dt = dt_;
	
	///////////////////////////////////////
	// Backward integration (cluster only)
    
	if(integrator==0){
// 		printf("%f\n", v[0]);
		dostep1(x,v,apar,potential,dt,back);
// 		printf("%f\n", v[0]);
		imin=1;
	}
	for(i=imin;i<N;i++){
		(*pt2dostep)(x,v,apar,potential,dt,back);
	}
	if(integrator==0)
		dostep1(x,v,apar,potential,dt,back);
	
	////////////////////////////////////////////
	// Forward integration (cluster and stream)
	
	// Initial step for the leapfrog integrator
	if (integrator==0){
		dostep1(x,v,apar,potential,dt,sign);
		for(j=0;j<3;j++)
			x[j]=x[j]-dt*v[j];
	
		dostep(x,v,apar,potential,dt,sign);
		imin=1;
		
		// Update output arrays
		t2n(x, xc1, xc2, xc3, 0);
		Rj[k]=jacobi(x, v, apar, potential, Mcl);	// Jacobi radius
		r=len(x);
		rm=(r-Rj[k])/r;
		rp=(r+Rj[k])/r;
		
		// Angular velocity
		omega[0]=x[1]*v[2]-x[2]*v[1];
		omega[1]=x[2]*v[0]-x[0]*v[2];
		omega[2]=x[0]*v[1]-x[1]*v[0];
		om=len(omega)/(r*r);
		vtot=len(v);
		vlead=(vtot-om*Rj[0])/vtot;
		vtrail=(vtot+om*Rj[0])/vtot;
		
		dvl[k]/=r;
		dvt[k]/=r;
		
		// Inner particle
		xm1[k]=x[0]*rm + dR*Rj[k]*gasdev(&s1);
		xm2[k]=x[1]*rm + dR*Rj[k]*gasdev(&s1);
		xm3[k]=x[2]*rm + dR*Rj[k]*gasdev(&s1);
		vm1[k]=v[0]*vlead - dvl[k]*x[0];
		vm2[k]=v[1]*vlead - dvl[k]*x[1];
		vm3[k]=v[2]*vlead - dvl[k]*x[2];
		
		// Outer particle
		xp1[k]=x[0]*rp + dR*Rj[k]*gasdev(&s1);
		xp2[k]=x[1]*rp + dR*Rj[k]*gasdev(&s1);
		xp3[k]=x[2]*rp + dR*Rj[k]*gasdev(&s1);
		vp1[k]=v[0]*vtrail + dvt[k]*x[0];
		vp2[k]=v[1]*vtrail + dvt[k]*x[1];
		vp3[k]=v[2]*vtrail + dvt[k]*x[2];
		k++;
	}

	// Subsequent steps
	for(i=imin;i<N;i++){
		Mcl-=dM;
		
		(*pt2dostep)(x,v,apar,potential,dt,sign);

        // Store cluster position
		t2n(x, xc1, xc2, xc3, i);
		
		// Propagate previously released stream particles
		for(j=0;j<k;j++){
			// Inner particle
			n2t(xs, xm1, xm2, xm3, j);
			n2t(vs, vm1, vm2, vm3, j);
			dostep_stream(x,xs,vs,apar,potential,Mcl,dt,sign);
// 			(*pt2dostep)(xs,vs,apar,potential,dt,sign);
			
			// Update
			t2n(xs, xm1, xm2, xm3, j);
			t2n(vs, vm1, vm2, vm3, j);
			
			// Outer particle
			n2t(xs, xp1, xp2, xp3, j);
			n2t(vs, vp1, vp2, vp3, j);
			dostep_stream(x,xs,vs,apar,potential,Mcl,dt,sign);
// 			(*pt2dostep)(xs,vs,apar,potential,dt,sign);
			
			// Update
			t2n(xs, xp1, xp2, xp3, j);
			t2n(vs, vp1, vp2, vp3, j);
		}
		
		if(i%M==0){
			// Release only at every Mth timestep
			// Jacobi tidal radius
			Rj[k]=jacobi(x, v, apar, potential, Mcl);
			r=len(x);
			rm=(r-Rj[k])/r;
			rp=(r+Rj[k])/r;
			
			// Angular velocity
			omega[0]=x[1]*v[2]-x[2]*v[1];
			omega[1]=x[2]*v[0]-x[0]*v[2];
			omega[2]=x[0]*v[1]-x[1]*v[0];
			om=len(omega)/(r*r);
			vtot=len(v);
			vlead=(vtot-om*Rj[k])/vtot;
			vtrail=(vtot+om*Rj[k])/vtot;
			
			dvl[k]/=r;
			dvt[k]/=r;
			
			// Generate 2 new stream particles at the tidal radius
			dRRj = dR*Rj[k];
			// Inner particle (leading tail)
			xm1[k]=x[0]*rm + dRRj*gasdev(&s1);
			xm2[k]=x[1]*rm + dRRj*gasdev(&s1);
			xm3[k]=x[2]*rm + dRRj*gasdev(&s1);
			vm1[k]=v[0]*vlead - dvl[k]*x[0];
			vm2[k]=v[1]*vlead - dvl[k]*x[1];
			vm3[k]=v[2]*vlead - dvl[k]*x[2];
			
			// Outer particle (trailing tail)
			xp1[k]=x[0]*rp + dRRj*gasdev(&s1);
			xp2[k]=x[1]*rp + dRRj*gasdev(&s1);
			xp3[k]=x[2]*rp + dRRj*gasdev(&s1);
			vp1[k]=v[0]*vtrail + dvt[k]*x[0];
			vp2[k]=v[1]*vtrail + dvt[k]*x[1];
			vp3[k]=v[2]*vtrail + dvt[k]*x[2];
			k++;
		}
	}
	
	if (integrator==0)
		dostep1(x,v,apar,potential,dt,back);
	
	// Free memory
	free(xc1);
	free(xc2);
	free(xc3);
	free(Rj);
	free(dvl);
	free(dvt);
	
	return 0;
} 

int orbit(double *x0, double *v0, double *x1, double *x2, double *x3, double *v1, double *v2, double *v3, double *par, int potential, int integrator, int N, double dt_, double direction)
{
	int i, Napar, imin=0;
	double x[3], v[3];
	
	// Initial position and velocities
	t2t(x0, x);
	t2t(v0, v);
	
    // Set up actual potential parameters;
    if(potential==1){
        Napar=3;
    }else if(potential==2 || potential==3){
        Napar=6;
    }else if(potential==4){
        Napar=11;
    }else if(potential==5){
        Napar=4;
    }else if(potential==6){
        Napar=15;
    }else if(potential==7){
        Napar=13;
    }else{
        Napar=1;
    }
	double apar[Napar];
	initpar(potential, par, apar);
	
	// Integrator switch
	void (*pt2dostep)(double*, double*, double*, int, double, double) = NULL;
	
	if(integrator==0){
		// Leapfrog
		pt2dostep=&dostep;
	}
	else if(integrator==1){
		// Runge-Kutta
		pt2dostep=&dostep_rk;
	}
	
	// Time step
	dt = dt_;
	
	///////////////////////////////////////
	// Backward integration (cluster only)
	
	if(integrator==0){
		dostep1(x,v,apar,potential,dt,direction);
		
		// Record
		t2n(x, x1, x2, x3, 0);
		t2n(v, v1, v2, v3, 0);
		imin=1;
	}
	for(i=imin;i<N;i++){
		(*pt2dostep)(x,v,apar,potential,dt,direction);
		
		// Record
		t2n(x, x1, x2, x3, i);
		t2n(v, v1, v2, v3, i);
	}
	if(integrator==0){
		dostep1(x,v,apar,potential,dt,direction);
		
		// Record
		t2n(x, x1, x2, x3, N-1);
		t2n(v, v1, v2, v3, N-1);
	}
	
	return 0;
}

void dostep(double *x, double *v, double *par, int potential, double deltat, double sign)
{	// Evolve point particle from x0, v0, for a time deltat in a given potential
	// evolve forward for sign=1, backwards for sign=-1
	// return final positions and velocities in x, v
	
	int i, j, Nstep;
	double xt[3], vt[3], at[3], dts;
	
	dts=sign*dt;			// Time step with a sign
	Nstep=(int) (deltat/dt);	// Number of steps to evolve

	for(i=0;i<Nstep;i++){
		// Forward the particle using the leapfrog integrator
		for(j=0;j<3;j++)
			xt[j]=x[j]+dts*v[j];
		force(xt, at, par, potential);
//         if(potential==7) par[12]+=dts;
		for(j=0;j<3;j++)
			vt[j]=v[j]+dts*at[j];
		
		// Update input vectors to current values
		for(j=0;j<3;j++){
			x[j]=xt[j];
			v[j]=vt[j];
		}
	}
	
}

void dostep1(double *x, double *v, double *par, int potential, double deltat, double sign)
{	// Make first step to set up the leapfrog integration
	
	double a[3], dts;
	
	dts=sign*dt;
	force(x, a, par, potential);
//     if(potential==7) par[12]+=dts;

	v[0]=v[0]+0.5*dts*a[0];
	v[1]=v[1]+0.5*dts*a[1];
	v[2]=v[2]+0.5*dts*a[2];
}

void dostep_rk(double *x, double *v, double *par, int potential, double deltat, double sign)
{	// Evolve point particle from x0, v0, for a time deltat in a given potential
	// evolve forward for sign=1, backwards for sign=-1
	// return final positions and velocities in x, v
	// prototype for Runge-Kutta integrator
	
	int i;
	double xt1[3], xt2[3], xt3[3], vt1[3], vt2[3], vt3[3], a[3], at1[3], at2[3], at3[3], dts, dt2;
	
	dts=sign*dt;			// Time step with a sign
// 	if(deltat<dt){
// 		dts=sign*deltat;
// 	}
	dt2=dts/2.;
	
	// Initial values
	force(x, a, par, potential);
	
	// First half-step
	for(i=0;i<3;i++){
		xt1[i]=x[i]+dt2*v[i];
		vt1[i]=v[i]+dt2*a[i];
	}
	force(xt1,at1,par,potential);
	
	// Second half-step
	for(i=0;i<3;i++){
		xt2[i]=x[i]+dt2*vt1[i];
		vt2[i]=v[i]+dt2*at1[i];
	}
	force(xt2,at2,par,potential);
	
	// Third step
	for(i=0;i<3;i++){
		xt3[i]=x[i]+dts*vt2[i];
		vt3[i]=v[i]+dts*at2[i];
	}
	force(xt3,at3,par,potential);
	
	// Final Runge-Kutta evaluation
	for(i=0;i<3;i++){
		x[i]+=dts/6.*(v[i]+2.*(vt1[i]+vt2[i])+vt3[i]);
		v[i]+=dts/6.*(a[i]+2.*(at1[i]+at2[i])+at3[i]);
	}
	
}

void dostep_stream(double *xc, double *x, double *v, double *par, int potential, double Mcl, double deltat, double sign)
{	// Same as dostep, except that stream particles also feel the Plummer potential from a cluster
	
	int i, j, Nstep;
	double xt[3], vt[3], at[3], xr[3], ar[3], dts;
	
	dts=sign*dt;			// Time step with a sign
	Nstep=(int) (deltat/dt);	// Number of steps to evolve

	for(i=0;i<Nstep;i++){
		// Forward the particle using the leapfrog integrator
		for(j=0;j<3;j++){
			xt[j]=x[j]+dts*v[j];
			xr[j]=xc[j]-xt[j];
		}
		force(xt, at, par, potential);
//         if(potential==7) par[12]+=dts;
        force_plummer(xr,ar,Mcl);
		for(j=0;j<3;j++)
			vt[j]=v[j]+dts*(at[j]+ar[j]);
		
		// Update input vectors to current values
		for(j=0;j<3;j++){
			x[j]=xt[j];
			v[j]=vt[j];
		}
	}
}

void force(double *x, double *a, double *par, int potential)
{
	int i;
	double r, aux, aux2;
	
	if(potential==0){
		// Point mass potential
		// par = [Mtot]
		r=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
		
		for(i=0;i<3;i++)
			a[i]=-par[0]*x[i]/(r*r*r);
		
	}else if(potential==1){
		// Logarithmic potential, as defined by Koposov et al. (2010)
		// par = [Vc, q, q^2]
		r=x[0]*x[0] + x[1]*x[1] + x[2]*x[2]/(par[2]);
		aux=-par[0]/r;
		
		a[0]=aux*x[0];
		a[1]=aux*x[1];
		a[2]=aux*x[2]/(par[1]);
		
	}else if(potential==2){
		// Triaxial logarithmic halo potential from Law & Majewski (2010)
		// par = [Vc^2, c1, c2, c3, c4, rhalo^2]
		r=par[1]*x[0]*x[0] + par[2]*x[1]*x[1] + par[3]*x[0]*x[1] + par[4]*x[2]*x[2] + par[5];
		aux=-par[0]/r;
		
		a[0]=aux*(2*par[1]*x[0] + par[3]*x[1]);
		a[1]=aux*(2*par[2]*x[1] + par[3]*x[0]);
		a[2]=aux*(2*par[4]*x[2]);
		
	}else if(potential==3){
		// Triaxial NFW halo potential, parameters similar to Law & Majewski (2010)
		// par = [GM, c1, c2, c3, c4, rhalo]
		r=sqrt(par[1]*x[0]*x[0] + par[2]*x[1]*x[1] + par[3]*x[0]*x[1] + par[4]*x[2]*x[2]);
		aux=0.5 * par[0]*pow(r,-3) * (1./(1.+par[5]/r)-log(1.+r/par[5]));
		
		a[0]=aux*(2*par[1]*x[0] + par[3]*x[1]);
		a[1]=aux*(2*par[2]*x[1] + par[3]*x[0]);
		a[2]=aux*(2*par[4]*x[2]);
		
	}else if(potential==4){
		// Composite Galactic potential featuring a disk, bulge, and flattened NFW halo (from Johnston/Law/Majewski/Helmi)
		// par = [GMb, ab, GMd, ad, bd^2, GM, c1, c2, c3, c4, rhalo]
		
		//Hernquist bulge
		r=sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
		aux=-par[0]/(r * (r+par[1]) * (r+par[1]));
		
		a[0]=aux*x[0];
		a[1]=aux*x[1];
		a[2]=aux*x[2];
		
		//Miyamoto-Nagai disk
		aux2=sqrt(x[2]*x[2] + par[4]);
		r=sqrt(x[0]*x[0] + x[1]*x[1] + (par[3] + aux2) * (par[3] + aux2));
		aux=-par[2]/(r*r*r);
		
		a[0]+=aux*x[0];
		a[1]+=aux*x[1];
		a[2]+=aux*x[2]*(par[3] + aux2)/aux2;
		
		//Triaxial NFW Halo
		r=sqrt(par[6]*x[0]*x[0] + par[7]*x[1]*x[1] + par[8]*x[0]*x[1] + par[9]*x[2]*x[2]);
		aux=0.5 * par[5]*pow(r,-3) * (1./(1.+par[10]/r)-log(1.+r/par[10]));
		
		a[0]+=aux*(2*par[6]*x[0] + par[8]*x[1]);
		a[1]+=aux*(2*par[7]*x[1] + par[8]*x[0]);
		a[2]+=aux*(2*par[9]*x[2]);
		
	}else if(potential==5){
		// Spherical NFW potential
		// par = [GM, Rh]
		r=sqrt(x[0]*x[0] + x[1]*x[1]*par[2] + x[2]*x[2]*par[3]);
		aux=par[0]*pow(r,-3) * (1./(1.+par[1]/r)-log(1.+r/par[1]));
		
		a[0]=aux*x[0];
		a[1]=aux*x[1]*par[2];
		a[2]=aux*x[2]*par[3];
	}else if(potential==6){
        // Galactic potential + LMC
        // par = [GMb, ab, GMd, ad, bd^2, GM, q^2, rhalo, GMlmc, Xlmc, Ylmc, Zlmc]
        
        //Hernquist bulge
        r=sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
        aux=-par[0]/(r * (r+par[1]) * (r+par[1]));
        
        a[0]=aux*x[0];
        a[1]=aux*x[1];
        a[2]=aux*x[2];
        
        //Miyamoto disk
        aux2=sqrt(x[2]*x[2] + par[4]);
        r=sqrt(x[0]*x[0] + x[1]*x[1] + (par[3] + aux2) * (par[3] + aux2));
        aux=-par[2]/(r*r*r);
        
        a[0]+=aux*x[0];
        a[1]+=aux*x[1];
        a[2]+=aux*x[2]*(par[3] + aux2)/aux2;
        
        //Triaxial NFW Halo
        r=sqrt(par[6]*x[0]*x[0] + par[7]*x[1]*x[1] + par[8]*x[0]*x[1] + par[9]*x[2]*x[2]);
        aux=0.5 * par[5]*pow(r,-3) * (1./(1.+par[10]/r)-log(1.+r/par[10]));
        
        a[0]+=aux*(2*par[6]*x[0] + par[8]*x[1]);
        a[1]+=aux*(2*par[7]*x[1] + par[8]*x[0]);
        a[2]+=aux*(2*par[9]*x[2]);
        
        // Point mass
        // added softening ~8pc, assuming X_LMC~-0.8kpc
        r=sqrt((x[0]-par[12])*(x[0]-par[12]) + (x[1]-par[13])*(x[1]-par[13]) + (x[2]-par[14])*(x[2]-par[14])) - 0.01*par[12];
        aux = par[11]*pow(r,-3);
        
        a[0]+=aux*x[0];
        a[1]+=aux*x[1];
        a[2]+=aux*x[2];
        
    }else if(potential==7){
        // Galactic potential + LMC on a string
        // par = [GMb, ab, GMd, ad, bd^2, GM, q^2, rhalo, Mlmc, t]
        
        //Hernquist bulge
        r=sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
        aux=-par[0]/(r * (r+par[1]) * (r+par[1]));
        
        a[0]=aux*x[0];
        a[1]=aux*x[1];
        a[2]=aux*x[2];
        
        //Miyamato disk
        aux2=sqrt(x[2]*x[2] + par[4]);
        r=sqrt(x[0]*x[0] + x[1]*x[1] + (par[3] + aux2) * (par[3] + aux2));
        aux=-par[2]/(r*r*r);
        
        a[0]+=aux*x[0];
        a[1]+=aux*x[1];
        a[2]+=aux*x[2]*(par[3] + aux2)/aux2;
        
        //Triaxial NFW Halo
        r=sqrt(par[6]*x[0]*x[0] + par[7]*x[1]*x[1] + par[8]*x[0]*x[1] + par[9]*x[2]*x[2]);
        aux=0.5 * par[5]*pow(r,-3) * (1./(1.+par[10]/r)-log(1.+r/par[10]));
        
        a[0]+=aux*(2*par[6]*x[0] + par[8]*x[1]);
        a[1]+=aux*(2*par[7]*x[1] + par[8]*x[0]);
        a[2]+=aux*(2*par[9]*x[2]);
        
        //LMC on a string
        double xlmc[3]={-2.509654716638902e19, -1.2653311505262738e21, -8.319850498177284e20};
        double vlmc[3]={-57, -226, 221};
        r=sqrt((x[0]-xlmc[0]-vlmc[0]*par[12])*(x[0]-xlmc[0]-vlmc[0]*par[12]) + (x[1]-xlmc[1]-vlmc[1]*par[12])*(x[1]-xlmc[1]-vlmc[1]*par[12]) + (x[2]-xlmc[2]-vlmc[2]*par[12])*(x[2]-xlmc[2]-vlmc[2]*par[12]));
        aux=-par[11] * pow(r,-3);
        
        a[0]+=aux*x[0];
        a[1]+=aux*x[1];
        a[2]+=aux*x[2];
    }
}

void force_plummer(double *x, double *a, double Mcl)
{	// Calculate acceleration a at a position x from a cluster with a Plummer profile
	// Assumes global definitions of cluster mass Mcl and radius Rcl
	int i;
	double r, raux;
	
	r=len(x);
	raux=pow(r*r+Rcl*Rcl, 1.5);
	
	for(i=0;i<3;i++)
		a[i]=G*Mcl*x[i]/raux;
}

void initpar(int potential, double *par, double *apar)
{
	if(potential==1){
		// Logarithmic potential, par = [Vc, q_phi]
		// apar = [Vc^2, q, q^2]
		apar[0]=par[0]*par[0];
		apar[1]=par[1];
		apar[2]=par[1]*par[1];
		
	}else if(potential==2){
		// Triaxial halo potential from Law & Majewski (2010)
		// par = [Vc, phi, q_1, q_2, q_z, rhalo]
		// apar = [Vc^2, c1, c2, c3, c4, rhalo^2]
		double cosphi, sinphi;
		
		cosphi=cos(par[1]);
		sinphi=sin(par[1]);
		
		apar[0]=par[0]*par[0];
		apar[1]=cosphi*cosphi/(par[2]*par[2]) + sinphi*sinphi/(par[3]*par[3]);
		apar[2]=cosphi*cosphi/(par[3]*par[3]) + sinphi*sinphi/(par[2]*par[2]);
		apar[3]=2*sinphi*cosphi*(1/(par[2]*par[2]) - 1/(par[3]*par[3]));
		apar[4]=1/(par[4]*par[4]);
		apar[5]=par[5]*par[5];
		
	}else if(potential==3){
		// Triaxial NFW halo potential from Law & Majewski (2010)
		// par = [V, rhalo, phi, q_1, q_2, q_z]
		// apar = [GM, c1, c2, c3, c4, rhalo]
		double cosphi, sinphi;
		
		cosphi=cos(par[2]);
		sinphi=sin(par[2]);
		
// 		apar[0]=G*Msun*pow(10,par[0]);
		apar[0]=par[0]*par[0]*par[1];
		apar[1]=cosphi*cosphi/(par[3]*par[3]) + sinphi*sinphi/(par[4]*par[4]);
		apar[2]=cosphi*cosphi/(par[4]*par[4]) + sinphi*sinphi/(par[3]*par[3]);
		apar[3]=2*sinphi*cosphi*(1/(par[3]*par[3]) - 1/(par[4]*par[4]));
		apar[4]=1/(par[5]*par[5]);
		apar[5]=par[1];

	}else if(potential==0){
		// Point mass potential, par = [Mtot]
		// apar = [G*Mtot]
		apar[0]=G*par[0];
		
	}else if(potential==4){
		// Composite Galactic potential featuring a disk, bulge, and triaxial NFW halo (from Johnston/Law/Majewski/Helmi)
		// par = [GMb, ab, GMd, ad, bd, V, rhalo, phi, q_1, q_2, q_z]
		// apar = [GMb, ab, GMd, ad, bd^2, GM, c1, c2, c3, c4, rhalo]
		double cosphi, sinphi; //, tq, tphi;
		
		apar[0]=G*par[0];
		apar[1]=par[1];
		apar[2]=G*par[2];
		apar[3]=par[3];
		apar[4]=par[4]*par[4];
		
		cosphi=cos(par[7]);
		sinphi=sin(par[7]);
		
		apar[5]=par[5]*par[5]*par[6];
		apar[6]=cosphi*cosphi/(par[8]*par[8]) + sinphi*sinphi/(par[9]*par[9]);
		apar[7]=cosphi*cosphi/(par[9]*par[9]) + sinphi*sinphi/(par[8]*par[8]);
		apar[8]=2*sinphi*cosphi*(1/(par[8]*par[8]) - 1/(par[9]*par[9]));
		apar[9]=1/(par[10]*par[10]);
		apar[10]=par[6];

    }else if(potential==5){
		apar[0]=G*Msun*pow(10,par[0]);
		apar[1]=par[1];
		apar[2]=par[5]*par[5];
		apar[3]=par[6]*par[6];

	}else if(potential==6){
        // Galactic potential + LMC
        // par = [GMb, ab, GMd, ad, bd, Vh, rhalo, phi, q1, q2, qz, Mlmc, Xlmc, Ylmc, Zlmc]
        // apar = [GMb, ab, GMd, ad, bd^2, GM, c1, c2, c3, c4, rhalo, GMlmc, Xlmc, Ylmc, Zlmc]
        double cosphi, sinphi;

        apar[0]=G*par[0];
        apar[1]=par[1];
        
        apar[2]=G*par[2];
        apar[3]=par[3];
        apar[4]=par[4]*par[4];
        
        cosphi=cos(par[7]);
        sinphi=sin(par[7]);
        apar[5]=par[5]*par[5]*par[6];
        apar[6]=cosphi*cosphi/(par[8]*par[8]) + sinphi*sinphi/(par[9]*par[9]);
        apar[7]=cosphi*cosphi/(par[9]*par[9]) + sinphi*sinphi/(par[8]*par[8]);
        apar[8]=2*sinphi*cosphi*(1/(par[8]*par[8]) - 1/(par[9]*par[9]));
        apar[9]=1/(par[10]*par[10]);
        apar[10]=par[6];
        
        apar[11]=G*par[11];
        apar[12]=par[12];
        apar[13]=par[13];
        apar[14]=par[14];
    }else if(potential==7){
        // Galactic potential + LMC on a string
        // par = [GMb, ab, GMd, ad, bd, V, rhalo, phi, q_1, q_2, q_z, Mlmc]
        // apar = [GMb, ab, GMd, ad, bd^2, GM, c1, c2, c3, c4, rhalo, Mlmc, t]
        double cosphi, sinphi; //, tq, tphi;
        
        apar[0]=G*par[0];
        apar[1]=par[1];
        apar[2]=G*par[2];
        apar[3]=par[3];
        apar[4]=par[4]*par[4];
        
        cosphi=cos(par[7]);
        sinphi=sin(par[7]);
        apar[5]=par[5]*par[5]*par[6];
        apar[6]=cosphi*cosphi/(par[8]*par[8]) + sinphi*sinphi/(par[9]*par[9]);
        apar[7]=cosphi*cosphi/(par[9]*par[9]) + sinphi*sinphi/(par[8]*par[8]);
        apar[8]=2*sinphi*cosphi*(1/(par[8]*par[8]) - 1/(par[9]*par[9]));
        apar[9]=1/(par[10]*par[10]);
        apar[10]=par[6];
        apar[11]=G*par[11];
        printf("%e\n", par[11]);
        apar[12]=0.;
    }
}

double jacobi(double *x, double *v, double *par, int potential, double Mcl)
{	// Jacobi radius of a cluster
	// at the position x, velocity v, and in a given potential
	int i;
	double R, om, dpot, delta, r;
	double omega[3], x1[3], x2[3], a1[3], a2[3], dx[3]; //, da[3];

	// Radial distance
	r=len(x);
	
	// Angular velocity
	omega[0]=x[1]*v[2]-x[2]*v[1];
	omega[1]=x[2]*v[0]-x[0]*v[2];
	omega[2]=x[0]*v[1]-x[1]*v[0];
	om=len(omega)/(r*r);
	
	// Potential derivative
	delta=0.02*kpc;
	for(i=0;i<3;i++){
		x1[i]=x[i]/r*(r-delta);
		x2[i]=x[i]/r*(r+delta);
	}
	force(x1, a1, par, potential);
	force(x2, a2, par, potential);
	for(i=0;i<3;i++){
		dx[i]=x1[i]-x2[i];
	}
	dpot=(len(a1)-len(a2))/len(dx);
	
	// Jacobi radius
	R=pow(G*Mcl/fabs(om*om+dpot),1./3.);
	
	return R;
}
