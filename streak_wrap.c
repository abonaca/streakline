#include <Python.h>
#include <numpy/arrayobject.h>
#include "streakline.h"

/* Docstrings */
static char module_docstring[] = "This module provides an interface for running a leapfrog orbit integrator using C.";
static char stream_docstring[] = "Calculate orbit in a point mass potential given the inital position and velocity.";
    
/* Available functions */
static PyObject *streakline_stream(PyObject *self, PyObject *args);

/* Module specification */
static PyMethodDef module_methods[] = {
	{"stream", streakline_stream, METH_VARARGS, stream_docstring},
	{NULL, NULL, 0, NULL}
};

/* Initialize the module */
PyMODINIT_FUNC initstreakline(void)
{
	PyObject *m = Py_InitModule3("streakline", module_methods, module_docstring);
	if (m == NULL)
		return;

	/* Load `numpy` functionality. */
	import_array();
}

double *pyvector_to_Carrayptrs(PyArrayObject *arrayin);

static PyObject *streakline_stream(PyObject *self, PyObject *args)
{
	int N, M, Ne, err, potential, integrator;
	double *x0, *v0, *par, *offset, mcli, mclf, rcl, dt_;
	PyObject *par_obj, *par_array, *x0_obj, *x0_array, *v0_obj, *v0_array, *offset_obj, *offset_array;

	// Parse the input tuple
	if (!PyArg_ParseTuple(args, "OOOOiiiidddd", &x0_obj, &v0_obj, &par_obj, &offset_obj, &potential, &integrator, &N, &M, &mcli, &mclf, &rcl, &dt_))	// reads in input parameters: initial position and velocity, choice of potential, potential parameters and the number of timesteps N
		return NULL;
	Ne=ceil((float)N/(float)M);
	
	// Interpret the input parameters as numpy arrays
	par_array = PyArray_FROM_OTF(par_obj, NPY_DOUBLE, NPY_IN_ARRAY);
	offset_array = PyArray_FROM_OTF(offset_obj, NPY_DOUBLE, NPY_IN_ARRAY);
	x0_array = PyArray_FROM_OTF(x0_obj, NPY_DOUBLE, NPY_IN_ARRAY);
	v0_array = PyArray_FROM_OTF(v0_obj, NPY_DOUBLE, NPY_IN_ARRAY);
	
	//If that didn't work, throw an exception
	if (par_array == NULL) {
		Py_XDECREF(par_array);
		return NULL;
	}
	if (offset_array == NULL) {
		Py_XDECREF(offset_array);
		return NULL;
	}
	if (x0_array == NULL) {
		Py_XDECREF(x0_array);
		return NULL;
	}
	if (v0_array == NULL) {
		Py_XDECREF(v0_array);
		return NULL;
	}
	// How many parameters are there?
// 	int Npar = (int)PyArray_DIM(par_array, 0);

	//Get pointers to the data as C-types. */
	par = (double*)PyArray_DATA(par_array);
	offset = (double*)PyArray_DATA(offset_array);
	x0 = (double*)PyArray_DATA(x0_array);
	v0 = (double*)PyArray_DATA(v0_array);
	
	// Set up return array pointers
	double *xm1, *xm2, *xm3, *xp1, *xp2, *xp3, *vm1, *vm2, *vm3, *vp1, *vp2, *vp3;
	int nd=1, dims[2];
	dims[0] = Ne;
	PyArrayObject *py_xm1, *py_xm2, *py_xm3, *py_xp1, *py_xp2, *py_xp3, *py_vm1, *py_vm2, *py_vm3, *py_vp1, *py_vp2, *py_vp3;
	
	// Python arrays
	py_xm1 = (PyArrayObject*) PyArray_FromDims(nd, dims, NPY_DOUBLE);
	py_xm2 = (PyArrayObject*) PyArray_FromDims(nd, dims, NPY_DOUBLE);
	py_xm3 = (PyArrayObject*) PyArray_FromDims(nd, dims, NPY_DOUBLE);
	py_xp1 = (PyArrayObject*) PyArray_FromDims(nd, dims, NPY_DOUBLE);
	py_xp2 = (PyArrayObject*) PyArray_FromDims(nd, dims, NPY_DOUBLE);
	py_xp3 = (PyArrayObject*) PyArray_FromDims(nd, dims, NPY_DOUBLE);
	py_vm1 = (PyArrayObject*) PyArray_FromDims(nd, dims, NPY_DOUBLE);
	py_vm2 = (PyArrayObject*) PyArray_FromDims(nd, dims, NPY_DOUBLE);
	py_vm3 = (PyArrayObject*) PyArray_FromDims(nd, dims, NPY_DOUBLE);
	py_vp1 = (PyArrayObject*) PyArray_FromDims(nd, dims, NPY_DOUBLE);
	py_vp2 = (PyArrayObject*) PyArray_FromDims(nd, dims, NPY_DOUBLE);
	py_vp3 = (PyArrayObject*) PyArray_FromDims(nd, dims, NPY_DOUBLE);
	
	// Pointers to C arrays
	xm1 = pyvector_to_Carrayptrs(py_xm1);
	xm2 = pyvector_to_Carrayptrs(py_xm2);
	xm3 = pyvector_to_Carrayptrs(py_xm3);
	xp1 = pyvector_to_Carrayptrs(py_xp1);
	xp2 = pyvector_to_Carrayptrs(py_xp2);
	xp3 = pyvector_to_Carrayptrs(py_xp3);
	vm1 = pyvector_to_Carrayptrs(py_vm1);
	vm2 = pyvector_to_Carrayptrs(py_vm2);
	vm3 = pyvector_to_Carrayptrs(py_vm3);
	vp1 = pyvector_to_Carrayptrs(py_vp1);
	vp2 = pyvector_to_Carrayptrs(py_vp2);
	vp3 = pyvector_to_Carrayptrs(py_vp3);

	// Call the external C function to calculate the geostationary orbit.
	err = stream(x0, v0, xm1, xm2, xm3, xp1, xp2, xp3, vm1, vm2, vm3, vp1, vp2, vp3, par, offset, potential, integrator, N, M, mcli, mclf, rcl, dt_);

	// Check if error raised
	if(err!=0) {
		PyErr_SetString(PyExc_RuntimeError, "Error occured in the leapfrog integrator.");
		return NULL;
	}
	
	// Store return array
	PyObject *out = Py_BuildValue("OOOOOOOOOOOO", py_xm1, py_xm2, py_xm3, py_xp1, py_xp2, py_xp3, py_vm1, py_vm2, py_vm3, py_vp1, py_vp2, py_vp3);
	
	// Clean up
	Py_XDECREF(par_array);
	Py_XDECREF(py_xm1);
	Py_XDECREF(py_xm2);
	Py_XDECREF(py_xm3);
	Py_XDECREF(py_xp1);
	Py_XDECREF(py_xp2);
	Py_XDECREF(py_xp3);
	Py_XDECREF(py_vm1);
	Py_XDECREF(py_vm2);
	Py_XDECREF(py_vm3);
	Py_XDECREF(py_vp1);
	Py_XDECREF(py_vp2);
	Py_XDECREF(py_vp3);
	
	// Return positions, velocities and energy as a function of time
	return out;
}

double *pyvector_to_Carrayptrs(PyArrayObject *arrayin)
{
// 	int n=arrayin->dimensions[0];
	return (double *) arrayin->data;  /* pointer to arrayin data as double */
}
