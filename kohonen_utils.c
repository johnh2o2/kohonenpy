// C functions to speed up Kohonen module
// NOT READY TO BE COMPILED!!

#include<stdio.h>
#include<stdlib.h>
#include<Python.h>

static PyMethodDef ck_methods[] = {
	{"find_bmu", ck_find_bmu, METH_VARARGS,
		"Finds the indices of the best matched unit in the grid"},
	{"imprint_vector_on_grid", ck_imprint_vector_on_grid, METH_VARARGS,
		"Imprints the vector onto the grid (used during training)"},
	{NULL, NULL, 0, NULL}        /* Sentinel */
};

double
ck_euclidean(double *a, double *b, int N){
	int i;
	double result = 0;
	for(i=0; i<N; i++) result += pow(a[i]-b[i],2)
	return sqrt(result)
}

double *
ck_phase_shift(double *a, int p, int N){
	double b[N];
	int i;
	for (i=0; i < N; i++) b[i] = a[(i+p)%N];
	return b
}

double
ck_min_euclidean(double *a, double *b, int N){
	int i, j;
	double min_dist = -1;
	double new_dist;
	for(i=0; i < N; i++){
		new_dist = ck_euclidean(a, ck_phase_shift(b, i, N));
		if (min_dist < 0 || new_dist < min_dist) min_dist = new_dist;
	}
	return min_dist;
}

void
ck_convert_grid_and_vector(PyObject *grid, PyObject *vec, 
							double ***cgrid, double *cvec, 
							npy_intp *grid_dims, npy_intp *vec_dims){
	int typenum = NPY_DOUBLE;
	PyArray_Descr *descr;
	descr = PyArray_DescrFromType(typenum);
	
	if (PyArray_AsCArray(&grid, (void ***)&cgrid, grid_dims, 3, descr) < 0){
		PyErr_SetString(PyExc_TypeError, "error converting grid to c array");
		return NULL;
	}
	if (PyArray_AsCArray(&vec, (void *)&cvec, vec_dims, 1, descr) < 0){
		PyErr_SetString(PyExc_TypeError, "error converting vector to c array");
		return NULL;
	}
}

/*
static PyObject *heirarchical_grid = NULL;

static PyObject *
ck_set_heirarchical_grid(PyObject *self, PyObject *args){

	
}
*/

//Find best-match unit (slow)
static PyObject *
ck_find_bmu(PyObject *self, PyObject *args){
	PyObject *pygrid, *pyvec;
	const char *dist_func_type;

	if (!PyArg_ParseTuple(args, "OOs", &pygrid, &pyvec, &dist_func_type)) return NULL;

	// Choose metric function
	if (dist_func_type == 'euc'){
		double (*dist)(double *a, double *b, int N) = ck_euclidean;
	}
	else if (dist_func_type == 'min_euc') {
		double (*dist)(double *a, double *b, int N) = ck_min_euclidean;
	}
	else {
		PyErr_SetString(PyExc_TypeError, "Do not understand dist_func_type");
		return NULL;
	}

	double ***cgrid;
	double *cvec;
	npy_intp griddims[3],vecdims[1];
	int i,j,bmu[2];
	double min_dist = -1, new_dist;

	//Create C arrays from numpy objects:
	ck_convert_grid_and_vector(pygrid, pyvec,  cgrid, cvec, griddims, vecdims);

	//Find the grid element with the smallest distance
	for (i=0; i < griddims[0]; i ++ ){
		for (j=0; j < griddims[1]; j ++){
			new_dist = dist(cgrid[i][j], cvec, vecdims[0])
			if (min_dist < 0 || new_dist < min_dist) {
				min_dist = new_dist;
				bmu = { i, j }
			}
		}
	}

	return Py_BuildValue("ii", bmu[0], bmu[1])
		
}	

static PyObject *
ck_imprint_vector_on_grid(PyPbject *self, PyObject *args){

	PyObject *pygrid, *pyvec, *Theta, *Alpha, *thargs;

	if (!PyArg_ParseTuple(args, "OOOO", &pygrid, &pyvec, &Theta, &Alpha)) return NULL;

	double ***cgrid;
	double *cvec;
	double Wvs, theta, alpha;
	alpha = PyFloat_AsDouble(Alpha)
	int i,j,k;
	npy_intp griddims[3],vecdims[1];

	//Create C arrays from numpy objects:
	ck_convert_grid_and_vector(pygrid, pyvec, cgrid, cvec, griddims, vecdims);

	for (i=0; i < griddims[0]; i++ ){
		for (j=0; j < griddims[1]; j++){
			Wvs = cgrid[i][j]
			thargs = Py_BuildValue("(ii)", i,j)
			theta = PyEval_CallObject(Theta, thargs);
			for (k=0; k < vecdims[0]; k++){
				cgrid[i][j][k] = Wvs[k] + theta * alpha * (D[k] - Wvs[k])
			}
		}
	}
	//TODO: figure out how to change the object's grid...?
	return NULL;
}	

PyMODINIT_FUNC
initkohonen_utils(void)
{
    (void) Py_InitModule("ckutils", ck_methods);
}

int
main(int argc, char *argv[])
{
    /* Pass argv[0] to the Python interpreter */
    Py_SetProgramName(argv[0]);

    /* Initialize the Python interpreter.  Required. */
    Py_Initialize();

    /* Add a static module */
    initkohonen_utils();
}
