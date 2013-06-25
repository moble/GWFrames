/**************************************************************************

    Copyright 2010-2012  Kevin M. Huffenberger & Benjamin D. Wandelt

    This file is part of spinsfast.

    Spinsfast is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Spinsfast is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with spinsfast.  If not, see <http://www.gnu.org/licenses/>.

***************************************************************************/

/* Code revision: 104, 2012-04-13 13:00:16 -0400 (Fri, 13 Apr 2012) */

#include <Python.h>
#include <numpy/ndarrayobject.h>
#include <fftw3.h>
#include <alm.h>
#include <wigner_d_halfpi.h>
#include <spinsfast_forward.h>
#include <spinsfast_backward.h>


///  Indexing

static PyObject *spinsfast_N_lm(PyObject *self, PyObject *args) {

  int lmax;
  if (!PyArg_ParseTuple(args, "i", &lmax)) 
    return(NULL);
    
  return Py_BuildValue("i", N_lm(lmax));
}


static PyObject *spinsfast_lm_ind(PyObject *self, PyObject *args) {
  int l,m, lmax;
  
  if (!PyArg_ParseTuple(args, "iii", &l,&m, &lmax)) 
    return NULL;
 
  return Py_BuildValue("i", lm_ind(l,m,lmax));

}


static PyObject *spinsfast_ind_lm(PyObject *self, PyObject *args) {

  int i, lmax;

  if (!PyArg_ParseTuple(args, "ii", &i, &lmax)) 
    return NULL;
  
  
  int nd = 1;
  npy_intp dims[1] = {2};
  
  int *lm = calloc(2,sizeof(int));
  
  ind_lm(i, &lm[0], &lm[1], lmax);
  
  PyObject *arr = PyArray_SimpleNewFromData(nd, dims, NPY_INT, lm);
  Py_INCREF(arr);
  
  return(arr);
}





///  Transform wrappers  ////////////////


static PyObject *spinsfast_mod_salm2map(PyObject *self, PyObject *args) {

  PyObject *input_array=NULL;
  int lmax = 0;
  int s = 0;
  int Ntheta = 0;
  int Nphi = 0;

  if (!PyArg_ParseTuple(args, "Oiiii", &input_array, &s, &lmax, &Ntheta, &Nphi)) 
    return NULL;

  fftw_complex *alm = PyArray_DATA(input_array);
  fftw_complex *f = calloc(Nphi*Ntheta,sizeof(fftw_complex));
  
  spinsfast_salm2map(alm, f, s, Ntheta, Nphi, lmax);
  
  npy_intp dims[2] = {Ntheta, Nphi};

  PyObject *arr = PyArray_SimpleNewFromData(2, dims, NPY_CDOUBLE, f);
  Py_INCREF(arr);
  
  return(arr);

}

static PyObject *spinsfast_mod_map2salm(PyObject *self, PyObject *args) {

  PyObject *input_array=NULL;
  int lmax = 0;
  int s = 0;

  if (!PyArg_ParseTuple(args, "Oii", &input_array, &s, &lmax)) 
    return NULL;


  npy_intp *dim = PyArray_DIMS(input_array);
  int Ntheta = dim[0];
  int Nphi = dim[1];


  npy_intp Nlm = N_lm(lmax);

  
  fftw_complex *alm = calloc(Nlm,sizeof(fftw_complex));
  fftw_complex *f = PyArray_DATA(input_array);

  spinsfast_map2salm(f, alm, s, Ntheta, Nphi, lmax);

  PyObject *arr = PyArray_SimpleNewFromData(1, &Nlm, NPY_CDOUBLE, alm);
  Py_INCREF(arr);
  
  return(arr);

}


///  Some helpers
static PyObject *spinsfast_quad_weights(PyObject *self, PyObject *args) {
  int Ntheta;

  if (!PyArg_ParseTuple(args, "i", &Ntheta)) 
    return NULL;

  int wsize = 2*(Ntheta-1);
  fftw_complex *W = calloc(wsize, sizeof(fftw_complex)); // fourier space weights

  spinsfast_quadrature_weights(W, wsize);
  
  npy_intp N = wsize;

  PyObject *arr = PyArray_SimpleNewFromData(1, &N, NPY_CDOUBLE, W);
  Py_INCREF(arr);
  
  return(arr);
}

static PyObject *spinsfast_mod_f_extend_MW(PyObject *self, PyObject *args) {

  PyObject *input_array=NULL;
  int s = 0;

  if (!PyArg_ParseTuple(args, "Oi", &input_array, &s)) 
    return NULL;

  
  fftw_complex *f = PyArray_DATA(input_array);
  npy_intp *dim = PyArray_DIMS(input_array);
  // int contig = PyArray_ISCONTIGUOUS(input_array);
  int Ntheta = dim[0];
  int Nphi = dim[1];

  int wsize = 2*(Ntheta-1);
  fftw_complex *F = fftw_malloc(wsize*Nphi*sizeof(fftw_complex));
  
  spinsfast_f_extend_MW(f, F, s, Ntheta, Nphi);

  npy_intp N[] = {wsize,Nphi};
  
  PyObject *arr = PyArray_SimpleNewFromData(2, N, NPY_CDOUBLE, F);
  Py_INCREF(arr);
  
  return(arr);

}
static PyObject *spinsfast_mod_f_extend_old(PyObject *self, PyObject *args) {

  PyObject *input_array=NULL;
  int s = 0;

  if (!PyArg_ParseTuple(args, "Oi", &input_array, &s)) 
    return NULL;

  
  fftw_complex *f = PyArray_DATA(input_array);
  npy_intp *dim = PyArray_DIMS(input_array);
  // int contig = PyArray_ISCONTIGUOUS(input_array);
  int Ntheta = dim[0];
  int Nphi = dim[1];

  int wsize = 2*(Ntheta-1);
  fftw_complex *F = fftw_malloc(wsize*Nphi*sizeof(fftw_complex));
  
  spinsfast_f_extend_old(f, F, s, Ntheta, Nphi);

  npy_intp N[] = {wsize,Nphi};
  
  PyObject *arr = PyArray_SimpleNewFromData(2, N, NPY_CDOUBLE, F);
  Py_INCREF(arr);
  
  return(arr);

}

static PyObject *spinsfast_mod_Imm(PyObject *self, PyObject *args) {

  PyObject *input_array=NULL;
  int lmax = 0;
  int s = 0;
  
  if (!PyArg_ParseTuple(args, "Oii", &input_array, &s, &lmax)) 
    return NULL;
  
  
  fftw_complex *f = PyArray_DATA(input_array);
  npy_intp *dim = PyArray_DIMS(input_array);
  // int contig = PyArray_ISCONTIGUOUS(input_array);
  int Ntheta = dim[0];
  int Nphi = dim[1];
  
  
  // int Npix = Nphi * Ntheta;
  int Nm = 2*lmax+1;
  fftw_complex *Imm = fftw_malloc(Nm*Nm*sizeof(fftw_complex));
  
  
  spinsfast_forward_multi_Imm (f, &s, 1, Ntheta, Nphi, lmax, Imm);
  
  npy_intp N[] = {Nm,Nm};
  
  PyObject *arr = PyArray_SimpleNewFromData(2, N, NPY_CDOUBLE, Imm);
  Py_INCREF(arr);
  
  return(arr);

}



/////////////////// Module info /////////////////////////////


static PyMethodDef spinsfastMethods[] = {
  {"N_lm",  spinsfast_N_lm, METH_VARARGS,"N_lm(lmax)"},
  {"ind_lm",  spinsfast_ind_lm, METH_VARARGS,"ind_lm(idx,lmax)"},
  {"lm_ind",  spinsfast_lm_ind, METH_VARARGS,"lm_ind(l,m,lmax)"},
  {"salm2map",  spinsfast_mod_salm2map, METH_VARARGS,"salm2map(alm,s,lmax,Ntheta,Nphi)"},
  {"map2salm",  spinsfast_mod_map2salm, METH_VARARGS,"map2salm(f,s,lmax)"},
  {"f_extend_MW",  spinsfast_mod_f_extend_MW, METH_VARARGS,"f_extend_MW(f,s)"},
  {"f_extend_old",  spinsfast_mod_f_extend_old, METH_VARARGS,"f_extend_MW(f,s)"},
  {"Imm",  spinsfast_mod_Imm, METH_VARARGS,"Imm(f,s,lmax)"},
  {"quadrature_weights",  spinsfast_quad_weights, METH_VARARGS,"quadrature_weights(Ntheta)"},
  {NULL, NULL, 0, NULL}        /* Sentinel */
};


PyMODINIT_FUNC initspinsfast(void)
{
  (void) Py_InitModule("spinsfast", spinsfastMethods);
  import_array();  // This is important for using the numpy_array api, otherwise segfaults!
}
