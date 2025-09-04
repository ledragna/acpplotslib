#https://cython.readthedocs.io/en/latest/src/userguide/numpy_tutorial.html#numpy-tutorial
# distutils: extra_compile_args=-fopenmp
# distutils: extra_link_args=-fopenmp
#cython: language_level=3

import numpy as np
cimport cython
from cython.parallel import prange, parallel
cimport openmp
from libc.math cimport sqrt, log, exp, pi, M_LN2

@cython.boundscheck(False)
@cython.wraparound(False)
def gaussian(double[:] array_x, double hwhm, double xzero, double yzero):
    cdef Py_ssize_t xlen = array_x.shape[0]
    array_y = np.zeros(xlen, dtype=np.double)
    cdef double[::1] array_y_view = array_y
    assert hwhm > 1e-5
    cdef double sigma2, norm, tmp_pow
    sigma2 = hwhm*hwhm/M_LN2
    norm = 1./sqrt(sigma2*pi)
    for tmp_x in range(xlen):
        tmp_pow = (array_x[tmp_x]-xzero)*(array_x[tmp_x]-xzero)
        array_y_view[tmp_x] = yzero*norm*exp(-tmp_pow/sigma2)
    return array_y

@cython.boundscheck(False)
@cython.wraparound(False)
def lorentzian(double[:] array_x, double hwhm, double xzero, double yzero):
    cdef Py_ssize_t xlen = array_x.shape[0]
    array_y = np.zeros(xlen, dtype=np.double)
    cdef double[::1] array_y_view = array_y
    assert hwhm > 1e-5
    cdef double norm, tmp_pow
    norm = hwhm/pi  # Normalization factor
    for tmp_x in range(xlen):
        tmp_pow = (array_x[tmp_x]-xzero)*(array_x[tmp_x]-xzero) + hwhm * hwhm
        array_y_view[tmp_x] = yzero*norm/tmp_pow
    return array_y

cdef inline double gauss_help(double ith_x, double x_zero, double y_zero,
                              double sigma2, double norm):
    """TODO"""
    return y_zero * norm * exp(-((ith_x - x_zero) * (ith_x - x_zero)) /sigma2)

cdef inline double loren_help(double ith_x, double x_zero, double y_zero,
                              double hwhm2, double norm):
    """Todo"""
    return y_zero * norm / (ith_x - x_zero) * (ith_x - x_zero) * hwhm2

@cython.boundscheck(False)
@cython.wraparound(False)
def broadgau(double[:] array_x_trans, double[:] array_y_trans, double [:] xaxis,
             double hwhm, int num_threads): # yfac=1.0, mul_x=False):
    """Broadens transitions stored in two separate arrays, x_trans, y_trans."""
    cdef:
        Py_ssize_t ymax = array_y_trans.shape[0]
        Py_ssize_t xlen = xaxis.shape[0]
        # pad local data to 64 byte avoid false sharing of cache-lines
        Py_ssize_t xaxis_padded = (((xlen - 1) // 8) + 1) * 8
        Py_ssize_t tmp_y, tmp_x, ind
        double[:] ylocal = np.zeros(xaxis_padded * num_threads, dtype=np.double)
        int tid
        double sigma2 = hwhm*hwhm/M_LN2
        double norm = 1./sqrt(sigma2*pi)

    assert hwhm > 1e-5

    yaxis = np.zeros(xlen, dtype=np.double)
    cdef double[::1] yaxis_view = yaxis

    with nogil, parallel(num_threads=num_threads):
        tid = openmp.omp_get_thread_num()
        for tmp_y in prange(ymax, schedule='static', chunksize=1):
            for tmp_x in range(xlen):
               ylocal[tid * xaxis_padded + tmp_x] += array_y_trans[tmp_y] * norm *\
                                                     exp(-((xaxis[tmp_x] - array_x_trans[tmp_y]) *
                                                           (xaxis[tmp_x] - array_x_trans[tmp_y])) / sigma2)
    for tid in range(num_threads):
        for ind in range(xlen):
            yaxis_view[ind] += ylocal[tid * xaxis_padded + ind]

    return yaxis


@cython.boundscheck(False)
@cython.wraparound(False)
def broadlor(double[:] array_x_trans, double[:] array_y_trans, double [:] xaxis,
             double hwhm, int num_threads): # yfac=1.0, mul_x=False):
    """Broadens transitions stored in two separate arrays, x_trans, y_trans."""
    cdef:
        Py_ssize_t ymax = array_y_trans.shape[0]
        Py_ssize_t xlen = xaxis.shape[0]
        # pad local data to 64 byte avoid false sharing of cache-lines
        Py_ssize_t xaxis_padded = (((xlen - 1) // 8) + 1) * 8
        Py_ssize_t tmp_y, tmp_x, ind
        double[:] ylocal = np.zeros(xaxis_padded * num_threads, dtype=np.double)
        int tid
        double hwhm2 = hwhm * hwhm
        double norm = hwhm/pi

    assert hwhm > 1e-5

    yaxis = np.zeros(xlen, dtype=np.double)
    cdef double[::1] yaxis_view = yaxis

    with nogil, parallel(num_threads=num_threads):
        tid = openmp.omp_get_thread_num()
        for tmp_y in prange(ymax, schedule='static', chunksize=1):
            for tmp_x in range(xlen):
               ylocal[tid * xaxis_padded + tmp_x] += array_y_trans[tmp_y] * norm /\
                                                    ((xaxis[tmp_x] - array_x_trans[tmp_y]) *
                                                     (xaxis[tmp_x] - array_x_trans[tmp_y]) + hwhm2)
    for tid in range(num_threads):
        for ind in range(xlen):
            yaxis_view[ind] += ylocal[tid * xaxis_padded + ind]

    return yaxis


@cython.boundscheck(False)
@cython.wraparound(False)
def broadening2(double[:] array_x_trans, double[:] array_y_trans, double [:] xaxis, double hwhm): # yfac=1.0, mul_x=False):
    """Broadens transitions stored in two separate arrays, x_trans, y_trans."""
    cdef:
        Py_ssize_t ymax = array_y_trans.shape[0]
        Py_ssize_t xlen = xaxis.shape[0]
        Py_ssize_t tmp_y, tmp_x, ind
        double sigma2 = hwhm*hwhm/M_LN2
        double norm = 1./sqrt(sigma2*pi)

    assert hwhm > 1e-5
    yaxis = np.zeros(xlen, dtype=np.double)
    cdef double[::1] yaxis_view = yaxis
    for tmp_y in range(ymax):
        for tmp_x in range(xlen):
           yaxis_view[tmp_x] += array_y_trans[tmp_y] * norm *\
                                exp(-((xaxis[tmp_x] - array_x_trans[tmp_y]) *
                                      (xaxis[tmp_x] - array_x_trans[tmp_y])) / sigma2)
    return yaxis



