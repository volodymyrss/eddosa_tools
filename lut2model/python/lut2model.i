%module lut2model

%{
    #define SWIG_FILE_WITH_INIT
    #include "lut2model.h"
%}

%include "numpy.i"

%init %{
    import_array();
%}

%apply (double* ARGIN_ARRAY1, int DIM1) {(double* inarr, int m)}
%apply (double* ARGOUT_ARRAY1, int DIM1) {(double* outarr, int n)}

%apply (double* IN_ARRAY2, int DIM1, int DIM2) {(double *inarr, int in1, int in2)}
%apply (double* INPLACE_ARRAY2, int DIM1, int DIM2) {(double *inplarr, int inpl1, int inpl2)}

%include "lut2model.h"
