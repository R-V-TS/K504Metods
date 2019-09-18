%module K504Metods

%{
#define SWIG_FILE_WITH_INIT
#include "src/DCT.h"
#include "src/RAWImage.h"
%}

%include "numpy.i" 

%init %{
	import_array();
%}

%apply (short* IN_ARRAY3, int DIM1, int DIM2, int DIM3) {(short* pixel_array, int image_width, int image_height, int channels_)}

%apply (float* ARGOUT_ARRAY1, int DIM1) {(float* DCT_Coff, int n)}
%apply (short* ARGOUT_ARRAY1, int DIM1) {(short* pixel_array, int n)}

%include src/RAWImage.h
%include src/DCT.h