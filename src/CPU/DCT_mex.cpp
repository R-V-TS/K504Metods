/*==========================================================
 * DCT_mex.c - simple DCT filter
 *
 * Filtered image
 *
 * The calling syntax is:
 *
 *      filtered = DCT_mex(channels, window_size, threshold, image_pixels)
 *
 * Please, use reshape(filtered, cols, rows, channels) after useng this function
 *
 * This is a MEX-file for MATLAB.
 * Copyright 2019, Rostyslav Tsekhmystro.
 *
 *========================================================*/

#include "mex.h"
#include "DCT.cpp"

void DCT_mex_filter(int cols, int rows, int channels, double* image_pixel, int window_size, double threshold, double* filtered_image)
{
    double *temp_array = new double[cols*rows];
    for(int z = 0; z < channels; z++)
    {
        for(int i = 0; i < cols; i++)
        {    
            for(int j = 0; j < rows; j++)
            {
                temp_array[(rows*i)+j] = image_pixel[(z*cols*rows)+(rows*i)+j];
            }
        }
        ImProcessing::DCT_FILT filter(temp_array, cols, rows, window_size, threshold); // JPEG = 0.012
        double *pixel_result = filter.imageFilter();
        for(int i = 0; i < rows; i++)
        {    
            for(int j = 0; j < cols; j++)
            {
                mexPrintf("%f ", pixel_result[(cols*i)+j]);
                filtered_image[(z*cols*rows)+(cols*i)+j] = pixel_result[(cols*i)+j];
            }
        }
    }
}

void mexFunction(int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    /*chech for 6 arguments
    if(nrhs < 6){
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Six inputs required.");
    }
    
    Check arguments type
    if(!mxIsUint8(prhs[0])){
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notScalar", "Input 1 must be a uint8 type.");
    }
    if(!mxIsUint8(prhs[1])){
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notScalar", "Input 2 must be a uint8 type.");
    }
    if(!mxIsUint8(prhs[2])){
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notScalar", "Input 3 must be a uint8 type.");
    }
    if(!mxIsUint8(prhs[3])){
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notScalar", "Input 4 must be a uint8 type.");
    }
    if(!mxIsDouble(prhs[4])){
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notScalar", "Input 5 must be a double type.");
    }
    if(!mxIsDouble(prhs[5])){
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notScalar", "Input 6 must be a double type.");
    }*/
       
    /*int cols = *(int*)mxGetData(prhs[0]);
    int rows = *(int*)mxGetData(prhs[1]);
    int channels = *(int*)mxGetData(prhs[2]);
    int window = *(int*)mxGetData(prhs[3]);
    double threshold = mxGetScalar(prhs[4]);
    double* pixels = mxGetDoubles(prhs[5]);
    */
    
    int channels = *(int*)mxGetData(prhs[0]);
    int window = *(int*)mxGetData(prhs[1]);
    double threshold = mxGetScalar(prhs[2]);
    double* pixels = mxGetDoubles(prhs[3]);
    int cols = mxGetN(prhs[3])/channels;
    int rows = mxGetM(prhs[3]);
    mexPrintf("%i %i\n", cols, rows);

    
    plhs[0] = mxCreateDoubleMatrix((mwSize)rows, (mwSize)(cols*channels), mxREAL);
    double* outMatrix = mxGetDoubles(plhs[0]);
    
    DCT_mex_filter(cols, rows, channels, pixels, window, threshold, outMatrix);
}







