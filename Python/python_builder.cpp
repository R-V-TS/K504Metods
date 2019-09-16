//
// Created by rostislav on 24.08.19.
//
// Do not include to CMake build target
#import <Python.h>
#import <iostream>
#include "../src/CPU/DCT.cpp"

double* image;
int width;
int height;
int channels;

extern "C" void DCT_frequency(short* image_array, int width, int height, int channels)
{
    int window_size = 8;
    int matrix_flag[8][8] = {
            {0,  1,  1,  1,  1,  2,  2,  2},
            {1,  1,  1,  1,  2,  2,  2,  3},
            {1,  1,  1,  2,  2,  2,  3,  3},
            {1,  1,  2,  2,  2,  3,  3,  3},
            {1,  2,  2,  2,  3,  3,  3,  4},
            {2,  2,  2,  3,  3,  3,  4,  4},
            {2,  2,  3,  3,  3,  4,  4,  4},
            {2,  3,  3,  3,  4,  4,  4,  4}
    };

    float* EST = new float[(width/8)*(height/8)*4];
    float* coefficients = new float[16];
    float* block = new float[window_size*window_size];
    using namespace ImProcessing;
    char *image_pixel = new char[width*height];
    for(int i = 0, j = 0; j <width*height*channels; i++, j+=3) image_pixel[i] = (char) image_array[j];

    float *DCT_ARRAY = DCT(image_pixel, width, height, window_size);
    int z = 0;

    float slices_sum[4] = {0,0,0,0};
    float mean[4] = {0,0,0,0};
    float variance[4] = {0,0,0,0};
    float skeweness[4] = {0,0,0,0};
    float kurtosis[4] = {0,0,0,0};
    int count = 0;
    float fulldct = 0;
    int pixel_p = 0;
    for(int i = 0; i < height; i+=window_size)
    {
        for(int j = 0; j < width; j+=window_size)
        {
            getImageBlock(DCT_ARRAY, i, j, width, window_size, block);
            block[0] = 0;

            for(int l = 0; l < window_size; l++) {
                for (int k = 0; k < window_size; k++) {
                    fulldct += block[(l * window_size) + k];
                    slices_sum[matrix_flag[l][k]-1] += block[(l * window_size) + k];
                }
            }

            for(int l = 0; l < 4; l++) {
                slices_sum[l] /= fulldct;
                EST[pixel_p] = slices_sum[l];
                pixel_p++;
            }

            count++;
        }
    }

    for(int l = 0; l < 4; l++)
    {
        pixel_p = 0;
        for(int i = 0; i < height/8; i++)
        {
            for(int j = 0; j < width/8; j++)
            {
                mean[l] += EST[pixel_p];
                pixel_p += 4;
            }
        }

        mean[l] /= count;

        pixel_p = 0;
        for(int i = 0; i < height/8; i++)
        {
            for(int j = 0; j < width/8; j++)
            {
                variance[l] += pow((EST[pixel_p] - mean[l]), 2);
                skeweness[l] += pow((EST[pixel_p]-mean[l]), 3);
                kurtosis[l] += pow((EST[pixel_p]-mean[l]), 4);
                pixel_p += 4;
            }
        }

        skeweness[l] = (skeweness[l]/count)/pow(sqrt(variance[l]/count), 3);
        kurtosis[l] = (kurtosis[l]/count)/pow(variance[l]/count, 2);
        variance[l] = variance[l]/(count-1);
        coefficients[z] = mean[l];
        coefficients[z+1] = variance[l];
        coefficients[z+2] = skeweness[l];
        coefficients[z+3] = kurtosis[l];
        z += 4;
    }

    for(int i = 0; i < 16; i++) printf("%f ", coefficients[i]);
    printf("\n");
}

