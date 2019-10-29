//
// Created by k504-r on 16.09.19.
//
#include <cuda.h>
#include "GPU_Functions.h"
#include "../DCT_Matrices.h"

namespace ImProcessing {
    __device__ void matrixMultiply(float *matrix1, float *matrix2, int window_size, float *result) {
        float sum = 0;
        for (int i = 0; i < window_size; i++)
            for (int j = 0; j < window_size; j++) {
                for (int k = 0; k < window_size; k++) {
                    sum += *(matrix1 + (window_size * i) + k) * (*(matrix2 + (window_size * k) + j));
                }
                *(result + (window_size * i) + j) = sum;
                sum = 0;
            }
    }

    __global__ void gDCT(unsigned char *image, int *image_width, int *pixel_delay, float *DCT_Creator, float *DCT_Creator_T, int *window_size, float *result) {
        int x = (blockIdx.x) * (*window_size) * (*pixel_delay) + (threadIdx.x);
        int y = (blockIdx.y) * (*window_size);

        float *image_block = new float[(*window_size) * (*window_size)];
        int p_dy = y;
        for (int k = 0; k < (*window_size); k++, p_dy++) {
            for (int t = 0, p_dx = x; t < (*window_size); t++, p_dx += (*pixel_delay)) {
                image_block[((*window_size) * (k)) + (t)] = image[((*image_width)*(*pixel_delay) * p_dy) + p_dx];
            }
        }

        float *temp_block = new float[(*window_size) * (*window_size)];
        matrixMultiply(DCT_Creator, image_block, (*window_size), temp_block);
        matrixMultiply(temp_block, DCT_Creator_T, (*window_size), image_block);

        p_dy = y;

        for (int k = 0; k < (*window_size); k++, p_dy++) {
            for (int t = 0, p_dx = x; t < (*window_size); t++, p_dx += (*pixel_delay)) {
                result[((*image_width)*(*pixel_delay) * p_dy) + p_dx] = image_block[((*window_size) * (k)) + (t)];
            }
        }

        free(temp_block);
        free(image_block);
    }


    __global__ void gADCT(float *image, int *image_width, int *pixel_delay, float *DCT_Creator, float *DCT_Creator_T, int *window_size, unsigned char *result) {
        int x = (blockIdx.x) * (*window_size) * (*pixel_delay) + (threadIdx.x);
        int y = (blockIdx.y) * (*window_size);

        float *image_block = new float[(*window_size) * (*window_size)];
        int p_dy = y;
        for (int k = 0; k < (*window_size); k++, p_dy++) {
            for (int t = 0, p_dx = x; t < (*window_size); t++, p_dx += (*pixel_delay)) {
                image_block[((*window_size) * (k)) + (t)] = image[((*image_width)*(*pixel_delay) * p_dy) + p_dx];
            }
        }

        float *temp_block = new float[(*window_size) * (*window_size)];

        matrixMultiply(image_block, DCT_Creator, (*window_size), temp_block);
        matrixMultiply(DCT_Creator_T, temp_block, (*window_size), image_block);

        p_dy = y;

        for (int k = 0; k < (*window_size); k++, p_dy++) {
            for (int t = 0, p_dx = x; t < (*window_size); t++, p_dx += (*pixel_delay)) {
                result[((*image_width)*(*pixel_delay) * p_dy) + p_dx] = image_block[((*window_size) * (k)) + (t)];
            }
        }

        free(temp_block);
        free(image_block);
    }

    float *DCT_GPU(unsigned char *image, int width, int height, int channels, int window_size) {
        float *result = new float[width * height * channels];

        unsigned char *image_dev;
        float *DCT_creator_dev;
        float *DCT_creator_T_dev;
        float *result_dev;
        int *window_size_dev;
        int *image_width_dev;
        int *pixel_delay_dev;


        cudaMalloc((void **) &image_dev, sizeof(unsigned char) * width * height * channels); // 256 is image size
        cudaMalloc((void **) &DCT_creator_dev, sizeof(float) * window_size * window_size);
        cudaMalloc((void **) &DCT_creator_T_dev, sizeof(float) * window_size * window_size);
        cudaMalloc((void **) &window_size_dev, sizeof(int));
        cudaMalloc((void **) &image_width_dev, sizeof(int));
        cudaMalloc((void **) &pixel_delay_dev, sizeof(int));
        cudaMalloc((void **) &result_dev, sizeof(float) * width * height * channels);

        cudaMemcpyAsync(image_dev, image, sizeof(unsigned char) * width * height * channels, cudaMemcpyHostToDevice);   //copy image to videomemory
        cudaMemcpy(window_size_dev, &window_size, sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(image_width_dev, &width, sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(pixel_delay_dev, &channels, sizeof(int), cudaMemcpyHostToDevice);

        switch (window_size) {   // copy DCT creator to video
            case 2:
                cudaMemcpyAsync(DCT_creator_dev, ImProcessing::DCT_Creator2, sizeof(float) * window_size * window_size,
                                cudaMemcpyHostToDevice);
                cudaMemcpyAsync(DCT_creator_T_dev, ImProcessing::DCT_Creator2_T,
                                sizeof(float) * window_size * window_size, cudaMemcpyHostToDevice);
                break;
            case 4:
                cudaMemcpyAsync(DCT_creator_dev, ImProcessing::DCT_Creator4, sizeof(float) * window_size * window_size,
                                cudaMemcpyHostToDevice);
                cudaMemcpyAsync(DCT_creator_T_dev, ImProcessing::DCT_Creator4_T,
                                sizeof(float) * window_size * window_size, cudaMemcpyHostToDevice);
                break;
            case 8:
                cudaMemcpyAsync(DCT_creator_dev, ImProcessing::DCT_Creator8, sizeof(float) * window_size * window_size,
                                cudaMemcpyHostToDevice);
                cudaMemcpyAsync(DCT_creator_T_dev, ImProcessing::DCT_Creator8_T,
                                sizeof(float) * window_size * window_size, cudaMemcpyHostToDevice);
                break;
            case 16:
                cudaMemcpyAsync(DCT_creator_dev, ImProcessing::DCT_Creator16, sizeof(float) * window_size * window_size,
                                cudaMemcpyHostToDevice);
                cudaMemcpyAsync(DCT_creator_T_dev, ImProcessing::DCT_Creator16_T,
                                sizeof(float) * window_size * window_size, cudaMemcpyHostToDevice);
                break;
            case 32:
                cudaMemcpyAsync(DCT_creator_dev, ImProcessing::DCT_Creator32, sizeof(float) * window_size * window_size,
                                cudaMemcpyHostToDevice);
                cudaMemcpyAsync(DCT_creator_T_dev, ImProcessing::DCT_Creator32_T,
                                sizeof(float) * window_size * window_size, cudaMemcpyHostToDevice);
                break;
        }

        dim3 grid((width) / (window_size), height / window_size);

        gDCT<<< grid, channels >>>(image_dev, image_width_dev, pixel_delay_dev, DCT_creator_dev, DCT_creator_T_dev, window_size_dev, result_dev);

        cudaMemcpyAsync(result, result_dev, width * height * channels * sizeof(float), cudaMemcpyDeviceToHost);

        cudaFree(result_dev);
        cudaFree(image_dev);
        cudaFree(image_width_dev);
        cudaFree(window_size_dev);
        cudaFree(DCT_creator_dev);
        cudaFree(DCT_creator_T_dev);

        return result;
    }

    unsigned char *ADCT_GPU(float *image, int width, int height, int channels, int window_size) {
        unsigned char *result = new unsigned char[width * height * channels];

        float *image_dev;
        float *DCT_creator_dev;
        float *DCT_creator_T_dev;
        unsigned char *result_dev;
        int *window_size_dev;
        int *image_width_dev;
        int *pixel_delay_dev;


        cudaMalloc((void **) &image_dev, sizeof(float) * width * height * channels); // 256 is image size
        cudaMalloc((void **) &DCT_creator_dev, sizeof(float) * window_size * window_size);
        cudaMalloc((void **) &DCT_creator_T_dev, sizeof(float) * window_size * window_size);
        cudaMalloc((void **) &window_size_dev, sizeof(int));
        cudaMalloc((void **) &image_width_dev, sizeof(int));
        cudaMalloc((void **) &pixel_delay_dev, sizeof(int));
        cudaMalloc((void **) &result_dev, sizeof(unsigned char) * width * height * channels);

        cudaMemcpyAsync(image_dev, image, sizeof(float) * width * height * channels,
                        cudaMemcpyHostToDevice);   //copy image to videomemory
        cudaMemcpy(window_size_dev, &window_size, sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(image_width_dev, &width, sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(pixel_delay_dev, &channels, sizeof(int), cudaMemcpyHostToDevice);

        switch (window_size) {   // copy DCT creator to video
            case 2:
                cudaMemcpyAsync(DCT_creator_dev, ImProcessing::DCT_Creator2, sizeof(float) * window_size * window_size,
                                cudaMemcpyHostToDevice);
                cudaMemcpyAsync(DCT_creator_T_dev, ImProcessing::DCT_Creator2_T,
                                sizeof(float) * window_size * window_size, cudaMemcpyHostToDevice);
                break;
            case 4:
                cudaMemcpyAsync(DCT_creator_dev, ImProcessing::DCT_Creator4, sizeof(float) * window_size * window_size,
                                cudaMemcpyHostToDevice);
                cudaMemcpyAsync(DCT_creator_T_dev, ImProcessing::DCT_Creator4_T,
                                sizeof(float) * window_size * window_size, cudaMemcpyHostToDevice);
                break;
            case 8:
                cudaMemcpyAsync(DCT_creator_dev, ImProcessing::DCT_Creator8, sizeof(float) * window_size * window_size,
                                cudaMemcpyHostToDevice);
                cudaMemcpyAsync(DCT_creator_T_dev, ImProcessing::DCT_Creator8_T,
                                sizeof(float) * window_size * window_size, cudaMemcpyHostToDevice);
                break;
            case 16:
                cudaMemcpyAsync(DCT_creator_dev, ImProcessing::DCT_Creator16, sizeof(float) * window_size * window_size,
                                cudaMemcpyHostToDevice);
                cudaMemcpyAsync(DCT_creator_T_dev, ImProcessing::DCT_Creator16_T,
                                sizeof(float) * window_size * window_size, cudaMemcpyHostToDevice);
                break;
            case 32:
                cudaMemcpyAsync(DCT_creator_dev, ImProcessing::DCT_Creator32, sizeof(float) * window_size * window_size,
                                cudaMemcpyHostToDevice);
                cudaMemcpyAsync(DCT_creator_T_dev, ImProcessing::DCT_Creator32_T,
                                sizeof(float) * window_size * window_size, cudaMemcpyHostToDevice);
                break;
        }

        dim3 grid((width) / (window_size), height / window_size);

        gADCT << < grid, channels >> >
                         (image_dev, image_width_dev, pixel_delay_dev, DCT_creator_dev, DCT_creator_T_dev, window_size_dev, result_dev);

        cudaMemcpyAsync(result, result_dev, width * height * channels * sizeof(unsigned char), cudaMemcpyDeviceToHost);

        cudaFree(result_dev);
        cudaFree(image_dev);
        cudaFree(image_width_dev);
        cudaFree(window_size_dev);
        cudaFree(DCT_creator_dev);
        cudaFree(DCT_creator_T_dev);

        return result;
    }
}




