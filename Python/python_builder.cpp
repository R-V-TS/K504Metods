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

extern "C" void DCT_filt(double* image_array, int width_, int height_, int channels_, int window_size, double threshold_)
{
    if(window_size <= 32){
        double *pixel_ride = new double[width_*height_];;
        for(int channel = 0; channel < channels_; channel++)
        {
            int pixel_flag = 0;
            std::cout << "Channel " << channel << " in process" << "\n";
            for (int i = 0; i < height_; i++) {
                for (int j = channel; j < width_*channels_; j += channels_) {
                    *(pixel_ride + pixel_flag) = *(image_array+(i*width_*channels_)+j);
                    pixel_flag++;
                }
            }

            ImProcessing::DCT_FILT filter(pixel_ride, width_, height_, window_size, threshold_); // JPEG = 0.012
            double *pixel_result = filter.imageFilter();

            std::cout << "Channel " << channel << " success" << "\n";

            pixel_flag = 0;
            for (int i = 0; i < height_; i++)
                for(int j = channel; j < width_*channels_; j+= channels_)
                {
                    *(image_array+(i*width_*channels_)+j) = *(pixel_result + pixel_flag);
                    pixel_flag++;
                }
        }
    } else{
        return;
    }
}
