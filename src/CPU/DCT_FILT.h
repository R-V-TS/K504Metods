//
// Created by rostislav on 12.08.19.
//

#include "Image.h"
#include <vector>
#include <thread>

#ifndef IMAGEPROCESSING_EXE_DCT_FILT_H
#define IMAGEPROCESSING_EXE_DCT_FILT_H

namespace ImProcessing {
    class DCT_FILT {
    public:
        DCT_FILT(double* image_, int width, int height, int wind_size_, double threshold_);
        double* imageFilterMultithread();
        double* imageFilter();

    private:
        double* image;
        double* DCT_creator_mtx;
        double* DCT_creator_mtx_T;
        int image_width;
        int image_height;
        double* getImageBlock(int i, int j);
        double* MultiplyMatrix(double *matrix1, double *matrix2);

        // void DCT_thread(int w_i, int h_j);

        double threshold;
        int window_size;
        double* filtered_image;

        const int MAX_WIDTH = 2048;
        const int MAX_HEIGHT = 2048;

        //void dct(int w_i, int h_j, double* dct);
        void odct();
        int thread_number = 0;
        struct MyThread;


    };
}


#endif //IMAGEPROCESSING_EXE_DCT_FILT_H
