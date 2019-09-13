//
// Created by rostislav on 12.08.19.
//

#include <thread>

#ifndef IMAGEPROCESSING_EXE_DCT_FILT_H
#define IMAGEPROCESSING_EXE_DCT_FILT_H

namespace ImProcessing {
    float* DCT(char* image_, int width_, int height_, int wind_size_);
    char* ADCT(float* image_dct, int width_, int height_, int wind_size_);
    void getImageBlock(float* image, int i_, int j_, int image_width, int window_size, float* block);
    float abs(float zn);
}


#endif //IMAGEPROCESSING_EXE_DCT_FILT_H
