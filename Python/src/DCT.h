//
// Created by rostislav on 12.08.19.
//

#ifndef IMAGEPROCESSING_EXE_DCT_FILT_H
#define IMAGEPROCESSING_EXE_DCT_FILT_H
#include <cstdint>

namespace ImProcessing {
    float* DCT(uint8_t* image_, int width_, int height_, int wind_size_);
    uint8_t* ADCT(float* image_dct, int width_, int height_, int wind_size_);
    void getImageBlock(float* image, int i_, int j_, int image_width, int window_size, float* block);
    float abs(float zn);
}


#endif //IMAGEPROCESSING_EXE_DCT_FILT_H
