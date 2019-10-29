//
// Created by rostislav on 24.08.19.
//

#ifndef IMAGEPROCESSING_EXE_IMAGEQUAOLITY_H
#define IMAGEPROCESSING_EXE_IMAGEQUAOLITY_H
#include <cstdint>

namespace ImProcessing{
    float* MSE(uint8_t* P_image, uint8_t* Q_image, int width, int height, int channels);
    float* PSNR(uint8_t* P_image, uint8_t* Q_image, int width, int height, int channels);
    float* PSNRHVSM(uint8_t* P_image, uint8_t* Q_image, int width, int height, int channels);
}

#endif //IMAGEPROCESSING_EXE_IMAGEQUAOLITY_H
