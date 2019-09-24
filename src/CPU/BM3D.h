//
// Created by rostislav on 24.09.19.
//

#ifndef IMAGEPROJCPU_BM3D_H
#define IMAGEPROJCPU_BM3D_H

#include <cstdint>
#include <cstdlib>
#include <string>

namespace ImProcessing
{
    uint8_t* BM3D_CPU(uint8_t* image, int im_width, int im_height, float noise_sigma, std::string metric, float* MSK, float* W);
}

#endif //IMAGEPROJCPU_BM3D_H
