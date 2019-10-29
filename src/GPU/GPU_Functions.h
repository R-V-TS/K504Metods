//
// Created by k504-r on 16.09.19.
//

#ifndef CUDA__GPU_FUNCTIONS_H
#define CUDA__GPU_FUNCTIONS_H

namespace ImProcessing {
    float *DCT_GPU(unsigned char *image, int width, int height, int channels, int window_size);
    unsigned char *ADCT_GPU(float* DCT_array, int width, int height, int channels, int window_size);
}

#endif //CUDA__GPU_FUNCTIONS_H
