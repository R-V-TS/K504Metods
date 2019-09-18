//
// Created by k504-r on 16.09.19.
//

#include "GPU_Functions.h"
#include <cstdio>

namespace ImProcessing {
    float *DCT_GPU(char *image, int width, int height, int channels, int window_size) {
        printf("GPU not support cuda\n");
        return nullptr;
    }

    char *ADCT_GPU(float *DCT_array, int width, int height, int channels, int window_size) {
        printf("GPU not support cuda\n");
        return nullptr;
    }
}
