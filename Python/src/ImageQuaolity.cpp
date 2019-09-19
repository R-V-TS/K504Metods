#include <cmath>
#include <cstdint>
#include "ImageQuaolity.h"

namespace ImProcessing{
    float* MSE(uint8_t* P_image, uint8_t* Q_image, int width, int height, int channels){
        float* MSE_res = new float[channels];
        for(int ch = 0; ch < channels; ch++)
        {
            float MSE_local = 0;
            for(int i = 0; i < height; i++)
                for(int j = ch; j < width; j+=channels)
                {
                    MSE_local += pow(P_image[(i*width*channels) + j] - Q_image[(i*width*channels)+j], 2);
                }
            MSE_local /= (width*height);
            MSE_res[ch] = MSE_local;
        }
        return MSE_res;
    }

    float* PSNR(uint8_t* P_image, uint8_t* Q_image, int width, int height, int channels){
        float* PSNR_res = new float[channels];
        float* MSE_im = MSE(P_image, Q_image, width, height, channels);
        for(int ch = 0; ch < channels; ch++)
        {
            uint8_t max = 0;
            for(int i = 0; i < height; i++)
                for(int j = ch; j < width; j += channels)
                {
                    if(max < P_image[(i*width*channels)+j]) max = P_image[(i*width*channels)+j];
                }
            PSNR_res[ch] = 20*log10f(max/MSE_im[ch]);
        }
        return PSNR_res;
    }
}