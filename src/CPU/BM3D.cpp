//
// Created by rostislav on 24.09.19.
//

#include "BM3D.h"
#include <cmath>
#include "../DCT_Matrices.h"

namespace ImProcessing{

    float* normalize(uint8_t* image, unsigned int length, float max_of){
        float *res = new float[length];
        for(int i = 0; i < length; i++)
            res[i] = image[i]/max_of;
        return res;
    }

    float* BM3D_thr(float* image, int im_width, int im_height, float* Trans2D, float* Trans3D, float MSK, float* W, float* Wwind2D, int bmax, float nSa, int Nstep, int stepSN, std::string metric, float otherThr, float alpha)
    {
        float *res = new float[2];
        return res;
    }

    uint8_t* BM3D_CPU(uint8_t* image, int im_width, int im_height, float noise_sigma, std::string metric, float* MSK, float* W){
        uint8_t *result = new uint8_t[im_width*im_height];
        float* normalizing_image = normalize(image, im_width*im_height, 255);
        /* Set parameters for BM3D filter */
        int bsize = 8;               // розмер блока
        int bstep = 3;               // шаг обработки
        int bmax = 16;               // кол-во блоков в подобии
        int asize = 49;              // розмер области поиска
        float lambda_thr3D = 2.7;    // Порог для обработки
        float beta = 3.0;            // Параметр для окна Кайзера
        float threshold = 0;
        if( metric == "Euclidean") threshold = 3500;
        else if(metric == "BrayCurtis") threshold = 500;
        else threshold = 35000;
        // Block parameters
        float otherThr = threshold*bsize*bsize/(255*255);
        int stepFS = 1;
        int decLevel = 0;
        float alpha = 1.0;
        /* End block set parameters */
        float *Trans2D = &DCT_Creator8[0][0];
        float hlog = round(log2(bmax));


        return result;
    }

}