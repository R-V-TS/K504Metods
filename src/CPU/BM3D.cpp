//
// Created by rostislav on 24.09.19.
//

#include "BM3D.h"
#include <cmath>
#include <cstring>
#include "../DCT_Matrices.h"
#include "DCT.h"

namespace ImProcessing{

    float Transform3D_2[2][2] = {
            {0.7071, 0.7071},
            {0.7071, -0.7071}
    };

    float Transform3D_4[4][4] = {
            {0.5000,    0.5000,    0.5000,    0.5000},
            {0.5000,    0.5000,   -0.5000,   -0.5000},
            {0.7071,   -0.7071,         0,         0},
            {0,         0,         0.7071,   -0.7071}
    };

    float Transform3D_8[8][8] = {
            {0.3536,    0.3536,    0.3536,    0.3536,    0.3536,    0.3536,    0.3536,    0.3536},
            {0.3536,    0.3536,    0.3536,    0.3536,   -0.3536,   -0.3536,   -0.3536,   -0.3536},
            {0.5000,    0.5000,   -0.5000,   -0.5000,         0,         0,         0,         0},
            {     0,         0,         0,         0,    0.5000,    0.5000,   -0.5000,   -0.5000},
            {0.7071,   -0.7071,         0,         0,         0,         0,         0,         0},
            {     0,         0,    0.7071,   -0.7071,         0,         0,         0,         0},
            {     0,         0,         0,         0,    0.7071,   -0.7071,         0,         0},
            {     0,         0,         0,         0,         0,         0,    0.7071,   -0.7071}
    };

    float Transform3D_16[16][16] = {
            {0.250, 0.250, 0.250, 0.250, 0.250, 0.250, 0.250, 0.250, 0.250, 0.250, 0.250, 0.250, 0.250, 0.250, 0.250, 0.250},
            {0.250, 0.250, 0.250, 0.250, 0.250, 0.250, 0.250, 0.250, -0.250, -0.250, -0.250, -0.250, -0.250, -0.250, -0.250, -0.250},
            {0.353553390593274,    0.353553390593274,    0.353553390593274,    0.353553390593274,    -0.353553390593274, -0.353553390593274, -0.353553390593274, -0.353553390593274,    0, 0, 0, 0, 0, 0, 0, 0},
            {0,     0,     0,     0,     0,     0,     0,     0,     0.353553390593274,    0.353553390593274,    0.353553390593274,    0.353553390593274, -0.353553390593274, -0.353553390593274, -0.353553390593274, -0.353553390593274},
            {0.50,   0.50,    -0.50,    -0.50,     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0,     0,     0,     0,     0.50,  0.50,     -0.50,    -0.50, 0, 0, 0, 0, 0, 0, 0, 0},
            {0,     0,     0,     0,     0,     0,     0,     0,     0.50,      0.50,    -0.50,    -0.50, 0, 0, 0, 0},
            {0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0.50,	0.50,	-0.50,	-0.50},
            {0.707106781186548,	-0.707106781186548,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
            {0,	0,	0.707106781186548,	-0.707106781186548,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
            {0,	0,	0,	0,	0.707106781186548,	-0.707106781186548,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
            {0,	0,	0,	0,	0,	0,	0.707106781186548,	-0.707106781186548,	0,	0,	0,	0,	0,	0,	0,	0},
            {0,	0,	0,	0,	0,	0,	0,	0,	0.707106781186548,	-0.707106781186548,	0,	0,	0,	0,	0,	0},
            {0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0.707106781186548,	-0.707106781186548,	0,	0,	0,	0},
            {0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0.707106781186548,	-0.707106781186548,	0,	0},
            {0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0.707106781186548,	-0.707106781186548}
    };

    float* normalize(uint8_t* image, unsigned int length, float max_of){
        float *res = new float[length];
        for(int i = 0; i < length; i++)
            res[i] = image[i]/max_of;
        return res;
    }

    float* getVector(float *Array_vectors, int start, int finish)
    {
        /*
         * Function for get vector from big array
         */
        float* vector = new float[finish-start];
        for(int i = start; i < finish; i++)
            vector[i-start] = Array_vectors[i];
        return vector;
    }

    int max(int first_, int second_)
    {
        if(first_ > second_) return first_;
        else return second_;
    }

    int min(int first_, int second_)
    {
        if(first_ < second_) return first_;
        else return second_;
    }

    float* matrixMultiplyBit(float* matrix1, float* matrix2, int length){
        /*
         * Multiply matrix bit by bit
         */
        float* res = new float[length];
        for(int i = 0; i < length; i++)
            res[i] = matrix1[i]*matrix2[i];
        return res;
    }

    float* repmat(float* vector, int vector_width, int width_rep)
    {
        /*
         *  copies vector to width_rep column
         */
        float* rep = new float[vector_width*width_rep];
        for(int i = 0; i < width_rep; i++)
        {
            for(int j = 0; j < vector_width; j++)
                rep[(i*width_rep)+j] = vector[j];
        }
        return rep;
    }

    float* BM3D_thr(float* image, int im_width, int im_height, int window_size, float *MSK, float* W, float* Wwind2D, int bmax, float nSa, int Nstep, int stepSN, std::string metric, float otherThr, float alpha)
    {
        float *res = new float[2];

        /* Преобразование массивов в вектора не требуется, так как работаем с указателями */
        int oP = Nstep;
        int nm1 = window_size;

        /* Adding rows and col to new image !Bne is new image! */
        float* Bne = new float[(im_width+16)*(im_height+16)];
        for(int i = 0, i_bne = nm1; i < im_height; i++, i_bne++)
            for(int j = 0, j_bne = nm1; j < im_width; j++, j_bne++)
                Bne[(i_bne*(im_width+16))+j_bne] = image[(i*im_width)+j];

        for(int l = 0, l_bne = nm1-1; l < nm1; l++, l_bne--)  // Add top row
            for(int k = 0, k_bne = nm1; k < im_width; k++, k_bne++)
                Bne[(l_bne*(im_width+16))+k_bne] = image[(l*im_width)+k];

        for(int l = im_height+16-nm1-1, l_bne = im_height+16-nm1; l > im_height+16-nm1*2-1; l--, l_bne++)  // Add bottom row
            for(int k = 0, k_bne = nm1; k < im_width; k++, k_bne++)
                Bne[(l_bne*(im_width+16))+k_bne] = image[(l*im_width)+k];

        for(int l = 0; l < im_height+16; l++)                                           // Add left coll
            for(int k = nm1, k_bne = nm1-1; k < nm1*2; k++, k_bne--)
                Bne[(l*(im_width+16))+k_bne] = Bne[(l*(im_width+16))+k];


        for(int l = 0; l < im_height+16; l++)                                          // Add right coll
            for(int k = im_width+16-nm1-1, k_bne = im_width+16-nm1; k > im_width+16-nm1*2-1; k--, k_bne++)
                Bne[(l*(im_width+16))+k_bne] = Bne[(l*(im_width+16))+k];

        int width = im_width + nm1*2;
        int height = im_height + nm1*2;
        /* End transform block */

        float* Buff = new float[width*height];
        float* weig = new float[width*height];

        for(int i = 0; i < width*height; i++)
        {
            Buff[i] = Bne[i]/32;
            weig[i] = 1/32;
        }

        float* temp = new float[nm1*nm1];
        float* dct_temp = new float[nm1*nm1];
        int AV_index = 0;
        int ss1 = height-nm1;
        int ss2 = width-nm1;
        float* AV = new float[ss1*ss2*nm1*nm1];
        // Create array, where dct block as vector
        for(int i = 0; i < ss1; i++){
            for(int j = 0; j < ss2;j++)
            {
                getImageBlock(Bne, i, j, width, nm1, temp);
                MultiplyMatrix(&DCT_Creator8[0][0], temp, dct_temp, nm1);
                MultiplyMatrix(dct_temp, &DCT_Creator8_T[0][0], temp, nm1);
                for(int l = AV_index; l < AV_index+(nm1*nm1); l++)
                {
                    AV[l] = temp[l-AV_index];
                }
                AV_index += (nm1*nm1);
            }
        }

        uint8_t iCtr = oP - 1;
        uint8_t cCtr = 0;
        uint8_t j = 0;
        uint8_t i = 0;
        uint8_t startUp, startDown, startLeft, startRight;
        for(int iii = 0; iii < (ss1-1); iii++)
        {
            iCtr++;
            if((iCtr < oP) && (iii < (ss1-1)))
                continue;
            else
                iCtr = 0;

            cCtr = oP-1;
            for(int jjj = 0; jjj < ss2-1; jjj++)
            {
                cCtr++;
                if((cCtr < oP) && (iii < (ss1-1)))
                    continue;
                else
                    cCtr = 0;

                j = jjj + iii*ss2;
                startUp = max(0, iii-int(nSa));
                startLeft = max(0, jjj-int(nSa));
                startDown = min(ss1, iii+int(nSa+1));
                startRight = min(ss2, jjj+int(nSa+1));

                float* curr_coll = getVector(AV, j*16, (j*16)+16);
                float* val_index = new float[bmax];
                uint8_t* min_index = new uint8_t[bmax];
                for(int z = 0; z < bmax; z++)
                {
                    val_index[z] = 10000000;
                    min_index[z] = 0;
                }

                float* rep_current_mat = repmat(matrixMultiplyBit(curr_coll, W, nm1*nm1), nm1*nm1, startRight-startLeft);

                for (uint8_t mmm = startUp; mmm < startDown; mmm++)
                {
                    uint8_t length = (startRight+(mmm*ss2))*16 - (startLeft+(mmm*ss2))*16;
                    float* rep_patch_mat = matrixMultiplyBit(getVector(AV, (startLeft+(mmm*ss2))*16, (startRight+(mmm*ss2))*16), repmat(W, nm1*nm1, startRight-startLeft), length);

                    for(int ir = 0; ir < length; ir++) printf("%f ", rep_patch_mat[ir]);
                }
            }
        }

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

        float Wwin2D[8][8] = {
                {0.0419777634056010,	0.102657454350536,	0.162602669831016,	0.199848511633810,	0.199848511633810,	0.162602669831016,	0.102657454350536,	0.0419777634056010},
                {0.102657454350536,	0.251050844036304,	0.397648059382450,	0.488733505447627,	0.488733505447627,	0.397648059382450,	0.251050844036304,	0.102657454350536},
                {0.162602669831016,	0.397648059382450,	0.629848426670745,	0.774121794899551,	0.774121794899551,	0.629848426670745,	0.397648059382450,	0.162602669831016},
                {0.199848511633810,	0.488733505447627,	0.774121794899551,	0.951442486736209,	0.951442486736209,	0.774121794899551,	0.488733505447627,	0.199848511633810},
                {0.199848511633810,	0.488733505447627,	0.774121794899551,	0.951442486736209,	0.951442486736209,	0.774121794899551,	0.488733505447627,	0.199848511633810},
                {0.162602669831016,	0.397648059382450,	0.629848426670745,	0.774121794899551,	0.774121794899551,	0.629848426670745,	0.397648059382450,	0.162602669831016},
                {0.102657454350536,	0.251050844036304,	0.397648059382450,	0.488733505447627,	0.488733505447627,	0.397648059382450,	0.251050844036304,	0.102657454350536},
                {0.0419777634056010,	0.102657454350536,	0.162602669831016,	0.199848511633810,	0.199848511633810,	0.162602669831016,	0.102657454350536,	0.0419777634056010}
        };

        for(int i = 0; i < bsize*bsize; i++)
            MSK[i] /= 255*lambda_thr3D*noise_sigma;


        float* res = BM3D_thr(normalizing_image, im_width, im_height, bsize, MSK, W, &Wwin2D[0][0], bmax, (asize-1)/2, bstep, stepFS, metric, otherThr, alpha);

        return result;
    }

}