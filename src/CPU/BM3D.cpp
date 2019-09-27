//
// Created by rostislav on 24.09.19.
//

#include "BM3D.h"
#include <cmath>
#include <cstring>
#include "../DCT_Matrices.h"
#include "DCT.h"
#include <vector>
#include <thread>
#include <ctime>

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

    float* matrixDivisionBit(float* matrix1, float* matrix2, int length)
    {
        /*
         * Division matrix bit by bit
         */
        float* res = new float[length];
        for(int i = 0; i < length; i++)
            res[i] = matrix1[i]/matrix2[i];
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
                rep[(i*vector_width)+j] = vector[j];
        }
        return rep;
    }

    float* represent2vector(float* matrix, int n, int m)
    {
        /*
         * Only for matrix where n==m
         */
        float* res = new float[n*m];
        for(int i = 0; i < n; i++)
            for(int j = 0; j < m; j++)
                res[(j*m)+i] = matrix[(i*n)+j];
        return res;
    }

    void BM3D_coreThread(int i_start, int i_finish, int oP, int ss1, int ss2, int window_size, float* AV, float *MSK, float* W, int bmax, float nSa, std::string metric, float otherThr)
    {
        /*
         * Main thread function to run code in many thread
         */
        int oT = otherThr;
        unsigned int start_time = clock();
        int iCtr = (oP - 1);
        int cCtr = 0;
        int j = 0;
        int i = 0;
        int startUp, startDown, startLeft, startRight;
        int nm1 = window_size;
        int nm1_vec = nm1*nm1;

        for(int iii = i_start; iii < i_finish; iii++)
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
                startUp = max(0, iii-nSa);
                startLeft = max(0, jjj-nSa);
                startDown = min(ss1, iii+nSa+1);
                startRight = min(ss2, jjj+nSa+1);

                int current_coll = j*nm1_vec;
                //float* curr_coll = getVector(AV, j*nm1_vec, (j*nm1_vec)+(nm1_vec));

                for (int mmm = startUp; mmm < startDown; mmm++){
                    float* val_index = new float[bmax];
                    uint8_t* min_index = new uint8_t[bmax];
                    for(int z = 0; z < bmax; z++)
                    {
                        val_index[z] = 10000000;
                        min_index[z] = 0;
                    }

                    uint8_t length = startRight-startLeft;

                    int current_patch_ind = (startLeft+(mmm*ss2))*nm1_vec;
                    //float *current_patch = getVector(AV, (startLeft+(mmm*ss2))*nm1_vec, (startRight+(mmm*ss2))*nm1_vec);

                    double* dvec = new double[length];
                    int k_d = 0;
                    float BC_sum2 = 0;

                    for(int i_d = 0; i_d < nm1_vec; i_d++)   // нахождение вектора дальности
                    {
                        if (metric == "Euclidean") {
                            //for(int j_d = 0; j_d < length; j_d++) dvec[j_d] += pow((rep_current_mat[k_d + j_d] - rep_patch_mat[k_d+j_d]), 2);
                        } else if(metric == "BrayCurtis")
                        {
                            /*if(i_d != 0)
                            {
                                dvec[i_d] += abs((rep_current_mat[k_d] - rep_patch_mat[k_d]));
                                BC_sum2 += abs((rep_current_mat[k_d] + rep_patch_mat[k_d]));
                            }*/
                        } else if(metric == "Canberra") {  // Работает хорошо, за остальных не уверне
                            if(i_d != 0) {
                                for (int j_d = 0; j_d < length; j_d++) {
                                    dvec[j_d] += abs((AV[current_coll + i_d] - AV[current_patch_ind + (j_d*nm1_vec) + i_d])*W[i_d]) / (abs(AV[current_coll + i_d]*W[i_d]) + abs(AV[current_patch_ind + (j_d*nm1_vec)+i_d]*W[i_d]));
                                }
                            }
                        }
                        if(metric == "BrayCurtis") dvec[i_d] /= BC_sum2;
                    }

                    i = startLeft + mmm*ss2;
                    // Добавляем веса в массив и сортируем их в порядке убывания
                    for(int d_i=0; d_i < length; d_i++)
                    {
                        if(dvec[d_i] < val_index[0] && dvec[d_i] < oT)
                        {
                            if((i+d_i) == j)
                                continue;

                            val_index[0] = dvec[d_i];  // помещаем в первую ячейку
                            min_index[0] = i+d_i;

                            int jj = 0;
                            while(jj < bmax-2 && (val_index[jj] < val_index[jj+1]))
                            {
                                float temp_val = val_index[jj];
                                int temp_ind = min_index[jj];
                                val_index[jj] = val_index[jj+1];
                                min_index[jj] = min_index[jj+1];
                                min_index[jj+1] = temp_ind;
                                val_index[jj+1] = temp_val;
                                jj += 1;
                            }
                        }
                    }

                    for(int i_all = 0; i_all < length; i_all++) {
                        dvec[i_all] = 0;
                    }

                    //delete[] current_patch;
                    delete[] dvec;
                }

                //delete[] curr_coll;
            }
        }
        unsigned int finish_time = clock();
        //printf("\nThread %i end. Time = %f s\n", (ss1/i_start), (float) (finish_time-start_time)/CLOCKS_PER_SEC);
    }

    float* BM3D_thr(float* image, int im_width, int im_height, int window_size, float *MSK, float* W, float* Wwind2D, int bmax, float nSa, int Nstep, int stepSN, std::string metric, float otherThr, float alpha)
    {
        float *res = new float[2];
        /* Преобразование массивов в вектора не требуется, так как работаем с указателями */
        int oP = Nstep;
        int nm1 = window_size;
        int nm1_vec = nm1*nm1;

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

        float* temp = new float[nm1_vec];
        float* dct_temp = new float[nm1_vec];
        int AV_index = 0;
        int ss1 = height-nm1;
        int ss2 = width-nm1;
        float* AV = new float[ss1*ss2*nm1_vec];

        // Create array, where dct block as vector
        for(int i = 0; i < ss1; i++){
            for(int j = 0; j < ss2;j++)
            {
                getImageBlock(Bne, i, j, width, nm1, temp);
                MultiplyMatrix(&DCT_Creator8[0][0], temp, dct_temp, nm1);
                MultiplyMatrix(dct_temp, &DCT_Creator8_T[0][0], temp, nm1);
                for(int l = AV_index; l < AV_index+(nm1_vec); l++)
                {
                    AV[l] = temp[l-AV_index];
                }
                AV_index += (nm1_vec);
            }
        }

        std::vector<std::thread> thread_array;
        int tcols = 70;
        int thread_cols = (ss1-1)/tcols; // Делим все изображение для количества потоков
        int i = 0; // Указатель на край потока
        for(; i < tcols-1; i++)
            thread_array.push_back(std::thread(BM3D_coreThread, i*thread_cols, (i+1)*thread_cols-1, oP, ss1, ss2, window_size, AV, MSK, W, bmax, nSa, metric, otherThr));
        thread_array.push_back(std::thread(BM3D_coreThread, i*thread_cols, ss1-1, oP, ss1, ss2, window_size, AV, MSK, W, bmax, nSa, metric, otherThr));

        for(i = 0; i < thread_array.size(); i++)
            thread_array[i].join();

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