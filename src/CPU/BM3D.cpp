//
// Created by rostislav on 16.10.19.
//

//
// Created by rostislav on 04.10.19.
//

#include <thread>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <cstdint>
#include <cmath>
#include <iostream>
#include "../DCT_Matrices.h"
#include "../Matrix.h"

void MultiplyMatrix(const float *matrix1, const float *matrix2, float *result, int window_size){
    const float *a, *b;
    for(int i = 0; i < window_size; i++) {
        float* c = result + (i * window_size);
        a = matrix1 + (i * window_size);
        for (int j = 0; j < window_size; j++) {
            b = matrix2 + j;
            c[j] = 0;
            for (int k = 0; k < window_size; k++) {
                c[j] += a[k] * b[k*window_size];
            }
        }
    }
}

void getImageBlock(float* image, int i_, int j_, int image_width, int window_size, float* block) {
    for (int i = 0; i < window_size; i++) {
        float *b = image + image_width * (i + i_) + j_;
        float *a = block + window_size * i;
        for (int j = 0; j < window_size; j++) {
            a[j] = b[j];
        }
    }
}

inline int max(const int one, const int two)
{
    if(one > two) return one;
    else return two;
}

inline int min(const int one, const int two)
{
    if(one < two) return one;
    else return two;
}

inline float sign_f(const float one)
{
    if(one > 0) return 1;
    else if(one == 0) return 0;
    else return -1;
}

uint8_t* BM3D(uint8_t* image_, int width_, int height_)
{
    using namespace ImProcessing;
    srand(time(NULL));
    uint16_t im_size = width_;
    uint16_t im_otr_size = width_+16;

    float* image = new float[im_size*im_size];
    float* im_otr = new float[(im_size+16)*(im_size+16)];
    float* Buff = new float[(im_size+16)*(im_size+16)];
    float* weig = new float[(im_size+16)*(im_size+16)];

    float W2v[64] ={0,0,0,0,0,0,0,0,
                    0,0,0,0,0,0,0,0,
                    0,0,0,0,0,0,0,0,
                    0,0,0,0,0,0,0,0,
                    0,0,0,0,0,0,0,0,
                    0,0,0,0,0,0,0,0,
                    0,0,0,0,0,0,0,0,
                    0,0,0,0,0,0,0,0};
    float Wwin2D[64] = {
            1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1
    };

    for(int i = 0; i < im_size*im_size; i++) image[i] = image_[i];
    unsigned int start = clock();

    for(int i = 0; i < im_size; i++)
    {
        float *a = im_otr + ((i+8)*im_otr_size) + 8;
        float *b = image + (i*im_size);
        for(int j = 0; j < im_size; j++){
            a[j] = b[j]/255;
        }
    }

    for(int i = 0; i < im_size; i++)
    {
        float *a = im_otr + i + 8;
        float *b = im_otr + i + 8;
        for(int j = 0; j < 8; j++)
        {
            a[(j*im_otr_size)] = b[(im_otr_size*(15-j))];
            a[((im_otr_size-j-1)*im_otr_size)] = b[((im_otr_size-16+j)*(im_otr_size))];
        }
    }

    for(int i = 0; i < 8; i++)
    {
        float *a = im_otr + i;
        float *a2 = im_otr + im_otr_size - i - 1;
        float *b = im_otr + 15 - i;
        float *c = im_otr + im_otr_size - 16 + i;
        for(int j = 0; j < im_otr_size; j++)
        {
            a[im_otr_size*j] = b[im_otr_size*j];
            a2[im_otr_size*j] = c[im_otr_size*j];
        }
    }

    for(int i = 0; i < im_otr_size*im_otr_size; i++)
    {
        Buff[i] = im_otr[i]/32;
        weig[i] = (float)1/32;
    }

    int AV_stride = 8*8;
    int ss1 = im_otr_size-8;
    int ss2 = im_otr_size-8;
    float* AV = new float[im_otr_size * im_otr_size * AV_stride];

    float* AV_pointer = AV;
    float* DCT_creator_mtx = &DCT_Creator8[0][0];
    float* DCT_creator_mtx_T = &DCT_Creator8_T[0][0];

    float *block = new float[64];
    float *temp = new float[64];

    for (int i = 0; i < ss1; i++) {
        for (int j = 0; j < ss2; j++) {
            getImageBlock(im_otr, i, j, im_otr_size, 8, block); // get block from image
            MultiplyMatrix(DCT_creator_mtx, block, temp, 8);
            MultiplyMatrix(temp, DCT_creator_mtx_T, block, 8);   // D*A*D'

            for(int z = 0; z < AV_stride; z++) {
                AV_pointer[z] = block[z];
            }
            AV_pointer = AV_pointer + AV_stride;
        }
    }
    free(image);
    //free(im_otr);

    int P = 3;
    int bmax = 16;
    int aSize = 49;
    int i, j;
    float alpha = 1;
    uint16_t startUp, startDown, startLeft, startRight;
    int nSa = 24; // (aSize-1)/2
    float otherThr = 34.44, tmp;
    int minVal_pointer, tmpi;
    float *currentColumn, *currentPatch;
    int* minIndex = new int[bmax];
    float* minValue = new float[bmax];
    float* dVec = new float[(nSa+1)];
    auto* dct3_vecTmp = new float[bmax*AV_stride];
    int elem_size; // size for Haar
    float* DCT3 = new float[AV_stride*bmax]; // Temp block for Haar

    float* bmax_dct3 = DCT3+AV_stride*bmax;
    float* DCT3_pointer = DCT3;
    float* one_mtx_pointer = dct3_vecTmp;
    float* two_mtx_pointer;

    for(int iii = 0, iCtr = P; iii < ss1; iii += P, iCtr++)
    {
        //if(iCtr < P) continue;
        //else iCtr = 0;
        for(int jjj = 0, jCtr = P; jjj < ss2; jjj+=P, jCtr++)
        {
            //if(jCtr < P) continue;
            //else jCtr = 0;
            int j = jjj + iii*ss2;

            int startUp = max(0, iii-nSa);
            int startLeft = max(0, jjj-nSa);
            int startDown = min(ss1, iii+nSa);
            int startRight = min(ss2, jjj+nSa);

            currentColumn = &AV[j*AV_stride];
            for(int z = 0; z < bmax; z++) minValue[z] = 0;

            for(int mmm = startUp; mmm < startDown; mmm++){
                currentPatch = &AV[(startLeft + (mmm*ss2))*AV_stride];

                for(int zzz = 0; zzz < startRight - startLeft; zzz++){
                    dVec[zzz] = 0;
                    for(int lll = 1; lll < AV_stride; lll++) {
                        dVec[zzz] += fabsf(currentColumn[lll] - currentPatch[lll]) / (fabsf(currentColumn[lll]) + fabsf(currentPatch[lll]));
                        //printf("%f (%f) / %f (%f) = %f\n", currentColumn[lll], fabsf(currentColumn[lll] - currentPatch[lll]), currentPatch[lll], (fabsf(currentColumn[lll]) + fabsf(currentPatch[lll])), dVec[zzz]);
                    }
                    //std::cout << dVec[zzz] << " ";
                    currentPatch += AV_stride;
                }
                //std::cout << "\n";


                minVal_pointer = 1;
                i = startLeft + mmm*ss2;
                for(int lll = 0; lll < nSa+1; lll++){
                    if(dVec[lll] < otherThr){
                        minValue[minVal_pointer] = dVec[lll];
                        minIndex[minVal_pointer] = i + lll;

                        for(int zzz = minVal_pointer; zzz > 0; zzz--)
                        {
                            if(minValue[zzz-1] > minValue[zzz])
                            {
                                tmp = minValue[zzz];
                                tmpi = minIndex[zzz];
                                minValue[zzz] = minValue[zzz-1];
                                minIndex[zzz] = minIndex[zzz-1];
                                minValue[zzz-1] = tmp;
                                minIndex[zzz-1] = tmpi;
                            }
                        }
                        minVal_pointer++;
                    }
                }
                minValue[0] = 0;
                minIndex[0] = j;

                /*for(int zzz = 0; zzz < minVal_pointer; zzz++)
                    printf("%f ", minValue[zzz]);
                printf("\n");*/

                --minVal_pointer;
                float *Trans3D = nullptr;
                float* Trans3D_T;
                /*if(minVal_pointer > 63){
                    elem_size = 64;
                }
                else if(minVal_pointer > 31) {
                    elem_size = 32;
                }
                else */
                if(minVal_pointer > 15) {
                    elem_size = 16;
                    Trans3D = &Trans3D_16[0][0];
                    Trans3D_T = &Trans3D_16_T[0][0];
                }
                else if(minVal_pointer > 8)
                {
                    elem_size = 8;
                    Trans3D = &Trans3D_8[0][0];
                    Trans3D_T = &Trans3D_8_T[0][0];
                }
                else if(minVal_pointer > 3){
                    elem_size = 4;
                    Trans3D = &Trans3D_4[0][0];
                    Trans3D_T = &Trans3D_4_T[0][0];
                }
                else if(minVal_pointer > 1) {
                    elem_size = 2;
                    Trans3D = &Trans3D_2[0][0];
                    Trans3D_T = &Trans3D_2_T[0][0];
                }
                else {
                    elem_size = 1;
                    Trans3D = &Trans3D_1;
                    Trans3D_T = &Trans3D_1;
                }


                for(int lll = 0; lll < elem_size; lll++) {
                    for(int ckk = 0; ckk < AV_stride; ckk++)
                    {
                        dct3_vecTmp[(lll*AV_stride) + ckk] = AV[minIndex[lll]*AV_stride + ckk];
                    }
                }

                // Haar transformation
                for(int h_a = 0; h_a < AV_stride; h_a++){
                    one_mtx_pointer = dct3_vecTmp + (h_a*elem_size); // Выбираем строку
                    DCT3_pointer = DCT3 + h_a*elem_size;
                    for(int w_a = 0; w_a < elem_size; w_a++){
                        two_mtx_pointer = Trans3D_T + w_a; //Выбираем столбец
                        DCT3_pointer[w_a] = 0;
                        for(int z_a = 0; z_a < elem_size; z_a++)
                        {
                            DCT3_pointer[w_a] += (one_mtx_pointer[z_a]*two_mtx_pointer[z_a*elem_size]);
                        }
                    }
                }
/*
                for(int fff = 0; fff < AV_stride; fff++){
                    for(int zzz = 0; zzz < elem_size; zzz++)
                    {
                        printf("%f ", DCT3[fff*elem_size + zzz]);
                    }
                    printf("\n");
                }
                printf("-----\n");
*/
                //printf("%i %i\n", minIndex[0], minIndex[1]);
                float* SURVIORS = new float[AV_stride*elem_size];
                float* S_pointer = SURVIORS;
                float W2v_pointer;
                int nhar = 0;
                for(int h_a = 0; h_a < AV_stride; h_a++)
                {
                    W2v_pointer = W2v[h_a];
                    DCT3_pointer = DCT3 + h_a*elem_size;
                    S_pointer = SURVIORS + h_a*elem_size;
                    for(int w_a = 0; w_a < elem_size; w_a++)
                    {
                        if(fabsf(DCT3_pointer[w_a]) > W2v_pointer)
                        {
                            nhar++;
                            S_pointer[w_a] = 1;
                        }
                        else S_pointer[w_a] = 0;
                    }
                }

                float temp;
                for(int h_a = 1; h_a < AV_stride*elem_size; h_a++)
                {
                    temp = DCT3[h_a];
                    DCT3[h_a] = sign_f(temp)*pow(fabsf(temp), alpha)*SURVIORS[h_a];
                }

/*                for(int fff = 0; fff < AV_stride; fff++){
                    for(int zzz = 0; zzz < elem_size; zzz++)
                    {
                        printf("%f ", DCT3[fff*elem_size + zzz]);
                    }
                    printf("\n");
                }
                printf("-----\n");
*/

                // Haar reverce transform
                for(int h_a = 0; h_a < AV_stride; h_a++){
                    one_mtx_pointer = DCT3 + (h_a*elem_size); // Выбираем строку
                    DCT3_pointer = dct3_vecTmp + h_a*elem_size;
                    for(int w_a = 0; w_a < elem_size; w_a++){
                        two_mtx_pointer = Trans3D + w_a; //Выбираем столбец
                        DCT3_pointer[w_a] = 0;
                        for(int z_a = 0; z_a < elem_size; z_a++)
                        {
                            DCT3_pointer[w_a] += (one_mtx_pointer[z_a]*two_mtx_pointer[z_a*elem_size]);
                        }
                    }
                }

/*                for(int fff = 0; fff < AV_stride; fff++){
                    for(int zzz = 0; zzz < elem_size; zzz++)
                    {
                        printf("%f ", dct3_vecTmp[fff*elem_size + zzz]);
                    }
                    printf("\n");
                }
                printf("-----\n");
*/

                // TODO:: Adaptive step
                if(nhar == 0) nhar = 1;
                float* win = new float[AV_stride];
                for(int h_a = 0; h_a < AV_stride; h_a++) win[h_a] = (1/(float)nhar)*Wwin2D[h_a];

                int i_in, j_in;
                float* temp_RDCT = new float[AV_stride];
                float* temp_RDCT_2 = new float[AV_stride];
                for(int immapped = 0; immapped < elem_size; immapped++)
                {
                    i_in = floorf(minIndex[immapped]/ss2);
                    j_in = minIndex[immapped] - (i_in)*ss2;

                    memcpy(temp_RDCT, dct3_vecTmp+(immapped*AV_stride), AV_stride* sizeof(float)); //Copy block to temp

                    MultiplyMatrix(DCT_creator_mtx_T, temp_RDCT, temp_RDCT_2, 8);
                    MultiplyMatrix(temp_RDCT_2, DCT_creator_mtx, temp_RDCT, 8);

/*                    for(int fff = 0; fff < 8; fff++){
                        for(int zzz = 0; zzz < 8; zzz++)
                        {
                            printf("%f ", im_otr[fff*im_otr_size + zzz]);
                        }
                        printf("\n");
                    }
                    printf("-----\n");
                    for(int fff = 0; fff < 8; fff++){
                        for(int zzz = 0; zzz < 8; zzz++)
                        {
                            printf("%f ", temp_RDCT[fff*8 + zzz]);
                        }
                        printf("\n");
                    }
                    printf("-----\n");
*/
                    for(int ii_b = i_in, i_b = 0; ii_b < 8; ii_b++, i_b++)
                    {
                        S_pointer = Buff + ii_b*im_otr_size;
                        minVal_pointer = i_b*8;
                        two_mtx_pointer = weig + ii_b*im_otr_size;
                        for(int jj_b = j_in, j_b = 0; jj_b < 8; jj_b++, j_b++)
                        {
                            S_pointer[jj_b] += max(min(temp_RDCT[minVal_pointer+j_b], 1), 0)*win[minVal_pointer+j_b];
                            two_mtx_pointer[jj_b] += win[minVal_pointer+j_b];
                        }
                    }
                }

                /*for(int h_a = 0; h_a < AV_stride; h_a++){
                    for(int w_a = 0; w_a < elem_size; w_a++) {
                        printf("%f ", DCT3[h_a*elem_size + w_a]);
                    }
                    printf("\n");
                }
                */

            }
        }
    }

    uint8_t* im_res = new uint8_t[im_size*im_size];
    int st = 8*im_otr_size+8;
    for(int ll = 0; ll < im_size; ll++)
        for(int ff = 0; ff < im_size; ff++)
            im_res[ll*im_size+ff] = (Buff[st+ll*im_otr_size+ff] / weig[st+ll*im_otr_size+ff])*255;
            //im_res[ll*im_size+ff] = im_otr[st+ll*im_otr_size+ff]*255;
        //im_res[ll] = (Buff[st+ll] / weig[st+ll])*255;
    unsigned int finish = clock();
    float time = ((float)(finish-start)/CLOCKS_PER_SEC);

    std::cout << "Time = " << std::to_string(time) << " s\n";

    for(int ll = 0; ll < 16; ll++)
    {
        for(int ff = 0; ff < 16; ff++)
        {
            printf("%i ", im_res[ll*im_size+ff]);
        }
        printf("\n");
    }

    return im_res;
}
