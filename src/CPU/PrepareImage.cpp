//
// Created by rostislav on 20.07.19.
//

#include "PrepareImage.h"
#include "SETUP_FILE.h"
#include <cmath>

namespace ImProcessing{
    PrepareImage::PrepareImage(ImProcessing::Image* image_) {
        image = image_;
        MatrixOperation = new Matrix(BLOCK_OP_W, BLOCK_OP_H);
    }

    double* PrepareImage::calculateCoefficient() {
        DCT_frequency = new double[int(image->getDimension()->getHeight()/8)*int(image->getDimension()->getWidth()/8)*4];
        double* block;
        double* temp;
        int flag_dct_freq = 0;
        double freq_sum = 0;

        for(int i = 0; i < image->getDimension()->getWidth(); i+=BLOCK_OP_W)
            for(int j = 0; j < image->getDimension()->getHeight(); j+=BLOCK_OP_H)
            {
                block = image->getBlockPixel(i, j, BLOCK_OP_W, BLOCK_OP_H);
                temp = MatrixOperation->matrixMultiply(DCT_CREATOR[0], block);
                block = MatrixOperation->matrixMultiply(temp, DCT_CREATOR_Transpon[0]); //Block DCT
                MatrixOperation->matrixPower(block); // DCT^2
                *(block) = 0; // block[0][0] = 0
                double fulldct = MatrixOperation->matrixFullSum(block); // full sum DCT 8*8

                if(fulldct != 0) {
                    // Calculate coefficient
                    for (int k = 0; k < 4; k++)
                    {
                        freq_sum = 0;
                        for(int t = 0; t < 8; t++)
                            for (int l = 0; l < 8; l++) {
                                if ((k + 1) == flags[t][l])
                                    freq_sum += *(block + (t * 8) + l);
                            }

                        *(DCT_frequency+flag_dct_freq) = freq_sum/fulldct;
                        flag_dct_freq++;
                    }
                }
            }
        return calculateStat(DCT_frequency, int(image->getDimension()->getHeight()/8)*int(image->getDimension()->getWidth()/8)*4);
    }


    double power(double ch, int sh){
        double res = 1;
        for(int i = 0; i < sh; i++)
            res *= ch;
        return res;
    }

    double* PrepareImage::calculateStat(double *DCT_freq, int length) {
        double* result = new double[16];
        double mean[4] = {0,0,0,0};
        int sizeDCT = (image->getDimension()->getWidth()-7)*(image->getDimension()->getHeight()-7);

        for(int i = 0; i < length; i+=4)
        {
            mean[2] += *(DCT_freq+i+2);
            mean[3] += *(DCT_freq+i+3);
            mean[1] += *(DCT_freq+i+1);
            mean[0] += *(DCT_freq+i);
        }

        for(int i = 0; i < 4;i++)
        {
            mean[i] /= sizeDCT;
            *(result+i) = mean[i];
        }

        double sum_variance[4] = {0,0,0,0};
        double sum_skeweness[4] = {0,0,0,0};
        double sum_kurtosis[4] = {0,0,0,0};

        for(int i = 0; i < length; i+=4){
            sum_variance[0] += power(*(DCT_freq+i), 2);
            sum_skeweness[0] += power(*(DCT_freq+i) - mean[0], 3);
            sum_kurtosis[0] += power(*(DCT_freq+i) - mean[0], 4);
            sum_variance[1] += power(*(DCT_freq+i+1), 2);
            sum_skeweness[1] += power(*(DCT_freq+i+1) - mean[1], 3);
            sum_kurtosis[1] += power(*(DCT_freq+i+1) - mean[1], 4);
            sum_variance[2] += power(*(DCT_freq+i+2), 2);
            sum_skeweness[2] += power(*(DCT_freq+i+2) - mean[2], 3);
            sum_kurtosis[2] += power(*(DCT_freq+i+2) - mean[2], 4);
            sum_variance[3] += power(*(DCT_freq+i+3), 2);
            sum_skeweness[3] += power(*(DCT_freq+i+3) - mean[3], 3);
            sum_kurtosis[3] += power(*(DCT_freq+i+3) - mean[3], 4);
        }

        for(int i = 0; i < 4; i++){
            *(result+1+(i*4)) = sum_variance[i]/(sizeDCT-1);
            *(result+2+(i*4)) = (sum_skeweness[i]/sizeDCT)/power(sqrt(sum_variance[i]/sizeDCT), 3);
            *(result+3+(i*4)) = (sum_kurtosis[i]/sizeDCT)/power(sum_variance[i]/sizeDCT, 2);
        }
        return result;
    }
}
