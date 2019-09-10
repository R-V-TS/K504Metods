//
// Created by rostislav on 12.08.19.
//

#include "DCT_FILT.h"
#include <iostream>

namespace ImProcessing{

    DCT_FILT::DCT_FILT(double* image_, int width_, int height_, int wind_size_, double threshold_)
    {
        image_width = width_;
        image_height = height_;
        image = image_;
        window_size = wind_size_;
        threshold = threshold_;
        switch (wind_size_) {
            case 2:
                DCT_creator_mtx = &DCT_Creator2[0][0];
                DCT_creator_mtx_T = &DCT_Creator2_T[0][0];
                break;
            case 4:
                DCT_creator_mtx = &DCT_Creator4[0][0];
                DCT_creator_mtx_T = &DCT_Creator4_T[0][0];
                break;
            case 8:
                DCT_creator_mtx = &DCT_Creator8[0][0];
                DCT_creator_mtx_T = &DCT_Creator8_T[0][0];
                break;
            case 16:
                DCT_creator_mtx = &DCT_Creator16[0][0];
                DCT_creator_mtx_T = &DCT_Creator16_T[0][0];
                break;
            case 32:
                DCT_creator_mtx = &DCT_Creator32[0][0];
                DCT_creator_mtx_T = &DCT_Creator32_T[0][0];
                break;
        }
        filtered_image = new double[image_width*image_height];
    }

    bool isThreadComplete = false;

    double* DCT_FILT::MultiplyMatrix(double *matrix1, double *matrix2){
        double *result = new double[window_size*window_size];
        double sum = 0;
        for(int i = 0; i < window_size; i++)
            for(int j = 0; j < window_size; j++) {
                for(int k = 0; k < window_size; k++){
                    sum += *(matrix1+(window_size*i)+k) * (*(matrix2+(window_size*k)+j));
                }
                *(result+(window_size*i)+j) = sum;
                sum = 0;
            }
        return result;
    }

    double* DCT_FILT::getImageBlock(int i_, int j_) {
        double *block = new double[window_size*window_size];
        for (int i = i_; i < window_size+i_; i++)
            for (int j = j_; j < window_size+j_; j++) {
                *(block + (window_size * (i-i_)) + (j-j_)) = *(image + (image_width * (i)) + (j));
            }
        return block;
    }

    double abs(double zn)
    {
        if(zn < 0)
            return zn * (-1);
        else
            return zn;
    }

    double* DCT_FILT::imageFilter() {
        double *block = new double[window_size]; // block from image
        double *temp = new double[window_size]; // temporary block

        double *im_temp = new double[image_width*image_height]; // temporary array for saved dct image coefficient
        double max = 0;

        for(int i = 0; i < image_height; i+=window_size) {
            for (int j = 0; j < image_width; j += window_size) {
                block = getImageBlock(i, j); // get block from image
                temp = MultiplyMatrix(DCT_creator_mtx, block);
                block = MultiplyMatrix(temp, DCT_creator_mtx_T);   // D*A*D'

                for (int k = i; k < window_size + i; k++) // Put DCT block to dct temp block
                    for (int t = j; t < window_size + j; t++) {
                        *(im_temp + (image_width * k) + t) = *(block + (window_size * (k - i)) + (t - j));
                        if(*(im_temp + (image_width*k) + t) > max)
                            max = *(im_temp + (image_width*k) + t);
                    }
            }
        }



        for(int i = 0; i < image_height; i++) {
            for (int j = 0; j < image_width; j++) {
                if(abs(*(im_temp + (image_width*i) + j)) < max * threshold)
                    *(im_temp + (image_width*i) + j) = 0;
            }
        }

        for(int i = 0; i < image_height; i+=window_size) {
            for (int j = 0; j < image_width; j += window_size) {
                double *block = new double[window_size*window_size];
                for (int k = i; k < window_size+i; k++)
                    for (int l = j; l < window_size+j; l++) {
                        *(block + (window_size * (k-i)) + (l-j)) = *(im_temp + (image_width * (k)) + (l));
                    }

                temp = MultiplyMatrix(block, DCT_creator_mtx);
                block = MultiplyMatrix(DCT_creator_mtx_T, temp);

                for(int k = i; k < window_size+i; k++) // Put DCT block to filtered
                    for (int t = j; t < window_size+j; t++)
                    {
                        *(filtered_image+(image_width*k)+t) = *(block+(window_size*(k-i))+(t-j));
                    }
            }
        }
        return filtered_image;
    }

    double* MultiplyMtxMulti(double *matrix1, double *matrix2, int window_size){
        double *result = new double[window_size*window_size];
        double sum = 0;
        for(int i = 0; i < window_size; i++)
            for(int j = 0; j < window_size; j++) {
                for(int k = 0; k < window_size; k++){
                    sum += *(matrix1+(window_size*i)+k) * (*(matrix2+(window_size*k)+j));
                }
                *(result+(window_size*i)+j) = sum;
                sum = 0;
            }
        return result;
    }

    static void DCT_thread(int w_i, int h_j, int block_size, int width, int threshold, double& DCT_mtx_, double& DCT_mtx_T_, double& image_)
    {
        double* image = (&image_);
        double* DCT_mtx = &DCT_mtx_;
        double* DCT_mtx_T = &DCT_mtx_T_;
        double* block = new double[block_size*block_size];

        for (int i = w_i; i < block_size+w_i; i++)
            for (int j = h_j; j < block_size+h_j; j++) {
                *(block + (block_size * (i-w_i)) + (j-h_j)) = *(image + (width * (i)) + (j));
            }

        block = MultiplyMtxMulti(DCT_mtx, block, block_size);
        block = MultiplyMtxMulti(block, DCT_mtx_T, block_size);   // D*A*D'

        for(int k = w_i; k < block_size+w_i; k++) // Put DCT block to filtered
            for (int t = h_j; t < block_size+h_j; t++)
            {
                if(*(block+(block_size*(k-w_i))+(t-h_j)) > threshold){
                    *(block+(block_size*(k-w_i))+(t-h_j)) = 0;
                }
            }

        block = MultiplyMtxMulti(block, DCT_mtx, block_size);
        block = MultiplyMtxMulti(DCT_mtx_T, block, block_size);

        for(int k = w_i; k < block_size+w_i; k++) // Put DCT block to filtered
            for (int t = h_j; t < block_size+h_j; t++)
            {
                *(image+(width*k)+t) = *(block+(block_size*(k-w_i))+(t-h_j));
            }
    }

    double* DCT_FILT::imageFilterMultithread() {
        /*
        std::vector<std::thread> dct_threads;

        for(int i = 0; i < image_height; i+=window_size) {
            for (int j = 0; j < image_width; j += window_size) {
                //dct_threads.push_back(std::thread(ImProcessing::thread_test, i, j));
                dct_threads.push_back(std::thread(ImProcessing::DCT_thread, i, j, window_size, image_width, threshold, std::ref(*DCT_creator_mtx), std::ref(*DCT_creator_mtx_T), std::ref(*image)));
            }
        }

        for(int i = 0; i < dct_threads.size(); i++)
            dct_threads[i].join();

        return filtered_image;
         */
        double z = 5;
        return &z;
    }
}