//
// Created by rostislav on 10.09.19.
//

#include "RAWImage.h"
#include "ERROR_Printer.h"
#include <opencv2/opencv.hpp>
#include <cmath>
#include <vector>
#include <random>
#include "CPU/ImageQuaolity.h"

// Metods
#include "CPU/DCT.h"
#include "GPU/GPU_Functions.h"

namespace ImProcessing
{
    RAWImage::RAWImage(std::string filename, ImProcessing::ImagePixelType type) {
        cv::Mat imCV = cv::imread(filename, 1);
        if(!imCV.data)
        {
            PrintError("Path don't have image!");
            return;
        } else{
            uint8_t *image_pxl = (uint8_t*)imCV.data;
            createImageArray(image_pxl, imCV.cols, imCV.rows, type);
        }

    }

    RAWImage::RAWImage(uint8_t *image_pixel, int image_width, int image_height, ImProcessing::ImagePixelType type) {
        createImageArray(image_pixel, image_width, image_height, type);
    }

    RAWImage::~RAWImage() {
        //delete image;
        width = 0;
        height = 0;
    }

    void RAWImage::createImageArray(uint8_t *image_pixel, int image_width, int image_height,
                                    ImProcessing::ImagePixelType type) {
        switch(type)
        {
            case 1:
                image = new uint8_t[image_height*image_width*3];
                original_image = new uint8_t[image_height*image_width*3];
                for(int i = 0; i < image_width*image_height*3; i++){
                    image[i] = image_pixel[i];
                    original_image[i] = image[i];
                }
                width = image_width;
                height = image_height;
                channels = 3;
                ImType = type;
                //delete image_pixel;
                break;
            case 4:
                width = image_width;
                height = image_height;
                channels = 3;
                ImType = type;
                B = new uint8_t[width*height];
                G = new uint8_t[width*height];
                R = new uint8_t[width*height];
                original_image = new uint8_t[image_height*image_width*3];
                int p = 0;
                for(int i = 0; i < image_width*image_height*channels; i+=channels)
                {
                    B[p] = image_pixel[i];
                    G[p] = image_pixel[i+1];
                    R[p] = image_pixel[i+2];
                    original_image[i] = B[p];
                    original_image[i+1] = G[p];
                    original_image[i+2] = R[p];
                    p++;
                }
                //free(image_pixel);
                break;
        }
    }

    void RAWImage::show() {
        if(image != NULL || (R != NULL && G != NULL && B != NULL))
        {
            cv::Mat image_show(width, height, CV_8UC3);
            switch (ImType)
            {
                case 1:
                    for(int i = 0; i < width*height*channels; i++)
                        image_show.data[i] = image[i];
                    break;
                case 4:
                    for(int i = 0; i < width*height*channels; i+=3)
                    {
                        image_show.data[i] = B[i/channels];
                        image_show.data[i+1] = G[i/channels];
                        image_show.data[i+2] = R[i/channels];
                    }
                    break;
            }
            cv::namedWindow("ImageShow");
            cv::imshow("ImageShow", image_show);
            cv::waitKey(0);
        }
    }

    void RAWImage::save(std::string destination_folder) {
        if(image != NULL || (R != NULL && G != NULL && B != NULL))
        {
            cv::Mat image_show(width, height, CV_8UC3);
            switch (ImType)
            {
                case 1:
                    for(int i = 0; i < width*height*channels; i++)
                        image_show.data[i] = image[i];
                    break;
                case 4:
                    for(int i = 0; i < width*height*channels; i+=3)
                    {
                        image_show.data[i] = B[i/channels];
                        image_show.data[i+1] = G[i/channels];
                        image_show.data[i+2] = R[i/channels];
                    }
                    break;
            }
            cv::imwrite(destination_folder + "result.png", image_show);
        }
    }

    void RAWImage::ApplyDCT(int window_size, double threshold, bool device) {
        if(device)
        {
            if(ImType != TYPE_3ARRAY)
            {
                transfer2OtherType(TYPE_3ARRAY);
            }

            float* DCT_ARRAY;
            float max = 0;
            DCT_ARRAY = DCT(B, width, height, window_size);
            free(B);
            for(int i = 0; i < width*height; i++)
            {
                if(DCT_ARRAY[i] > max) max = DCT_ARRAY[i];
            }
            for(int i = 0; i < width*height; i++)
            {
                if(abs(DCT_ARRAY[i]) < max*threshold) DCT_ARRAY[i] = 0;
            }
            B = ADCT(DCT_ARRAY, width, height, window_size);
            free(DCT_ARRAY);
            DCT_ARRAY = DCT(G, width, height, window_size);
            free(G);
            for(int i = 0; i < width*height; i++)
            {
                if(DCT_ARRAY[i] > max) max = DCT_ARRAY[i];
            }
            for(int i = 0; i < width*height; i++)
            {
                if(abs(DCT_ARRAY[i]) < max*threshold) DCT_ARRAY[i] = 0;
            }
            G = ADCT(DCT_ARRAY, width, height, window_size);
            free(DCT_ARRAY);
            DCT_ARRAY = DCT(R, width, height, window_size);
            free(R);
            for(int i = 0; i < width*height; i++)
            {
                if(DCT_ARRAY[i] > max) max = DCT_ARRAY[i];
            }
            for(int i = 0; i < width*height; i++)
            {
                if(abs(DCT_ARRAY[i]) < max*threshold) DCT_ARRAY[i] = 0;
            }
            R = ADCT(DCT_ARRAY, width, height, window_size);
            free(DCT_ARRAY);
        } else {
            if(ImType != TYPE_BGR)
            {
                transfer2OtherType(TYPE_BGR);
            }

            float *DCT_ARRAY = DCT_GPU(image, width, height, channels, window_size);
            float max = 0;
            for(int i = 0; i < width*height; i++)
            {
                if(DCT_ARRAY[i] > max) max = DCT_ARRAY[i];
            }
            for(int i = 0; i < width*height; i++)
            {
                if(abs(DCT_ARRAY[i]) < max*threshold) DCT_ARRAY[i] = 0;
            }
            image = ADCT_GPU(DCT_ARRAY, width, height, channels, window_size);
            free(DCT_ARRAY);
        }
    }

    float* RAWImage::DCTCoefficients(bool device) {
        if(ImType != TYPE_3ARRAY)
        {
            transfer2OtherType(TYPE_3ARRAY);
        }

        int window_size = 8;
        int matrix_flag[8][8] = {
                {0,  1,  1,  1,  1,  2,  2,  2},
                {1,  1,  1,  1,  2,  2,  2,  3},
                {1,  1,  1,  2,  2,  2,  3,  3},
                {1,  1,  2,  2,  2,  3,  3,  3},
                {1,  2,  2,  2,  3,  3,  3,  4},
                {2,  2,  2,  3,  3,  3,  4,  4},
                {2,  2,  3,  3,  3,  4,  4,  4},
                {2,  3,  3,  3,  4,  4,  4,  4}
        };

        float* EST = new float[(width/8)*(height/8)*4];
        float* coefficients = new float[16];
        float* block = new float[window_size*window_size];
        float* DCT_ARRAY;
        if(device) DCT_ARRAY = DCT(B, width, height, window_size);
        else DCT_ARRAY = DCT_GPU(B, width, height, 1, window_size);
        int z = 0;

        double slices_sum[4] = {0,0,0,0};
        double mean = 0;
        double variance = 0;
        double skeweness = 0;
        double kurtosis = 0;
        int count = 0;
        double fulldct = 0;
        int pixel_p = 0;
        for(int i = 0; i < height; i+=window_size)
        {
            for(int j = 0; j < width; j+=window_size)
            {
                fulldct = 0;
                getImageBlock(DCT_ARRAY, i, j, width, window_size, block);
                block[0] = 0;

                for(int l = 0; l < window_size; l++) {
                    for (int k = 0; k < window_size; k++) {
                        block[(l*window_size) + k] = block[(l*window_size) + k] * block[(l*window_size) + k];
                        fulldct += block[(l * window_size) + k];
                        slices_sum[matrix_flag[l][k]-1] += block[(l * window_size) + k];
                    }
                }

                for(int l = 0; l < 4; l++) {
                    slices_sum[l] /= fulldct;
                    EST[pixel_p] = slices_sum[l];
                    pixel_p++;
                    slices_sum[l] = 0;
                }
            }
        }

        for(int l = 0; l < 4; l++)
        {
            mean = 0;
            pixel_p = l;
            for(int i = 0; i < height/8; i++)
            {
                for(int j = 0; j < width/8; j++)
                {
                    mean += EST[pixel_p];
                    pixel_p += 4;
                    count++;
                }
            }
            count = width*height;
            mean /= count;

            pixel_p = l;
            variance = 0;
            skeweness = 0;
            kurtosis = 0;
            for(int i = 0; i < height/8; i++)
            {
                for(int j = 0; j < width/8; j++)
                {
                    variance += pow((EST[pixel_p] - mean), 2);
                    skeweness += pow((EST[pixel_p]-mean), 3);
                    kurtosis += pow((EST[pixel_p]-mean), 4);
                    pixel_p += 4;
                }
            }

            skeweness = (skeweness/count)/pow(sqrt(variance/count), 3);
            kurtosis = (kurtosis/count)/pow(variance/count, 2);
            variance = variance/(count-1);
            coefficients[z] = mean;
            coefficients[z+1] = variance;
            coefficients[z+2] = skeweness;
            coefficients[z+3] = kurtosis;
            z += 4;
        }
        free(EST);
        return coefficients;
    }

    void RAWImage::AddNoise(float mu, float sigma) {
        std::default_random_engine generator; //Random generator engine
        std::normal_distribution<float> distribution(mu, sigma); //distribution vector

        if(ImType != TYPE_3ARRAY)
            transfer2OtherType(TYPE_3ARRAY);

        for(int i = 0; i < height; i++)
            for(int j = 0; j < width; j++)
            {
                B[(i * width) + j] += distribution(generator);
                G[(i * width) + j] += distribution(generator);
                R[(i * width) + j] += distribution(generator);
            }
    }

    void RAWImage::printImageCharacteristics(){
        if(ImType != TYPE_BGR)
            transfer2OtherType(TYPE_BGR);

        printf("\nMetrics: \n");
        float* PSNR_ = PSNR(original_image, image, width, height, channels);
        printf("PSNR: ");
        for(int i = 0; i < 3; i++)
            printf("%f ", PSNR_[i]);
        printf("\n");

        float* PSNR_HVS_M = PSNRHVSM(original_image, image, width, height, channels);
        printf("PSNR-HVS / PSNR-HVS-M: ");
        for(int i = 0; i < 6; i+=2)
            printf(" %f / %f ", PSNR_HVS_M[i], PSNR_HVS_M[i+1]);
        printf("\n\n");
    }

    void RAWImage::transfer2OtherType(ImProcessing::ImagePixelType type) {
        if(ImType != type)
        {
            if(ImType == TYPE_BGR && type == TYPE_3ARRAY)
            {
                B = new uint8_t[width*height];
                G = new uint8_t[width*height];
                R = new uint8_t[width*height];
                int p = 0;
                for(int i = 0; i < width*height*channels; i+=channels, p++)
                {
                    B[p] = image[i];
                    G[p] = image[i+1];
                    R[p] = image[i+2];
                }
                delete image;
            }
            else if(ImType == TYPE_3ARRAY && type == TYPE_BGR)
            {
                image = new uint8_t[width*height*channels];
                int p = 0;
                for(int i = 0; i < width*height*channels; i+=channels, p++)
                {
                    image[i] = B[p];
                    image[i+1] = G[p];
                    image[i+2] = R[p];
                }
                delete B;
                delete G;
                delete R;
            }
            ImType = type;
        }
    }
}