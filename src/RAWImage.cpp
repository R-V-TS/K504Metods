//
// Created by rostislav on 10.09.19.
//

#include "RAWImage.h"
#include "ERROR_Printer.h"
#include <opencv2/opencv.hpp>
#include <cmath>
#include <vector>

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
            char *image_pxl = (char*)imCV.data;
            createImageArray(image_pxl, imCV.cols, imCV.rows, type);
        }

    }

    RAWImage::RAWImage(char *image_pixel, int image_width, int image_height, ImProcessing::ImagePixelType type) {
        createImageArray(image_pixel, image_width, image_height, type);
    }

    RAWImage::~RAWImage() {
        //delete image;
        width = 0;
        height = 0;
    }

    void RAWImage::createImageArray(char *image_pixel, int image_width, int image_height,
                                    ImProcessing::ImagePixelType type) {
        switch(type)
        {
            case 1:
                image = new char[image_height*image_width*3];
                for(int i = 0; i < image_width*image_height*3; i++){
                    image[i] = image_pixel[i];
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
                B = new char[width*height];
                G = new char[width*height];
                R = new char[width*height];
                int p = 0;
                for(int i = 0; i < image_width*image_height*channels; i+=channels)
                {
                    B[p] = image_pixel[i];
                    G[p] = image_pixel[i+1];
                    R[p] = image_pixel[i+2];
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

        float slices_sum[4] = {0,0,0,0};
        float mean[4] = {0,0,0,0};
        float variance[4] = {0,0,0,0};
        float skeweness[4] = {0,0,0,0};
        float kurtosis[4] = {0,0,0,0};
        int count = 0;
        float fulldct = 0;
        int pixel_p = 0;
        for(int i = 0; i < height; i+=window_size)
        {
            for(int j = 0; j < width; j+=window_size)
            {
                getImageBlock(DCT_ARRAY, i, j, width, window_size, block);
                block[0] = 0;

                for(int l = 0; l < window_size; l++) {
                    for (int k = 0; k < window_size; k++) {
                        fulldct += block[(l * window_size) + k];
                        slices_sum[matrix_flag[l][k]-1] += block[(l * window_size) + k];
                    }
                }

                for(int l = 0; l < 4; l++) {
                    slices_sum[l] /= fulldct;
                    EST[pixel_p] = slices_sum[l];
                    pixel_p++;
                }

                count++;
            }
        }

        for(int l = 0; l < 4; l++)
        {
            pixel_p = 0;
            for(int i = 0; i < height/8; i++)
            {
                for(int j = 0; j < width/8; j++)
                {
                    mean[l] += EST[pixel_p];
                    pixel_p += 4;
                }
            }

            mean[l] /= count;

            pixel_p = 0;
            for(int i = 0; i < height/8; i++)
            {
                for(int j = 0; j < width/8; j++)
                {
                    variance[l] += pow((EST[pixel_p] - mean[l]), 2);
                    skeweness[l] += pow((EST[pixel_p]-mean[l]), 3);
                    kurtosis[l] += pow((EST[pixel_p]-mean[l]), 4);
                    pixel_p += 4;
                }
            }

            skeweness[l] = (skeweness[l]/count)/pow(sqrt(variance[l]/count), 3);
            kurtosis[l] = (kurtosis[l]/count)/pow(variance[l]/count, 2);
            variance[l] = variance[l]/(count-1);
            coefficients[z] = mean[l];
            coefficients[z+1] = variance[l];
            coefficients[z+2] = skeweness[l];
            coefficients[z+3] = kurtosis[l];
            z += 4;
        }

        return coefficients;
    }

    void RAWImage::transfer2OtherType(ImProcessing::ImagePixelType type) {
        if(ImType != type)
        {
            if(ImType == TYPE_BGR && type == TYPE_3ARRAY)
            {
                B = new char[width*height];
                G = new char[width*height];
                R = new char[width*height];
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
                image = new char[width*height*channels];
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