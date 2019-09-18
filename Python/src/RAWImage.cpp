//
// Created by rostislav on 10.09.19.
//

#include "RAWImage.h"
#include <cmath>
#include <cstdio>


// Metods
#include "DCT.h"

using namespace std;

namespace ImProcessing
{
    RAWImage::RAWImage(short* pixel_array, int image_width, int image_height, int channels_, ImagePixelType type) {
        ImType = type;
        switch(ImType) {
            case TYPE_BGR: {
                int t = 0;
                image = new short[image_height * image_width * channels_];
                for (int i = 0; i < image_height; i++) {
                    for (int j = 0; j < image_width * channels_; j++) {
                        image[t] = pixel_array[(i*image_width*channels_) + j];
                        t++;
                    }
                }
                width = image_width;
                height = image_height;
                channels = channels_;
                //delete image_pixel
                break;
            }
            case 4: {
                width = image_width;
                height = image_height;
                channels = channels_;
                B = new short[width * height];
                G = new short[width * height];
                R = new short[width * height];
                int p = 0;
                for (int i = 0; i < image_height; i++) {
                    for (int j = 0; j < image_width * channels_; j += channels_) {
                        B[p] = pixel_array[(i*image_width*channels_) + j];
                        G[p] = pixel_array[(i*image_width*channels_) + j + 1];
                        R[p] = pixel_array[(i*image_width*channels_) + j + 2];
                        p++;
                    }
                }
                //free(image_pixel);
                break;
            }
        }
    }
    
    RAWImage::~RAWImage() {
        //delete image;
        width = 0;
        height = 0;
    }

    void RAWImage::printImage(){
        if(ImType != TYPE_BGR)
        {
            transfer2OtherType(TYPE_BGR);
        }
        for(int i = 0; i < height; i++)
        {
            for(int j = 0; j < width*channels; j++)
            {
                printf("%i ", image[(i*width*channels) + j]);
            }
            printf("\n");
        }
    };
    
    void RAWImage::getImage(short* pixel_array, int n){
        if(ImType != TYPE_BGR)
        {
            transfer2OtherType(TYPE_BGR);
        }
        for(int i = 0; i < width*height*channels; i++)
        {
            if(image[i] <= 255 && image[i] >=0)
                pixel_array[i] = image[i];
            else if(image[i] > 255)
                pixel_array[i] = 255;
            else
                pixel_array[i] = 0;
        }
    }

    void RAWImage::ApplyDCT(int window_size, double threshold) {
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
    }

    void RAWImage::DCTCoefficients(float* DCT_Coff, int n) {
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
        float* block = new float[window_size*window_size];
        float* DCT_ARRAY;
        DCT_ARRAY = DCT(B, width, height, window_size);
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
            pixel_p = l;
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
            DCT_Coff[z] = mean[l];
            DCT_Coff[z+1] = variance[l];
            DCT_Coff[z+2] = skeweness[l];
            DCT_Coff[z+3] = kurtosis[l];
            z += 4;
        }
    }

    void RAWImage::transfer2OtherType(ImProcessing::ImagePixelType type) {
        if(ImType != type)
        {
            if(ImType == TYPE_BGR && type == TYPE_3ARRAY)
            {
                B = new short[width*height];
                G = new short[width*height];
                R = new short[width*height];
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
                image = new short[width*height*channels];
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
