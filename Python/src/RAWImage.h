//
// Created by rostislav on 10.09.19.
//

#include <cstdlib>
#include <vector>
#include <cstdint>

#ifndef CUDA__RAWIMAGE_H
#define CUDA__RAWIMAGE_H
using namespace std;

namespace ImProcessing {
    enum ImagePixelType {TYPE_BGR = 1, TYPE_RGB, TYPE_GRAYSCALE, TYPE_3ARRAY};

    class RAWImage {
    public:
        RAWImage(uint8_t* pixel_array, int image_width, int image_height, int channels_, ImagePixelType type);
        ~RAWImage();
        void ApplyDCT(int window_size, double threshold);
        void DCTCoefficients(float* DCT_Coff, int n);
        void getImage(uint8_t* pixel_array, int n);
        void printImage();
        void AddNoise(float mu, float sigma);
        void printImageCharacteristics();

    private:
        // functions
        void transfer2OtherType(ImagePixelType type);

        // params
        uint8_t *original_image;
        uint8_t *image;
        uint8_t *B;
        uint8_t *G;
        uint8_t *R;
        int width;
        int height;
        int channels;
        ImagePixelType ImType;

    };
}


#endif //CUDA__RAWIMAGE_H
