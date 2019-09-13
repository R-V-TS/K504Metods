//
// Created by rostislav on 10.09.19.
//

#include <cstdlib>
#include <string>

#ifndef CUDA__RAWIMAGE_H
#define CUDA__RAWIMAGE_H

namespace ImProcessing {
    enum ImagePixelType {TYPE_BGR = 1, TYPE_RGB, TYPE_GRAYSCALE, TYPE_3ARRAY};

    class RAWImage {
    public:
        RAWImage(char *image_pixel, int image_width, int image_height, ImagePixelType type);
        RAWImage(std::string filename, ImagePixelType type);
        ~RAWImage();
        void show();
        void ApplyDCT(int window_size, double threshold, bool device); // device == true for CPU variant
        float* DCTCoefficients();

    private:
        // functions
        void createImageArray(char *image_pixel, int image_width, int image_height, ImagePixelType type);
        void transfer2OtherType(ImagePixelType type);

        // params
        char *image;
        char *B;
        char *G;
        char *R;
        int width;
        int height;
        int channels;
        ImagePixelType ImType;

    };
}


#endif //CUDA__RAWIMAGE_H
