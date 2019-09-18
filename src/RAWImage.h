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
        RAWImage(unsigned char *image_pixel, int image_width, int image_height, ImagePixelType type);
        RAWImage(std::string filename, ImagePixelType type);
        ~RAWImage();
        void show();
        void save(std::string destination_folder);
        void ApplyDCT(int window_size, double threshold, bool device); // device == true for CPU variant
        float* DCTCoefficients(bool device);

    private:
        // functions
        void createImageArray(unsigned char *image_pixel, int image_width, int image_height, ImagePixelType type);
        void transfer2OtherType(ImagePixelType type);

        // params
        unsigned char *image;
        unsigned char *B;
        unsigned char *G;
        unsigned char *R;
        int width;
        int height;
        int channels;
        ImagePixelType ImType;

    };
}


#endif //CUDA__RAWIMAGE_H
