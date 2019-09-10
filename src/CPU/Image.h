#ifndef IMAGEPROCCESSING_IMAGE_H
#define IMAGEPROCCESSING_IMAGE_H

#include <png.h>

namespace ImProcessing{
    // Dimension of image
    class Dimension{
    public:
        Dimension(int width_, int height_, int depth_);
        int getWidth(){
            return width;
        }
        int getHeight(){
            return height;
        }
        int getDepth(){
            return depth;
        }

    private:
        int width;
        int height;
        int depth;
    };

    class Image{
    public:
        // Constructor for class
        Image(int *image_pixel_, int width_, int height_, int depth_);
        Dimension* getDimension();
        double* getBlockPixel(int pixel_w_start, int pixel_h_start, int width_block, int height_block);
        double* frequency_calculate(); // Calculate frequency coefficients

    private:
        Dimension *image_dimension;
        int *image_pixel;

    };
}


#endif //IMAGEPROCCESSING_IMAGE_H