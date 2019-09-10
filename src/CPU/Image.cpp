#include "Image.h"

#include <iostream>

// create new namespace
namespace ImProcessing{

    // Constructor for class Dimension
    // @param width, height and depth
    Dimension::Dimension(int width_, int height_, int depth_) {
        width = width_;
        height = height_;
        depth = depth_;
    }

    //Constructor for class Image
    Image::Image(int *image_pixel_, int width_, int height_, int depth_) {
        image_dimension = new Dimension(width_, height_, depth_);
        image_pixel = image_pixel_;
    }

    Dimension* Image::getDimension() {
        return image_dimension;
    }

    double* Image::getBlockPixel(int pixel_w_start, int pixel_h_start, int width_block, int height_block) {
        double* block = new double[width_block*height_block];
        for(int i = 0; i < width_block; i++)
            for(int j = 0; j < height_block; j++)
            {
                *(block+((i*8)+j)) = *(image_pixel+(((i+pixel_w_start)*image_dimension->getWidth())+(j+pixel_h_start)));
            }
        return block;
    }

    double* Image::frequency_calculate() {

    }
}