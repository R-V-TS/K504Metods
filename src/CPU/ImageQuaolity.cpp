//
// Created by rostislav on 24.08.19.
//

#include "ImageQuaolity.h"
#include "cmath"

namespace ImProcessing{
    double mean(double *image_, int width, int height){
        double sum = 0;
        for(int i = 0; i < height; i++)
            for(int j = 0; j < width; j++)
            {
                sum += (double)(*(image_+(width*i)+j));
            }
        return sum/(width*height);
    }

    double SMD(double *image_, int width, int height)
    {
        double sum = 0;
        double mean_ = mean(image_, width, height);
        for(int i = 0; i < height; i++)
            for(int j = 0; j < width; j++)
            {
                sum += ((*(image_)-mean_)*(*(image_)-mean_));
            }
        return sqrt(sum/(width*height));
    }
}
