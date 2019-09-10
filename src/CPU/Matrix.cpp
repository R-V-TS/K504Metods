//
// Created by rostislav on 20.07.19.
//

#include "Matrix.h"

namespace ImProcessing{
    Matrix::Matrix(int width_, int height_) {
        width = width_;
        height = height_;
    }

    double* Matrix::matrixMultiply(double* matrix1, double* matrix2) {
        double* multi = new double[width*height];
        for(int i = 0; i < width; i++)
        {
            for(int j = 0; j < height; j++)
            {
                double summ = 0;
                for(int k = 0; k< width; k++)
                {
                    summ += (*(matrix1+((i*width)+k)) * (*(matrix2+((k*width)+j))));
                }
                *(multi+((i*width)+j)) = summ;
            }
        }
        return multi;
    }

    void Matrix::matrixPower(double* matrix){
        for(int i = 0; i < width; i++)
            for (int j = 0; j < height; ++j) {
                *(matrix+(i*width)+j) = (*(matrix+(i*width)+j)) * (*(matrix+(i*width)+j));
            }
    }

    double Matrix::matrixFullSum(double *matrix) {
        double sum = 0;
        for(int i = 0; i < width; i++)
            for (int j = 0; j < height; ++j) {
                sum += *(matrix+(i*width)+j);
            }
        return sum;
    }
}
