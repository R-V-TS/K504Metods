//
// Created by rostislav on 20.07.19.
//

#ifndef IMAGEPROCCESSING_MATRIX_H
#define IMAGEPROCCESSING_MATRIX_H

namespace ImProcessing {
    class Matrix{
    public:
        Matrix(int width_, int height_);
        double* matrixMultiply(double* matrix1, double* matrix2);
        void matrixPower(double* matrix);
        double matrixFullSum(double* matrix);

    private:
        // Parameters of matrix size
        int width;
        int height;

    };
}


#endif //IMAGEPROCCESSING_MATRIX_H
