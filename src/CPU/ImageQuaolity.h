//
// Created by rostislav on 24.08.19.
//

#ifndef IMAGEPROCESSING_EXE_IMAGEQUAOLITY_H
#define IMAGEPROCESSING_EXE_IMAGEQUAOLITY_H

namespace ImProcessing{
    double mean(double *image_, int width, int height);
    double SMD(double *image_, int width, int height); // square mean deviation
}

#endif //IMAGEPROCESSING_EXE_IMAGEQUAOLITY_H
