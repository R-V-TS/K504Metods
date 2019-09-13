//
// Created by rostislav on 10.09.19.
//

#include <cstdlib>
#include "RAWImage.h"

int main()
{
    using namespace ImProcessing;
    RAWImage image("../Images/baboon.png", TYPE_BGR);

    image.ApplyDCT(8, 0.012, true);

    float* coeff = image.DCTCoefficients();
    for(int i = 0; i < 16; i++) printf("%f ", coeff[i]);
    image.show();
    return 0;
}


