#ifndef GAUSSPIRAMYDE_H
#define GAUSSPIRAMYDE_H

#include "grayscalematrix.h"
#include "convolution.h"



class ScaleOperation
{
public:
    ScaleOperation();

    struct Layer{
        GrayScaleMatrix matrix;
        double currentSigma;
        double actualSigma;
        Layer(GrayScaleMatrix inputmatrix,double currSigma, double actSigma):matrix(inputmatrix)
        {
            currentSigma = currSigma;
            actualSigma = actSigma;
        }
    };

    struct Octave{
        QVector<Layer> layers;
        int level;
    };

    struct Pyramid{
        QVector<Octave> octaves;
    };

    static GrayScaleMatrix Downsample(GrayScaleMatrix inputMatrix);
    static Pyramid GetPyramid(GrayScaleMatrix inputMatrix, int octaveNum, int levelNum, double sigmaA, double sigma0);
    static double GetL(ScaleOperation::Pyramid inputPyramid, int y, int x, double sigma);

};

#endif // GAUSSPIRAMYDE_H
