#include "grayscalematrix.h"
#include "convolution.h"

class ScaleOperation
{
public:
    ScaleOperation();
    static GrayScaleMatrix Downsample(GrayScaleMatrix inputMatrix);

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


    static Pyramid GetPyramid(GrayScaleMatrix inputMatrix, int octaveNum, int levelNum, double sigmaA, double sigma0); 
    static double GetL(ScaleOperation::Pyramid inputPyramid, int y, int x, double sigma);

};

