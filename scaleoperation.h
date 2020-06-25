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
        int overlay = 0;
    };


    static Pyramid GetPyramid(GrayScaleMatrix inputMatrix,
                              int numOfOctaves, int numOfLayers,
                              double sigmaA, double sigma0
                              );
    static Pyramid GetPyramidWithOverlay(GrayScaleMatrix inputMatrix,
                                         int numOfOctaves, int numOfLayers,
                                         double sigmaA, double sigma0,
                                         int overlay = 3
                                         ) ;
    static Pyramid GetDoGPyramid(Pyramid inputPyramid);
    static double GetL(ScaleOperation::Pyramid inputPyramid, int y, int x, double sigma);

};

