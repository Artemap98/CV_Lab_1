#ifndef KEYPOINTS_H
#define KEYPOINTS_H

#include <QVector>
#include "grayscalematrix.h"
#include "convolution.h"
#include <QDebug>

class KeyFeatures
{
public:
    KeyFeatures();
    struct KeyPoint
    {
        int x,y;
        double sLocal;
        KeyPoint(int x, int y, double sLocal)
        {
            this->x = x;
            this->y = y;
            this->sLocal = sLocal;
        }
    };

    struct KeyPointSet
    {
        QVector<KeyPoint> keyPoints;
    };

    static GrayScaleMatrix GetMoravecMatrix(GrayScaleMatrix inputGSMatrix, int windowRadius);
    static GrayScaleMatrix GetPointsMoravec(GrayScaleMatrix inputGSMatrix, GrayScaleMatrix moravecMatrix, int windowRadius, int resultPointsNum);
    static GrayScaleMatrix GetHarrisMatrix(GrayScaleMatrix inputGSMatrix, int windowRadius);
    static GrayScaleMatrix GetPointsHarris(GrayScaleMatrix inputGSMatrix, GrayScaleMatrix moravecMatrix, int windowRadius, int resultPointsNum);
    static KeyPointSet GetLocalMaximums(GrayScaleMatrix harrisMatrix, int windowRadius, bool isHarris);
    static KeyPointSet ReducePoints(KeyPointSet points, int resultPointsNum, int maxRadius);


};

#endif // KEYPOINTS_H
