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
        double angle = 0;
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

    //получить отклик Моравека
    static GrayScaleMatrix GetMoravecMatrix(GrayScaleMatrix inputGSMatrix, int windowRadius);
    //отобразить найденные интересные точки на изображении на основе отклика
    static GrayScaleMatrix GetResultMoravec(GrayScaleMatrix inputGSMatrix, GrayScaleMatrix moravecMatrix, int windowRadius, int resultPointsNum);

    //отклик Харриса
    static GrayScaleMatrix GetHarrisMatrix(GrayScaleMatrix inputGSMatrix, int windowRadius);
    //получить набор интересных точек на основе отклика
    static KeyPointSet GetPointsHarris(GrayScaleMatrix inputGSMatrix, GrayScaleMatrix moravecMatrix, int windowRadius, int resultPointsNum);
    //отобразить точки на исходной картинке
    static GrayScaleMatrix GetResultHarrisMatrix(GrayScaleMatrix inputGSMatrix, KeyPointSet interestingPoints);

    static KeyPointSet GetLocalMaximums(GrayScaleMatrix harrisMatrix, int windowRadius, bool isHarris);
    static KeyPointSet ReducePoints(KeyPointSet points, int resultPointsNum, int maxRadius);


};

#endif // KEYPOINTS_H
