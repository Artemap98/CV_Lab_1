#ifndef CONVOLUTION_H
#define CONVOLUTION_H
#include "grayscalematrix.h"

#include <QVector>
#include <QtMath>

class Convolution
{
public:
    static GrayScaleMatrix Convolute(GrayScaleMatrix inputGSMatrix, QVector<QVector<double>> convCore);
    static GrayScaleMatrix DerivateX(GrayScaleMatrix inputGSMatrix); //получить массивы с частными производными
    static GrayScaleMatrix DerivateY(GrayScaleMatrix inputGSMatrix);
    static GrayScaleMatrix SobelOperator(GrayScaleMatrix inputGSMatrix);
    static GrayScaleMatrix SobelOperator(GrayScaleMatrix derivateX,GrayScaleMatrix derivateY);
    static GrayScaleMatrix GradientDirection(GrayScaleMatrix derivateXMatrix, GrayScaleMatrix derivateYMatrix);

    static GrayScaleMatrix GaussianFilter(GrayScaleMatrix inputGSMatrix,double sigma);//Сепарабельный фильтр Гаусса

};

#endif // CONVOLUTION_H
