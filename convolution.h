#ifndef CONVOLUTION_H
#define CONVOLUTION_H
#include "grayscalematrix.h"

#include <QVector>
#include <QtMath>

class Convolution
{
public:
    static GrayScaleMatrix Convolute(GrayScaleMatrix inputGSMatrix, QVector<QVector<double>> convCore);
    static QVector<QVector<double>> MatrixMult(QVector<QVector<double>> matrix1, QVector<QVector<double>> matrix2);

    static GrayScaleMatrix GetDerivateX(GrayScaleMatrix inputGSMatrix); //получить массивы с частными производными
    static GrayScaleMatrix GetDerivateY(GrayScaleMatrix inputGSMatrix);
    static GrayScaleMatrix SobelOperator(GrayScaleMatrix inputGSMatrix);
    static GrayScaleMatrix GaussianFilter(GrayScaleMatrix inputGSMatrix,double sigma);
};

#endif // CONVOLUTION_H
