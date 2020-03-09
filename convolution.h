#ifndef CONVOLUTION_H
#define CONVOLUTION_H
#include "grayscalematrix.h"

#include <QVector>

class Convolution
{
public:
    Convolution();
    static GrayScaleMatrix Convolute(GrayScaleMatrix inputGSMatrix, QVector<QVector<double>> convCore);
    QVector<QVector<double>> MatrixMult(QVector<QVector<double>> matrix1, QVector<QVector<double>> matrix2);
};

#endif // CONVOLUTION_H
