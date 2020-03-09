#include "grayscalematrix.h"

GrayScaleMatrix::GrayScaleMatrix(int width, int height)
{
    QVector<QVector<double>> tempMatrix;
    for(int i=0; i < height; i++)
    {
        QVector<double> tempLine;
        for(int j = 0; j < width; j++)
        {
            tempLine.append(0);
        }
        tempMatrix.append(tempLine);
    }
    imageMatrix = tempMatrix;
}

GrayScaleMatrix::GrayScaleMatrix(QVector<QVector<double>> newMatrix)
{
    imageMatrix = newMatrix;
}

GrayScaleMatrix::GrayScaleMatrix(QVector<QVector<unsigned char>> newMatrix)
{
    SetMatrixDoubleFrom255(newMatrix);
}

QVector<QVector<unsigned char>> GrayScaleMatrix::GetMatrix255()
{
    QVector<QVector<unsigned char>> resultMatrix;
    double width = imageMatrix[0].size();
    double height = imageMatrix.size();
    for (int i = 0; i < height; i++) {
        QVector<unsigned char> resiltLine;
        for (int j = 0; j < width; j++) {
            resiltLine.append(imageMatrix[j][i]*255);
        }
        resultMatrix.append(resiltLine);
    }
    return resultMatrix;
}


QVector<QVector<double>> GrayScaleMatrix::GetMatrixDouble()
{
    return imageMatrix;
}


void GrayScaleMatrix::SetMatrixDoubleFrom255(QVector<QVector<unsigned char>> inputMatrix)
{
    QVector<QVector<double>> resultMatrix;
    double width = inputMatrix[0].size();
    double height = inputMatrix.size();
    for (int i = 0; i < height; i++) {
        QVector<double> resiltLine;
        for (int j = 0; j < width; j++) {
            resiltLine.append(inputMatrix[j][i]/255);
        }
        resultMatrix.append(resiltLine);
    }
    imageMatrix = resultMatrix;
}


double GrayScaleMatrix::GetValue(int x, int y)
{
    return imageMatrix.at(y).at(x);
}

void GrayScaleMatrix::SetValue(int x, int y, double value)
{
    imageMatrix[y][x] = value;
}

