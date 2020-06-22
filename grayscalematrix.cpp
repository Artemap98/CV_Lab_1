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

    double min=1000000,max=-1000000;

    //output data normalization
    for(int i=0; i<imageMatrix.size(); i++)
    {
        for(int j=0; j<imageMatrix[0].size(); j++)
        {
            if(imageMatrix[i][j]>max)
                max=imageMatrix[i][j];
            if(imageMatrix[i][j]<min)
                min=imageMatrix[i][j];
        }
    }


    for (int i = 0; i < height; i++) {
        QVector<unsigned char> resultLine;
        for (int j = 0; j < width; j++) {
            resultLine.append((double)(imageMatrix[i][j]-min)*255/(max-min));
        }
        resultMatrix.append(resultLine);
    }
    return resultMatrix;
}


void GrayScaleMatrix::NormalizeDouble()
{
    int width = imageMatrix[0].size();
    int height = imageMatrix.size();

    double min=1,max=0;

    for(int i=0; i<imageMatrix.size(); i++)
    {
        for(int j=0; j<imageMatrix[0].size(); j++)
        {
            if(imageMatrix[i][j]>max)
                max=imageMatrix[i][j];
            if(imageMatrix[i][j]<min)
                min=imageMatrix[i][j];
        }
    }

    //нормализация выходных данных
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            imageMatrix[i][j] = (imageMatrix[i][j]-min)/(max-min);
        }
    }
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
        QVector<double> resultLine;
        for (int j = 0; j < width; j++) {
            resultLine.append((double)inputMatrix[i][j]/255);
        }
        resultMatrix.append(resultLine);
    }
    imageMatrix = resultMatrix;
}


double GrayScaleMatrix::GetValue(int x, int y)
{
    int i, j;
    i = x <= 0 ? 0 : (x >= imageMatrix[0].size() - 1 ? imageMatrix[0].size() - 1 : x);
    j = y <= 0 ? 0 : (y >= imageMatrix.size() - 1 ? imageMatrix.size() - 1 : y);
    return imageMatrix[j][i];
}

void GrayScaleMatrix::SetValue(int x, int y, double value)
{
    if(x < 0 || y < 0 || x >= imageMatrix[0].size() || y >= imageMatrix.size())
    {
        return;
    }
    imageMatrix[y][x] = value;
}

int GrayScaleMatrix::GetWidth()
{
    return imageMatrix[0].size();
}

int GrayScaleMatrix::GetHeight()
{
    return imageMatrix.size();
}

void GrayScaleMatrix::SetMatrixDouble(QVector<QVector<double>> inputMatrix)
{
    imageMatrix = inputMatrix;
}

