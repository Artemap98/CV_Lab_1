#ifndef IMAGEMATRIX_H
#define IMAGEMATRIX_H
#include <QVector>


class GrayScaleMatrix
{
public:
    GrayScaleMatrix(int width, int height);
    GrayScaleMatrix(QVector<QVector<double>> newMatrix);
    GrayScaleMatrix(QVector<QVector<unsigned char>> newMatrix);
    void NormalizeDouble();
    QVector<QVector<unsigned char>> GetMatrix255();
    QVector<QVector<double>> GetMatrixDouble();

    void SetMatrixDoubleFrom255(QVector<QVector<unsigned char>> inputMatrix);
    void SetMatrixDouble(QVector<QVector<double>> inputMatrix);

    double GetValue(int x, int y);
    void SetValue(int x, int y, double value);

    int GetWidth();
    int GetHeight();


private:
    QVector<QVector<double>> imageMatrix;
};

#endif // IMAGEMATRIX_H
