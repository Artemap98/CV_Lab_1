#ifndef IMAGEACCESSOR_H
#define IMAGEACCESSOR_H
#include "grayscalematrix.h"

#include <QImage>
#include <QVector>



class ImageAccessor
{
public:
    static GrayScaleMatrix GetMatrixFromImage(QString path);
    static void DrawImageFromMatrix(GrayScaleMatrix gsMatrix, QString path);
    static QImage GetImageFromMatrix(GrayScaleMatrix gsMatrix);
};

#endif // IMAGEACCESSOR_H
