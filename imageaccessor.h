#ifndef IMAGEACCESSOR_H
#define IMAGEACCESSOR_H
#include "grayscalematrix.h"

#include <QImage>
#include <QVector>



class ImageAccessor
{
public:
    ImageAccessor();
    static GrayScaleMatrix GetMatrixFromImage(QString path);
    static void DrawImageFromMatrix(GrayScaleMatrix imageMatrix, QString path);
};

#endif // IMAGEACCESSOR_H
