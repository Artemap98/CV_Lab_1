#ifndef IMAGEACCESSOR_H
#define IMAGEACCESSOR_H
#include "descriptor.h"
#include "grayscalematrix.h"

#include <QImage>
#include <QVector>



class ImageAccessor
{
public:
    static GrayScaleMatrix GetMatrixFromImage(QString path);
    static void DrawImageFromMatrix(GrayScaleMatrix gsMatrix, QString path);
    static QImage GetImageFromMatrix(GrayScaleMatrix gsMatrix);
    static GrayScaleMatrix PrintDescriptorOnImage(GrayScaleMatrix inputGSMatrix, QVector<Descriptor> descriptors);
};

#endif // IMAGEACCESSOR_H
