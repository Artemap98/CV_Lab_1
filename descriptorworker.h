#ifndef DESCRIPTORWORKER_H
#define DESCRIPTORWORKER_H
#include "descriptor.h"
#include "grayscalematrix.h"

class DescriptorWorker
{
public:
    DescriptorWorker();
    //получить дескриптор изображения(для лр 4)
    static QVector<Descriptor> GetDescriptorsFromImage(GrayScaleMatrix inputGSMatrix,
                                                       int harrisRadius,
                                                       int harrisPointsNum,
                                                       int basketNum,
                                                       int histogramGridSize,
                                                       int descriptorSize);

    //получить дескриптор, стойкий к наклонам
    static QVector<Descriptor> GetDescriptorsWithRotation( GrayScaleMatrix inputGSMatrix,
                                                           int harrisRadius,
                                                           int harrisPointsNum,
                                                           int basketNum,
                                                           int histogramGridSize,
                                                           int descriptorSize);

    //определить угол каждой интересной точки
    static KeyFeatures::KeyPointSet OrientPoints(KeyFeatures::KeyPointSet inputPoints,
                                                 GrayScaleMatrix gradientDirection,
                                                 GrayScaleMatrix gradientMagnitude,
                                                 int histogramGridSize,
                                                 int descriptorSize);

    //сравнить два дескриптора и вернуть объединенную картинку
    static QImage GetTwoImageDescriptorComparsion(
            GrayScaleMatrix inputGSMatrix1, QVector<Descriptor> descriptors1,
            GrayScaleMatrix inputGSMatrix2, QVector<Descriptor> descriptors2);

    //склеить две матрицы
    static GrayScaleMatrix MergeTwoMatrix(GrayScaleMatrix inputGSMatrix1, GrayScaleMatrix inputGSMatrix2);

    static double GetDistanceBetweenDescriptors(Descriptor d1,Descriptor d2);
    static double GetDistanceBetweenPoints(int x1, int x2, int y1, int y2);

};

#endif // DESCRIPTORWORKER_H
