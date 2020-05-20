#ifndef DESCRIPTORWORKER_H
#define DESCRIPTORWORKER_H
#include "descriptor.h"
#include "grayscalematrix.h"

class DescriptorWorker
{
public:
    DescriptorWorker();
    static QVector<Descriptor> GetDescriptorsFromImage(GrayScaleMatrix inputGSMatrix, 
                                                       int harrisRadius, 
                                                       int harrisPointsNum,
                                                       int basketNum, 
                                                       int histogramGridSize,
                                                       int descriptorSize);


    static QImage GetTwoImageDescriptorComparsion(
            GrayScaleMatrix inputGSMatrix1, QVector<Descriptor> descriptors1,
            GrayScaleMatrix inputGSMatrix2, QVector<Descriptor> descriptors2);

    static GrayScaleMatrix MergeTwoMatrix(GrayScaleMatrix inputGSMatrix1, GrayScaleMatrix inputGSMatrix2);

    static double GetDistance(Descriptor d1,Descriptor d2);

};

#endif // DESCRIPTORWORKER_H
