//Реализовать распределение значений градиентов по смежным гистограммам,
//добавить весовые коэффициенты исходя из расстояния до соответствующих центров.


#include "descriptorworker.h"
#include "imageaccessor.h"
#include <iostream>
#include <QDir>

void lab7(QString path,
          QString fileName1, QString extension1,
          QString fileName2, QString extension2,
          int harrisRadius,
          int harrisPointsNum,
          int basketNum,
          int histogramGridSize,
          int descriptorSize)
{
    QDir dir;
    QString labPath = path+"lab6"+fileName1+fileName2+"\\";
    dir.mkdir(labPath);
    std::cout<<"load img1..."<<std::endl;
    GrayScaleMatrix inputMatrix1 = ImageAccessor::GetMatrixFromImage(path+fileName1+extension1);
    std::cout<<"load img2..."<<std::endl;
    GrayScaleMatrix inputMatrix2 = ImageAccessor::GetMatrixFromImage(path+fileName2+extension2);

    std::cout<<"compute descriptor1..."<<std::endl;
    QVector<Descriptor> descriptors1 = DescriptorWorker::GetDescriptorsBlob(inputMatrix1, harrisRadius, harrisPointsNum, basketNum, histogramGridSize, descriptorSize);
    GrayScaleMatrix descrImage1 = ImageAccessor::PrintDescriptorOnImage(inputMatrix1,descriptors1);
    ImageAccessor::DrawImageFromMatrix(descrImage1,labPath+fileName1+"points"+extension1);

    std::cout<<"compute descriptor2..."<<std::endl;
    QVector<Descriptor> descriptors2 = DescriptorWorker::GetDescriptorsBlob(inputMatrix2, harrisRadius, harrisPointsNum, basketNum, histogramGridSize, descriptorSize);
    GrayScaleMatrix descrImage2 = ImageAccessor::PrintDescriptorOnImage(inputMatrix2,descriptors2);
    ImageAccessor::DrawImageFromMatrix(descrImage2,labPath+fileName2+"points"+extension2);

    QImage result = DescriptorWorker::GetTwoImageDescriptorComparsion(inputMatrix1,descriptors1,inputMatrix2,descriptors2);
    result.save(labPath+"result.png");
    std::cout<<"done"<<std::endl;
}


