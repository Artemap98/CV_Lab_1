//◦ Реализовать вычисление дескрипторов окрестностей заданных точек путем вычисления
//  градиентов в каждой точки изображения и разбиения окрестности на сетку
//◦ Реализовать вычисление гистограмм градиентов в ячейках сетки и нормализацию полученных дескрипторов
//◦ Реализовать визуализацию результатов поиска ближайших дескрипторов в двух изображениях.

#include "descriptorworker.h"
#include "imageaccessor.h"
#include <iostream>
#include <QDir>

void lab4(QString path,
          QString fileName1, QString extension1,
          QString fileName2, QString extension2,
          int harrisRadius,
          int harrisPointsNum,
          int basketNum,
          int histogramGridSize,
          int descriptorSize)
{
    QDir dir;
    QString labPath = path+"lab4"+fileName1+fileName2+"\\";
    dir.mkdir(labPath);
    std::cout<<"load img1..."<<std::endl;
    GrayScaleMatrix inputMatrix1 = ImageAccessor::GetMatrixFromImage(path+fileName1+extension1);
    std::cout<<"load img2..."<<std::endl;
    GrayScaleMatrix inputMatrix2 = ImageAccessor::GetMatrixFromImage(path+fileName2+extension2);
    std::cout<<"compute descriptor1..."<<std::endl;
    QVector<Descriptor> descriptors1 = DescriptorWorker::GetDescriptorsFromImage(inputMatrix1, harrisRadius, harrisPointsNum, basketNum, histogramGridSize, descriptorSize);
    std::cout<<"compute descriptor2..."<<std::endl;
    QVector<Descriptor> descriptors2 = DescriptorWorker::GetDescriptorsFromImage(inputMatrix2, harrisRadius, harrisPointsNum, basketNum, histogramGridSize, descriptorSize);
    QImage result = DescriptorWorker::GetTwoImageDescriptorComparsion(inputMatrix1,descriptors1,inputMatrix2,descriptors2);
    result.save(labPath+"result.png");
    std::cout<<"done"<<std::endl;
}
