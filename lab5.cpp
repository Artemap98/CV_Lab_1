//Реализовать относительную инвариантность вычисления дескрипторов
//    к вращению изображений на основе подхода SIFT.
//Реализовать этап оценки ориентации интересной точки и поворота сетки,
//    в которой вычисляются гистограммы градиентов.
//Оценить полученный алгоритм с точки зрения реакции на соответствующие искажения изображений,
//сравнить с полученным в четвертой работе.


#include "descriptorworker.h"
#include "imageaccessor.h"
#include <iostream>
#include <QDir>

GrayScaleMatrix PrintDescriptorOnImage(GrayScaleMatrix inputGSMatrix, QVector<Descriptor> descriptors)
{
    foreach(Descriptor descr, descriptors)
    {
        for(int i=-1; i<=1; i++)
        {
            for(int j=-1; j<=1; j++)
            {
                if(i==0 || j==0)
                {
                    try
                    {
                        inputGSMatrix.SetValue(descr.GetX()+j,descr.GetY()+i,1);
                    } catch(_exception e){}
                }
            }
        }
    }
    return  inputGSMatrix;
}

void lab5(QString path,
          QString fileName1, QString extension1,
          QString fileName2, QString extension2,
          int harrisRadius,
          int harrisPointsNum,
          int basketNum,
          int histogramGridSize,
          int descriptorSize)
{
    QDir dir;
    QString labPath = path+"lab5"+fileName1+fileName2+"\\";
    dir.mkdir(labPath);
    std::cout<<"load img1..."<<std::endl;
    GrayScaleMatrix inputMatrix1 = ImageAccessor::GetMatrixFromImage(path+fileName1+extension1);
    std::cout<<"load img2..."<<std::endl;
    GrayScaleMatrix inputMatrix2 = ImageAccessor::GetMatrixFromImage(path+fileName2+extension2);

    std::cout<<"compute descriptor1..."<<std::endl;
    QVector<Descriptor> descriptors1 = DescriptorWorker::GetDescriptorsWithRotation(inputMatrix1, harrisRadius, harrisPointsNum, basketNum, histogramGridSize, descriptorSize);
    GrayScaleMatrix descrImage1 = PrintDescriptorOnImage(inputMatrix1,descriptors1);
    ImageAccessor::DrawImageFromMatrix(descrImage1,labPath+fileName1+"points"+extension1);

    std::cout<<"compute descriptor2..."<<std::endl;
    QVector<Descriptor> descriptors2 = DescriptorWorker::GetDescriptorsWithRotation(inputMatrix2, harrisRadius, harrisPointsNum, basketNum, histogramGridSize, descriptorSize);
    GrayScaleMatrix descrImage2 = PrintDescriptorOnImage(inputMatrix2,descriptors2);
    ImageAccessor::DrawImageFromMatrix(descrImage2,labPath+fileName2+"points"+extension2);

    QImage result = DescriptorWorker::GetTwoImageDescriptorComparsion(inputMatrix1,descriptors1,inputMatrix2,descriptors2);
    result.save(labPath+"result.png");
    std::cout<<"done"<<std::endl;
}


