#include "imageaccessor.h"
#include "grayscalematrix.h"
#include "convolution.h"
#include <QDir>
#include <QString>
#include <iostream>

/*
Лабораторная работа 1. Свертка
-Написать основу для представления изображений и их обработки свертками
-Реализовать вычисление частных производных и оператора Собеля
-Реализовать фильтр Гаусса
-Реализовать отображение полученных результатов
NB! Нормирование выходных данных
*/
void lab1(QString path, QString fileName, double sigma)
{

    std::cout<<"load img..."<<std::endl;
    GrayScaleMatrix origMatr = ImageAccessor::GetMatrixFromImage(path+fileName);

    QDir dir;
    path = path+"lab1"+fileName+"\\";
    dir.mkdir(path);

    std::cout<<"save orig grayscale img..."<<std::endl;
    ImageAccessor::DrawImageFromMatrix(origMatr,path+"orig_"+fileName);

    std::cout<<"compute derivateX..."<<std::endl;
    GrayScaleMatrix derivateX = Convolution::GetDerivateX(origMatr);
    std::cout<<"save derivateX img..."<<std::endl;
    ImageAccessor::DrawImageFromMatrix(derivateX,path+"derivateX_"+fileName);

    std::cout<<"compute derivateY..."<<std::endl;
    GrayScaleMatrix derivateY = Convolution::GetDerivateY(origMatr);
    std::cout<<"save derivateY img..."<<std::endl;
    ImageAccessor::DrawImageFromMatrix(derivateY,path+"derivateY_"+fileName);

    std::cout<<"compute Sobel operator..."<<std::endl;
    GrayScaleMatrix sobel = Convolution::SobelOperator(origMatr);
    std::cout<<"save sobel img..."<<std::endl;
    ImageAccessor::DrawImageFromMatrix(sobel,path+"sobel_"+fileName);

    std::cout<<"compute Gauss filter..."<<std::endl;
    GrayScaleMatrix gauss = Convolution::GaussianFilter(origMatr,sigma);
    std::cout<<"save gauss img..."<<std::endl;
    ImageAccessor::DrawImageFromMatrix(gauss,path+"gauss_"+fileName);

    std::cout<<"Done!"<<std::endl;
}
