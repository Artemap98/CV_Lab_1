#include "imageaccessor.h"
#include "grayscalematrix.h"
#include "convolution.h"
#include <QString>
#include <iostream>


void lab1(QString path, QString fileName, double sigma)
{
    std::cout<<"load img..."<<std::endl;
    GrayScaleMatrix origMatr = ImageAccessor::GetMatrixFromImage(path+fileName);

    std::cout<<"compute derivateX..."<<std::endl;
    GrayScaleMatrix derivateX = Convolution::GetDerivateX(origMatr);
    std::cout<<"compute derivateY..."<<std::endl;
    GrayScaleMatrix derivateY = Convolution::GetDerivateY(origMatr);
    std::cout<<"compute Sobel operator..."<<std::endl;
    GrayScaleMatrix sobel = Convolution::SobelOperator(origMatr);
    std::cout<<"compute Gauss filter..."<<std::endl;
    GrayScaleMatrix gauss = Convolution::GaussianFilter(origMatr,sigma);

    std::cout<<"save orig grayscale img..."<<std::endl;
    ImageAccessor::DrawImageFromMatrix(origMatr,path+"orig_"+fileName);
    std::cout<<"save derivateX img..."<<std::endl;
    ImageAccessor::DrawImageFromMatrix(derivateX,path+"derivateX_"+fileName);
    std::cout<<"save derivateY img..."<<std::endl;
    ImageAccessor::DrawImageFromMatrix(derivateY,path+"derivateY_"+fileName);
    std::cout<<"save sobel img..."<<std::endl;
    ImageAccessor::DrawImageFromMatrix(sobel,path+"sobel_"+fileName);
    std::cout<<"save gauss img..."<<std::endl;
    ImageAccessor::DrawImageFromMatrix(gauss,path+"gauss_"+fileName);
    std::cout<<"Done!"<<std::endl;
}
