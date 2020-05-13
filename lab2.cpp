//Из заданного изображения построить гауссову пирамиду
//    ◦ Устанавливается количество октав
//    ◦ (можно вычислять исходя из размера изображения)
//    ◦ Устанавливается число уровней в октаве
//    ◦ Устанавливается 𝜎𝜎0и 𝜎𝜎1
//Реализовать функциюL(x,y,𝜎𝜎)
//    ◦ Поиск ближайшего изображения
//    ◦ Преобразование координат
//Реализовать отображение результатов
//    ◦ 𝜎𝜎на каждом масштабе в октаве
//    ◦ эффективная 𝜎𝜎для каждого масштаба
//    (Пов. сложность) реализовать построение пирамиды с -1 октавы

#include "imageaccessor.h"
#include "scaleoperation.h"

#include <QDir>
#include <QString>
#include <iostream>

void lab2(QString path, QString fileName, QString extension, int numOfOctaves, int numOfLayers, double sigmaA, double sigma0, double sigmaL)
{
    QDir dir;
    QString pyramidPath = path+"lab2"+fileName+"\\";
    dir.mkdir(pyramidPath);
    std::cout<<"load img..."<<std::endl;
    GrayScaleMatrix inputmatrix = ImageAccessor::GetMatrixFromImage(path+fileName+extension);
    std::cout<<"computing pyramid..."<<std::endl;
    ScaleOperation::Pyramid myPyramid = ScaleOperation::GetPyramid(inputmatrix,numOfOctaves, numOfLayers, sigmaA, sigma0);

    std::cout<<"saving pyramid..."<<std::endl;
    int octaveCount=0, layerCount=0;
    foreach(ScaleOperation::Octave currOctave, myPyramid.octaves)
    {
        layerCount=0;
        foreach(ScaleOperation::Layer currLayer, currOctave.layers)
        {
            ImageAccessor::DrawImageFromMatrix(currLayer.matrix,pyramidPath+"octave"+QString::number(octaveCount)+"_layer"+QString::number(layerCount++)+extension);
        }
        octaveCount++;
    }

    std::cout<<"Computing image with L(x,y,sigma)..."<<std::endl;
    GrayScaleMatrix lMatrix(inputmatrix.GetWidth(),inputmatrix.GetHeight());
    for(int i=0; i<lMatrix.GetHeight(); i++)
    {
        for(int j=0; j<lMatrix.GetWidth(); j++)
        {
            lMatrix.SetValue(j,i,ScaleOperation::GetL(myPyramid,j,i,sigmaL));
        }

    }
    std::cout<<"Saving image with L(x,y,sigma)..."<<std::endl;
    ImageAccessor::DrawImageFromMatrix(lMatrix,pyramidPath+"image_L(x,y,sigma)"+extension);
    std::cout<<"Done!"<<std::endl;
}
