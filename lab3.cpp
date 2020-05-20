//◦ Реализовать операторы Моравека и Харриса для поиска интересных точек в изображении

//◦ Реализовать фильтрацию интересных точек методом Adaptive Non-Maximum Suppression
    //для заданного количества необходимых точек

//◦ Оценить повторяемость результата при некоторых искажениях оригинального изображения –
    //сдвиг, поворот, шум, контрастность и яркость

//◦ Сравнить выдачу операторов Моравека и Харриса по повторяемости

//◦ (Пов. сложность) Реализовать алгоритм поиска краев Кэнни


#include "imageaccessor.h"
#include "KeyFeatures.h"

#include <QDir>
#include <QString>
#include <iostream>

void lab3(QString path, QString fileName, QString extension,int windowRadius, int resultPointsNum)
{
    QDir dir;
    QString labPath = path+"lab3"+fileName+"\\";
    dir.mkdir(labPath);

    std::cout<<"load img..."<<std::endl;
    GrayScaleMatrix inputMatrix = ImageAccessor::GetMatrixFromImage(path+fileName+extension);

    std::cout<<"compute harris response..."<<std::endl;
    GrayScaleMatrix harrisMatrix = KeyFeatures::GetHarrisMatrix(inputMatrix, windowRadius);
    std::cout<<"draw harris matrix..."<<std::endl;
    ImageAccessor::DrawImageFromMatrix(harrisMatrix,labPath+"harrisResponse"+extension);

    std::cout<<"compute harris points (may take a couple of minutes)..."<<std::endl;
    KeyFeatures::KeyPointSet harrisPoints = KeyFeatures::GetPointsHarris(inputMatrix,harrisMatrix, windowRadius, resultPointsNum);

    GrayScaleMatrix harrisPointsMatrix = KeyFeatures::GetResultHarrisMatrix(inputMatrix, harrisPoints);

    std::cout<<"draw harris points on img..."<<std::endl;
    ImageAccessor::DrawImageFromMatrix(harrisPointsMatrix,labPath+"harrisPoints"+extension);

    std::cout<<"compute moravec response..."<<std::endl;
    GrayScaleMatrix moravecMatrix = KeyFeatures::GetMoravecMatrix(inputMatrix, windowRadius);
    std::cout<<"draw moravec matrix..."<<std::endl;
    ImageAccessor::DrawImageFromMatrix(moravecMatrix,labPath+"moravecResponse"+extension);

    std::cout<<"compute moravec points (may take a couple of minutes)..."<<std::endl;
    GrayScaleMatrix moravecPoints = KeyFeatures::GetResultMoravec(inputMatrix,moravecMatrix, windowRadius, resultPointsNum);
    std::cout<<"draw moravec points on img..."<<std::endl;
    ImageAccessor::DrawImageFromMatrix(moravecPoints,labPath+"moravecPoints"+extension);




    std::cout<<"done!"<<std::endl;
}
