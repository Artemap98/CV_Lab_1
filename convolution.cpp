#include "convolution.h"
#include <QDebug>


GrayScaleMatrix Convolution::Convolute(GrayScaleMatrix inputGSMatrix, QVector<QVector<double>> convCore)
{

    int xRadius = convCore[0].size() / 2;
    int yRadius = convCore.size() / 2;

    int width = inputGSMatrix.GetWidth();
    int height = inputGSMatrix.GetHeight();

    GrayScaleMatrix outputMatrix(width,height);
    for (int i = 0; i < height; i++)
        for (int j = 0; j < width; j++)
        {
            double sum = 0; //результат свертки для одной точки
            for (int u = -yRadius; u <= yRadius; u++)//для каждого ряда в ядре
                for (int v = -xRadius; v <= xRadius; v++)  //для каждого значения в ряду
                    sum += inputGSMatrix.GetValue(j - v, i - u) * convCore[u+yRadius][v+xRadius];
            outputMatrix.SetValue(j,i,sum);
        }
    return outputMatrix;
}


GrayScaleMatrix Convolution::DerivateX(GrayScaleMatrix inputGSMatrix) //получить массивы с частными производными
{
    QVector<QVector<double>> core;
    core.append(QVector<double>({1,0,-1}));
    core.append(QVector<double>({2,0,-2}));
    core.append(QVector<double>({1,0,-1}));

    return Convolute(inputGSMatrix, core);
}



GrayScaleMatrix Convolution::DerivateY(GrayScaleMatrix inputGSMatrix)
{
    QVector<QVector<double>> core;
    core.append(QVector<double>({1,2,1}));
    core.append(QVector<double>({0,0,0}));
    core.append(QVector<double>({-1,-2,-1}));
    return Convolute(inputGSMatrix, core);
}



GrayScaleMatrix Convolution::SobelOperator(GrayScaleMatrix derivateX,GrayScaleMatrix derivateY)
{

    int     width = derivateX.GetWidth(),
            height = derivateX.GetHeight();

    GrayScaleMatrix gradientMatrix(width,height);

    double derXij, derYij, gradij;

    for (int i = 0; i < height; i++)
    {//все строки

        for (int j = 0; j < width; j++)
        {
            derXij = derivateX.GetValue(j,i);
            derYij = derivateY.GetValue(j,i);
            gradij = sqrt(derXij*derXij + derYij*derYij);
            gradientMatrix.SetValue(j,i,gradij);
        }
    }
    return gradientMatrix;
}


GrayScaleMatrix Convolution::SobelOperator(GrayScaleMatrix inputGSMatrix)
{
    GrayScaleMatrix derivateX = DerivateX(inputGSMatrix);
    GrayScaleMatrix derivateY = DerivateY(inputGSMatrix);

    return SobelOperator(derivateX,derivateY);
}


GrayScaleMatrix Convolution::GaussianFilter(GrayScaleMatrix inputGSMatrix,double sigma)
{
//    QVector<QVector<double> > core; //ядро свертки

//    int sigmaInt = static_cast<int>(sigma) * 3;
//    if (sigmaInt  % 2 == 0)
//        sigmaInt++;
//    double coeff = 2 * sigma * sigma;

//    for (int i = -sigmaInt; i <= sigmaInt; i++)
//    {
//        QVector<double> str;
//        for (int j = -sigmaInt; j <= sigmaInt; j++)
//        {
//            str.append(exp( -(i * i + j * j) / coeff) / (M_PI * coeff));
//        }
//        core.append(str);
//    }
//    return Convolute(inputGSMatrix, core);

    //qDebug()<< "--Gauss. sigma = " << sigma;
    QVector<QVector<double> > core; //ядро свертки
    QVector<QVector<double> > core1;

    int sigmaInt = static_cast<int>(sigma) * 3;
    if (sigmaInt  % 2 == 0)
        sigmaInt++;
    double coeff = 2 * sigma * sigma;

    QVector<double> str;
    double coreSum = 0;
    for (int j = -sigmaInt; j <= sigmaInt; j++)
    {
        double gaussValue = exp( -(j * j) / coeff) / sqrt((M_PI * coeff));
        coreSum+= gaussValue;
        str.append(gaussValue);
        QVector<double> str1;
        str1.append(gaussValue);
        core1.append(str1);
    }
    core.append(str);

    for (int j = -sigmaInt; j <= sigmaInt; j++)
    {
        core[0][j+sigmaInt]/= coreSum;
        core1[j+sigmaInt][0]/= coreSum;
    }


    GrayScaleMatrix workMatr = Convolute(inputGSMatrix, core);
    return Convolute(workMatr, core1); //непосредственно вычисляем

}

GrayScaleMatrix Convolution::GradientDirection(GrayScaleMatrix derivateXMatrix, GrayScaleMatrix derivateYMatrix)
{
    int width = derivateXMatrix.GetWidth(), height = derivateXMatrix.GetHeight();
    GrayScaleMatrix outputGSMatrix(width, height);

    for(int i=0; i<height; i++)
    {
        for(int j=0; j<width; j++)
        {
            outputGSMatrix.SetValue(j,i,
                                    atan2(derivateYMatrix.GetValue(j,i),
                                          derivateXMatrix.GetValue(j,i))
                                    * 180 / M_PI + 180);
        }

    }
    return outputGSMatrix;
}

