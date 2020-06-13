#include "convolution.h"


GrayScaleMatrix Convolution::Convolute(GrayScaleMatrix inputGSMatrix, QVector<QVector<double>> convCore)
{
//    int coreW = static_cast<int>(core[0].count() / 2);
//    int coreH = static_cast<int>(core.count() / 2);

    int coreWidth = convCore[0].size() / 2;
    int coreHeight = convCore.size() / 2;

    int workingWidth = inputGSMatrix.GetWidth() + coreWidth * 2;    //ширина и высота расширенного, рабочего изображения
    int workingHeight = inputGSMatrix.GetHeight() + coreHeight * 2;

    //расширяем изображение на размер ядра
    //расширенное рабочее изображение
    GrayScaleMatrix workMatrix(workingWidth,workingHeight);


    //Заполняем рабочее изображение
    for (int i = -coreHeight; i < workingHeight - coreHeight; i++)
    {
        for (int j = 0-coreWidth; j < workingWidth-coreWidth; j++)
        {
            //точки за пределами исходного изображения приравниваем граничным
            workMatrix.SetValue(j+coreWidth,i+coreHeight,inputGSMatrix.GetValue(j,i));
        }
    }

    int xSize = convCore[0].size();
    int ySize = convCore.size();

    //вычисляем сумму элементов ядра
    double coreSum = 0;
    for (int x = 0; x < xSize; ++x) {
        for (int y = 0; y < ySize; ++y) {
            coreSum += convCore[y][x];
        }
    }

    GrayScaleMatrix outputMatrix(inputGSMatrix.GetWidth(),inputGSMatrix.GetHeight());
    int     width = outputMatrix.GetWidth(),
            height = outputMatrix.GetHeight();
    //применяем свертку
    for (int i = 0; i < height; i++)
    {//все строки
        for (int j = 0; j < width; j++)
        {
            double sum = 0; //результат свертки для одной точки

            for (int u = -coreHeight; u <= coreHeight; u++)//для каждого ряда в ядре
            {
                for (int v = -coreWidth; v <= coreWidth; v++)  //для каждого значения в ряду
                {
                    sum += workMatrix.GetValue(j - v + coreWidth, i - u + coreHeight) * convCore[u+coreHeight][v+coreWidth];
                }
            }
            outputMatrix.SetValue(j,i,sum);
        }
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

    QVector<QVector<double> > core; //ядро свертки
    QVector<QVector<double> > core1;

    int sigmaInt = static_cast<int>(sigma) * 3;
    if (sigmaInt  % 2 == 0)
        sigmaInt++;
    double coeff = 2 * sigma * sigma;

    QVector<double> str;
    for (int j = -sigmaInt; j <= sigmaInt; j++)
    {
        str.append(exp( -(j * j) / coeff) / sqrt((M_PI * coeff)));
    }
    core.append(str);

    GrayScaleMatrix workMatr = Convolute(inputGSMatrix, core);

    for (int j = -sigmaInt; j <= sigmaInt; j++)
    {
        QVector<double> str1;
        str1.append(exp( -(j * j) / coeff) / sqrt((M_PI * coeff)));
        core1.append(str1);
    }

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

