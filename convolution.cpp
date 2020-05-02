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
//    double *imageWorking = new double [widthWorking * heightWorking]; //расширенное рабочее изображение
    GrayScaleMatrix workMatrix(workingWidth,workingHeight);


    //Заполняем рабочее изображение
    double curOffsetX = 0, curOffsetY = 0;
    for (int i = 0; i < workingHeight; i++)
    {
        for (int j = 0; j < workingWidth; j++)
        {
            //точки за пределами исходного изображения приравниваем граничным
            curOffsetY = i-coreHeight;
            curOffsetX = j-coreWidth;
            if(i < coreHeight)
            {
                curOffsetY = 0;
            } else if(workingHeight - i - 1 < coreHeight * 2)
            {
                curOffsetY = workingHeight - coreHeight * 2 -1;
            }

            if(j < coreWidth)
            {
                curOffsetX = 0;
            } else if(workingWidth - j - 1 < coreWidth * 2)
            {
                curOffsetX = workingWidth - coreWidth * 2 -1;
            }

            workMatrix.SetValue(j,i,inputGSMatrix.GetValue(curOffsetX,curOffsetY));
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

            for (int u = 0; u < ySize; u++)//для каждого ряда в ядре
            {
                for (int v = 0; v < xSize; v++)  //для каждого значения в ряду
                {
                    sum += workMatrix.GetValue(j-v+xSize-1, i - u +ySize-1) * convCore[u][v];
                }
            }
            outputMatrix.SetValue(j,i,sum);
        }
    }


//    double min=0, max=1;
//    for(int i=0; i<outputMatrix.GetHeight(); i++)
//    {
//        for(int j=0; j<outputMatrix.GetWidth(); j++)
//        {
//            if(outputMatrix.GetValue(j,i)>max)
//                max=outputMatrix.GetValue(j,i);
//            if(outputMatrix.GetValue(j,i)<min)
//                outputMatrix.GetValue(j,i);
//        }
//    }
//    for(int i=0; i<outputMatrix.GetHeight(); i++)
//    {
//        for(int j=0; j<outputMatrix.GetWidth(); j++)
//        {
//            outputMatrix.SetValue(j,i,(outputMatrix.GetValue(j,i)-min)/(max-min));
//        }
//    }

    return outputMatrix;
}


QVector<QVector<double>> Convolution::MatrixMult(QVector<QVector<double>> matrix1, QVector<QVector<double>> matrix2)
{

    QVector<QVector<double>> resultMatrix;
    for(int i = 0; i<matrix1.size(); i++)
    {
        QVector<double> matrixRow;
        for(int j=0; j<matrix2[0].size(); j++)
        {
            double elementSum=0;
            for(int k =0; k<matrix1.size(); k++)
            {
                elementSum += matrix1[i][k]*matrix2[k][j];
            }
        }
    }
    return resultMatrix;
}



GrayScaleMatrix Convolution::GetDerivateX(GrayScaleMatrix inputGSMatrix) //получить массивы с частными производными
{
    QVector<QVector<double>> core;
    core.append(QVector<double>({1,0,-1}));
    core.append(QVector<double>({2,0,-2}));
    core.append(QVector<double>({1,0,-1}));

    return Convolute(inputGSMatrix, core);
}



GrayScaleMatrix Convolution::GetDerivateY(GrayScaleMatrix inputGSMatrix)
{
    QVector<QVector<double>> core;
    core.append(QVector<double>({1,2,1}));
    core.append(QVector<double>({0,0,0}));
    core.append(QVector<double>({-1,-2,-1}));
    return Convolute(inputGSMatrix, core);
}



GrayScaleMatrix Convolution::SobelOperator(GrayScaleMatrix inputGSMatrix)
{
    GrayScaleMatrix derivateX = GetDerivateX(inputGSMatrix);
    GrayScaleMatrix derivateY = GetDerivateY(inputGSMatrix);

    int     width = inputGSMatrix.GetWidth(),
            height = inputGSMatrix.GetHeight();

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



GrayScaleMatrix Convolution::GaussianFilter(GrayScaleMatrix inputGSMatrix,double sigma)
{
    QVector<QVector<double> > core; //ядро свертки

    int sigmaHalf = static_cast<int>(sigma) * 3;
    if (sigmaHalf  % 2 == 0)
        ++sigmaHalf;
    double coeff = 2 * sigma * sigma;

    for (int i = -sigmaHalf; i <= sigmaHalf; i++)
    {
        QVector<double> str;
        for (int j = -sigmaHalf; j <= sigmaHalf; j++)
        {
            str.append(exp( -(i * i + j * j) / coeff) / (M_PI * coeff));
        }
        core.append(str);
    }

    return Convolute(inputGSMatrix, core); //непосредственно вычисляем
}


