#include "descriptorworker.h"
#include "imageaccessor.h"
#include <QPainter>
#define NNDR_NUM 0.5
#define VAL_MAX_COEFF 0.3


DescriptorWorker::DescriptorWorker()
{

}

GrayScaleMatrix DescriptorWorker::MergeTwoMatrix(GrayScaleMatrix inputGSMatrix1, GrayScaleMatrix inputGSMatrix2)
{
    //создаем матрицу с высотой бОльшего изображения и шириной равной сумме ширин этих иозображений
    GrayScaleMatrix resultGSMatrix(inputGSMatrix1.GetWidth()+inputGSMatrix2.GetWidth(),
                                   std::max(inputGSMatrix1.GetHeight(), inputGSMatrix2.GetHeight()));
    //из первого изображения пишем пиксели в матрицу как обычно
    for(int i=0; i<inputGSMatrix1.GetHeight(); i++)
    {
        for(int j=0; j< inputGSMatrix1.GetWidth(); j++)
        {
            resultGSMatrix.SetValue(j,i,inputGSMatrix1.GetValue(j,i));
        }
    }
    //из второго изображения пишем пиксели, отступая по Х на ширину первого изображения
    for(int i=0; i<inputGSMatrix2.GetHeight(); i++)
    {
        for(int j=0; j< inputGSMatrix2.GetWidth(); j++)
        {
            resultGSMatrix.SetValue(j+inputGSMatrix1.GetWidth(),i,inputGSMatrix2.GetValue(j,i));
        }
    }
    return resultGSMatrix;
}

double DescriptorWorker::GetDistanceBetweenDescriptors(Descriptor d1,Descriptor d2)
{
    double result = 0;
    int histogramsCount = d1.GetHistogramGridSize();
    int basketCount = d1.GetBasketNum();
    for(int i = 0; i < histogramsCount; i++){
        for(int j = 0; j < basketCount; j++){
            result += pow((d1.GetBasket(i,j) - d2.GetBasket(i,j)), 2);
        }
    }
    return sqrt(result);
}


double DescriptorWorker::GetDistanceBetweenPoints(int x1, int x2, int y1, int y2)
{
    return sqrt(pow((x1 - x2), 2) + pow((y1 - y2), 2));
}


QImage DescriptorWorker::GetTwoImageDescriptorComparsion(
            GrayScaleMatrix inputGSMatrix1, QVector<Descriptor> descriptors1,
            GrayScaleMatrix inputGSMatrix2, QVector<Descriptor> descriptors2)
{
    qDebug() << "--GetTwoImageDescriptorComparsion()";

    GrayScaleMatrix resultGSMatrix = MergeTwoMatrix(inputGSMatrix1,inputGSMatrix2);


    QImage resultImage(ImageAccessor::GetImageFromMatrix(resultGSMatrix));

    QPainter painter (&resultImage);
    painter.setPen(QColor(255, 255, 255, 200));
    int w = inputGSMatrix1.GetWidth();

    QVector<QVector<double>> distMatrix;
    double minValue = std::numeric_limits<double>::max();
    double maxValue = std::numeric_limits<double>::min();
    double middleValue = 0;

    //создаем матрицу расстояний между дескрипторов двух изображений
    for(int i = 0; i < descriptors1.size(); i++)
    {
        QVector<double> distRow;
        for(int j = 0; j < descriptors2.size(); j++)
        {
            double dist = GetDistanceBetweenDescriptors(descriptors1[i], descriptors2[j]);
            distRow.append(dist);

            middleValue += dist;
            if(dist < minValue){
                minValue = dist;
            }
            if(dist > maxValue){
                maxValue = dist;
            }
        }
        distMatrix.append(distRow);
    }
    middleValue /= descriptors1.size() * descriptors2.size();
    qDebug()<<"MinDist: " << minValue <<" MaxDist: " << maxValue << "middle: " << middleValue;



    qDebug() << "--Searching for pairs of near descriptors";

    //для каждого дескриптора первой картинки ищем соответствующий из второго, т.е. до которого наименьшее расстояние
    for(int i = 0; i < descriptors1.size(); i++){
        double firstMinValue = std::numeric_limits<double>::max();
        int firstMinIndex = 0;
        double secondMinValue = std::numeric_limits<double>::max();


        //для текущего дескриптора из первой картинки ищем ближайшие два дескриптора
        for(int j = 0; j < descriptors2.size(); j++){
            double dist = distMatrix[i][j];
            if(dist < firstMinValue){
                secondMinValue = firstMinValue;

                firstMinValue = dist;
                firstMinIndex = j;
            } else {
                if(dist < secondMinValue){
                    secondMinValue = dist;
                }
            }
        }

        //для найденного ближайшего дескриптора из второй матрицы ищем ближайшие из первой,
        //чтобы проверить, нет ли ложных срабатываний
        double firstMinValue2 = std::numeric_limits<double>::max();
        double secondMinValue2 = std::numeric_limits<double>::max();
        for(int j = 0; j < descriptors1.size(); j++){
            double dist = distMatrix[j][firstMinIndex];
            if(dist < firstMinValue2){
                secondMinValue2 = firstMinValue2;

                firstMinValue2 = dist;
            } else {
                if(dist < secondMinValue2){
                    secondMinValue2 = dist;
                }
            }
        }

        qDebug()<<"firstMinValue: " << firstMinValue <<" secondMinValue: " << secondMinValue;
        qDebug()<<"firstMinValue2: " << firstMinValue2 <<" secondMinValue2: " << secondMinValue2;

        //проверяем, есть ли другой достаточно близкий дескриптор,
        //а также является ли найденное расстояние достаточно близким
        if(firstMinValue / secondMinValue < NNDR_NUM
                && firstMinValue2 / secondMinValue2 < NNDR_NUM
                && firstMinValue < maxValue * VAL_MAX_COEFF)
        {
            //рисуем найденные точки на изображении
            for(int ii=-1; ii<=1; ii++)
            {
                for(int jj=-1; jj<=1; jj++)
                {
                    if(ii==0 || jj==0)
                    {
                        painter.drawPoint(descriptors1[i].GetX()+ii, descriptors1[i].GetY()+jj);
                        painter.drawPoint(descriptors2[firstMinIndex].GetX() + w +ii, descriptors2[firstMinIndex].GetY()+jj);

                    }
                }
            }

            //рисуем прямую между двумя найденными точками
            painter.drawLine(descriptors1[i].GetX(), descriptors1[i].GetY(), descriptors2[firstMinIndex].GetX() + w, descriptors2[firstMinIndex].GetY());
        }
    }

    return resultImage;

}


QVector<Descriptor> DescriptorWorker::GetDescriptorsFromImage(GrayScaleMatrix inputGSMatrix, 
                                                            int harrisRadius, 
                                                            int harrisPointsNum,
                                                            int basketNum, 
                                                            int histogramGridSize,
                                                            int descriptorSize)
{
    QVector<Descriptor> imageDescriptor;
    GrayScaleMatrix harrisMatrix = KeyFeatures::GetHarrisMatrix(inputGSMatrix, harrisRadius);
    KeyFeatures::KeyPointSet harrisPoints = KeyFeatures::GetPointsHarris(inputGSMatrix,harrisMatrix, harrisRadius, harrisPointsNum);

    GrayScaleMatrix dx = Convolution::DerivateX(inputGSMatrix);
    GrayScaleMatrix dy = Convolution::DerivateY(inputGSMatrix);
    GrayScaleMatrix gradientDirection = Convolution::GradientDirection(dx, dy);
    GrayScaleMatrix gradientMagnitude = Convolution::SobelOperator(dx, dy);
    
    double basketSize = 360. / basketNum;
    int descriptorRadius = descriptorSize / 2;
    int histogramRadius = histogramGridSize / 2;

    double sigma = static_cast<double>(histogramGridSize) / 6;
    QVector<QVector<double>> gaussKernel;

    double coeff = 1 / (2 * M_PI * sigma * sigma);
    double delitel = 2 * sigma * sigma;

    //формируем ядро Гаусса
    for (int u = -histogramRadius; u <= histogramRadius; u++)
    {
        QVector<double> gaussRow;
        for (int v = -histogramRadius; v <= histogramRadius; v++)
        {
            gaussRow.append( coeff * exp(- (u * u + v * v) / delitel));
        }
         gaussKernel.append(gaussRow);
    }
    

    qDebug() << "--computing histograms";
    int pointCount=0;

    foreach(KeyFeatures::KeyPoint keyPoint, harrisPoints.keyPoints)
    {
        pointCount++;
        int     x = keyPoint.x,
                y = keyPoint.y;
        Descriptor pointDescriptor(basketNum,histogramGridSize,descriptorSize,x, y);
        int currHist = 0;
        //для каждой гистограммы в дескрипторе
        for(int ih=-descriptorRadius; ih<descriptorRadius; ih++)
        {
            for(int jh=-descriptorRadius; jh<descriptorRadius; jh++)
            {
                //для каждого дескриптора в сетке гистограммы
                for(int ig=-histogramRadius; ig<histogramRadius; ig++)
                {
                    for(int jg=-histogramRadius; jg<histogramRadius; jg++)
                    {
                        int     currX=x+jh*histogramGridSize+jg, //абсолютное положение относительных координат в сетке
                                currY=y+ih*histogramGridSize+ig;

                        double currDirection = gradientDirection.GetValue(currX,currY);
                        currDirection = (currDirection < 0) ? currDirection + 360 : currDirection;
                        currDirection = (currDirection >= 360) ? currDirection - 360 : currDirection;
                        double currMagnitude = gradientMagnitude.GetValue(currX,currY);

                        //определяем, в какую корзиночку класть пирожок
                        double basketBetw = currDirection / basketSize;
                        int basket1 = floor(basketBetw);
                        double b1Weight = 1; //какую часть пирожка кладем в одну корзину

                        int basket2 = ceil(basketBetw);
                        double b2Weight = 0; //и какую часть в соседнюю


                        if(basketBetw < basket1 + 0.5)
                        {
                            basket2 = basket1 - 1;
                            if(basket2 < 0) basket2 = basketNum - 1;

                            b1Weight = abs(basketBetw - floor(basketBetw) + 0.5);
                        }
                        else
                        {
                            basket2 = basket1 + 1;
                            if(basket2 > basketNum - 1) basket2 = 0;

                            b1Weight = abs(basketBetw - floor(basketBetw) - 0.5);

                        }
                        b2Weight = 1. - b1Weight;

                        pointDescriptor.addValueToBasket(currHist, basket1, currMagnitude*b1Weight * gaussKernel[histogramRadius+ig][histogramRadius+jg]);
                        pointDescriptor.addValueToBasket(currHist, basket2, currMagnitude*b2Weight * gaussKernel[histogramRadius+ig][histogramRadius+jg]);
                    }
                }
                currHist++;
            }
        }
        pointDescriptor.NormalizeDescriptor();//приводим в божеский вид
        imageDescriptor.append(pointDescriptor);
    }
    return imageDescriptor;

}


QVector<Descriptor> DescriptorWorker::GetDescriptorsWithRotation(GrayScaleMatrix inputGSMatrix,
                                                                 int harrisRadius,
                                                                 int harrisPointsNum,
                                                                 int basketNum,
                                                                 int histogramGridSize,
                                                                 int descriptorSize)
{
    QVector<Descriptor> imageDescriptor;

    GrayScaleMatrix dx = Convolution::DerivateX(inputGSMatrix);
    GrayScaleMatrix dy = Convolution::DerivateY(inputGSMatrix);
    GrayScaleMatrix gradientDirection = Convolution::GradientDirection(dx, dy);
    GrayScaleMatrix gradientMagnitude = Convolution::SobelOperator(dx, dy);

    GrayScaleMatrix harrisMatrix = KeyFeatures::GetHarrisMatrix(inputGSMatrix, harrisRadius);
    KeyFeatures::KeyPointSet harrisPoints = KeyFeatures::GetPointsHarris(inputGSMatrix,harrisMatrix, harrisRadius, harrisPointsNum);
    //вычисляем уголы поворота интересных точек
    harrisPoints = DescriptorWorker::OrientPoints(harrisPoints, gradientDirection, gradientMagnitude, histogramGridSize, descriptorSize);


    double basketSize = 360. / basketNum;
    int descriptorRadius = descriptorSize / 2;
    int descriptorGridRadius = descriptorRadius * histogramGridSize;

    //находим ядро Гаусса
    double sigma = static_cast<double>(descriptorGridRadius) / 6;
    QVector<QVector<double>> gaussKernel;
    double coeff = 1 / (2 * M_PI * sigma * sigma);
    double delitel = 2 * sigma * sigma;
    for (int u = -descriptorGridRadius; u <= descriptorGridRadius; u++)
    {
        QVector<double> gaussRow;
        for (int v = -descriptorGridRadius; v <= descriptorGridRadius; v++)
        {
            gaussRow.append( coeff * exp(- (u * u + v * v) / delitel));
        }
         gaussKernel.append(gaussRow);
    }


    qDebug() << "--computing histograms";
    int pointCount=0;

    foreach(KeyFeatures::KeyPoint keyPoint, harrisPoints.keyPoints)
    {
        pointCount++;
        int     x = keyPoint.x,//координаты и угол текущей точки
                y = keyPoint.y,
                angle = keyPoint.angle;
        Descriptor pointDescriptor(basketNum,histogramGridSize,descriptorSize,x, y);
        for(int i=-descriptorGridRadius; i<=descriptorGridRadius; i++)
        {
            for(int j=-descriptorGridRadius; j<=descriptorGridRadius; j++)
            {
                        int     currX=x+i,//смещение координат в сетке
                                currY=y+j;

                        //повернутые координаты для определения, в какую гистограмму писать значения
                        int rotatedX = i * cos((360.0 - angle) * M_PI / 180.0) - j * sin((360.0 - angle) * M_PI / 180.0) + 0.5;
                        int rotatedY = i * sin((360.0 - angle) * M_PI / 180.0) + j * cos((360.0 - angle) * M_PI / 180.0) + 0.5;

                        double pointAngle = atan2(rotatedY, rotatedX ); //угол, образуемый отрезком
                        rotatedY -= sin(pointAngle) * qMax(rotatedX - 8, 0) + 0.5;
                        rotatedX -= cos(pointAngle) * qMax(rotatedY - 8, 0) + 0.5;

                            //За границей?
                            if(rotatedX < -descriptorGridRadius)
                            {
                                rotatedX = 0;
                            }
                            else
                            {
                                if(rotatedX >= descriptorGridRadius)
                                {
                                    rotatedX = descriptorGridRadius;
                                }
                                else
                                {
                                    rotatedX += descriptorGridRadius;
                                }
                            }

                            if(rotatedY < -descriptorGridRadius)
                            {
                                rotatedY = 0;
                            }
                            else
                            {
                                if(rotatedY >= descriptorGridRadius)
                                {
                                    rotatedY = descriptorGridRadius;
                                }
                                else
                                {
                                    rotatedY += descriptorGridRadius;
                                }
                            }
                        double currDirection = gradientDirection.GetValue(currX,currY) - angle;
                        currDirection = (currDirection < 0) ? currDirection + 360 : currDirection;
                        currDirection = (currDirection >= 360) ? currDirection - 360 : currDirection;
                        double currMagnitude = gradientMagnitude.GetValue(currX,currY);

                        //в какие корзины пишем значение
                        double basketBetw = currDirection / basketSize;

                        int basket1 = floor(basketBetw);
                        double b1Weight = 1;

                        int basket2 = ceil(basketBetw);
                        double b2Weight = 0;

                        if(basketBetw < basket1 + 0.5)
                        {
                            basket2 = basket1 - 1;
                            if(basket2 < 0) basket2 = basketNum - 1;

                            b1Weight = abs(basketBetw - floor(basketBetw) + 0.5);
                        }
                        else
                        {
                            basket2 = basket1 + 1;
                            if(basket2 > basketNum - 1) basket2 = 0;

                            b1Weight = abs(basketBetw - floor(basketBetw) - 0.5);

                        }
                        b2Weight = 1. - b1Weight;

                        //определяем гистограмму, в которую записываем значение
                        int histCol = static_cast<int>(static_cast<double>(rotatedX) / (histogramGridSize) - 0.2);
                        int histRow = static_cast<int>(static_cast<double>(rotatedY) / (histogramGridSize) - 0.2);
                        int currHist =  descriptorSize * histRow + histCol;

                        pointDescriptor.addValueToBasket(currHist, basket1, currMagnitude * b1Weight);
                        pointDescriptor.addValueToBasket(currHist, basket2, currMagnitude * b2Weight);
            }
        }
        pointDescriptor.NormalizeDescriptor();
        imageDescriptor.append(pointDescriptor);
    }
    qDebug() << "--histograms computed";
    return imageDescriptor;

}


KeyFeatures::KeyPointSet DescriptorWorker::OrientPoints(KeyFeatures::KeyPointSet inputPoints,
                                                        GrayScaleMatrix gradientDirection,
                                                        GrayScaleMatrix gradientMagnitude,
                                                        int histogramGridSize,
                                                        int descriptorSize)
{

    KeyFeatures::KeyPointSet orientedPoints;

    int localBasketCount = 36;  //число корзин
    double localBasketSize = 360.0 / localBasketCount;  //охват корзины

    int descriptorRadius = descriptorSize / 2;
    //int histogramRadius = histogramGridSize / 2;
    int descriptorGridRadius = descriptorRadius * histogramGridSize;

    //находим ядро Гаусса
    double sigma = static_cast<double>(descriptorGridRadius) / 6;
    QVector<QVector<double>> gaussKernel;
    double coeff = 1 / (2 * M_PI * sigma * sigma);
    double delitel = 2 * sigma * sigma;
    for (int u = -descriptorGridRadius; u <= descriptorGridRadius; u++)
    {
        QVector<double> gaussRow;
        for (int v = -descriptorGridRadius; v <= descriptorGridRadius; v++)
        {
            gaussRow.append( coeff * exp(-(u * u + v * v) / delitel));
        }
         gaussKernel.append(gaussRow);
    }

    for(int index = 0; index < inputPoints.keyPoints.size(); index++) {

        double localBaskets[localBasketCount];
        for (int i = 0; i < localBasketCount; i++){
            localBaskets[i] = 0;
        }

        KeyFeatures::KeyPoint currPoint = inputPoints.keyPoints[index];

        for(int i = -descriptorGridRadius; i <= descriptorGridRadius; i++){
            for(int j = -descriptorGridRadius; j <= descriptorGridRadius; j++){

                //В пределах ?
                if(GetDistanceBetweenPoints(i,j, 0.0,0.0) < sqrt(pow(descriptorGridRadius,2) + pow(descriptorGridRadius,2))){

                    //Направление Фи
                    double direction =  gradientDirection.GetValue(currPoint.x + i, currPoint.y + j);

                    //в какую корзину пишем
                    double basketBetw = direction / localBasketSize;

                    int basket1 = floor(basketBetw);
                    double b1Weight = 1;

                    int basket2 = ceil(basketBetw);
                    double b2Weight = 0;

                    if(basketBetw < basket1 + 0.5)
                    {
                        basket2 = basket1 - 1;
                        if(basket2 < 0) basket2 = localBasketCount - 1;

                        b1Weight = abs(basketBetw - floor(basketBetw) + 0.5);
                    }
                    else
                    {
                        basket2 = basket1 + 1;
                        if(basket2 > localBasketCount - 1) basket2 = 0;

                        b1Weight = abs(basketBetw - floor(basketBetw) - 0.5);

                    }
                    b2Weight = 1. - b1Weight;

                    double currMagnitude = gradientMagnitude.GetValue(currPoint.x + i,currPoint.y + j);

                    localBaskets[basket1] += currMagnitude
                            * b1Weight
                            * gaussKernel[(i + descriptorGridRadius)][(j + descriptorGridRadius)];

                    localBaskets[basket2] += currMagnitude
                            * b2Weight
                            * gaussKernel[(i + descriptorGridRadius)][(j + descriptorGridRadius)];
                }
            }
        }

        double firstMaxValue = -1;
        int firstMaxIndex = -1;
        double secondMaxValue = -1;
        int secondMaxIndex = -1;

        //ищем первую и вторую максимальную
        for(int i = 0; i < localBasketCount; i++){
            if(localBaskets[i] > firstMaxValue){
                secondMaxValue = firstMaxValue;
                secondMaxIndex = firstMaxIndex;

                firstMaxValue = localBaskets[i];
                firstMaxIndex = i;
            } else {
                if(localBaskets[i] > secondMaxValue){
                    secondMaxValue = localBaskets[i];
                    secondMaxIndex = i;
                }
            }
        }

        //добавляем первую
        KeyFeatures::KeyPoint firstPoint(inputPoints.keyPoints[index]);
        firstPoint.angle = (firstMaxIndex * localBasketSize);
        orientedPoints.keyPoints.push_back(firstPoint);

        //если вторая корзина >= 0.8 от макс значения, то добваляем ее тоже
        if(secondMaxValue >= (firstMaxValue * 0.8)){
            KeyFeatures::KeyPoint otherPoint(inputPoints.keyPoints[index]);
            otherPoint.angle = (secondMaxIndex * localBasketSize);
            orientedPoints.keyPoints.push_back(otherPoint);
        }

    }
    return orientedPoints;
}


