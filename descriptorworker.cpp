#include "descriptorworker.h"
#include "imageaccessor.h"
#include "scaleoperation.h"
#include <QPainter>
#include <QtAlgorithms>
#define NNDR_NUM 0.8
#define VAL_MAX_COEFF 0.3
#define EPSILON 0.15


DescriptorWorker::DescriptorWorker()
{

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

                        pointDescriptor.addValueToBasket(currHist, basket1, currMagnitude*b1Weight);
                        pointDescriptor.addValueToBasket(currHist, basket2, currMagnitude*b2Weight);
                    }
                }
                currHist++;
            }
        }
        pointDescriptor.NormalizeDescriptor();
        imageDescriptor.append(pointDescriptor);
    }
    return imageDescriptor;

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
    painter.setPen(QColor(40, 255, 240, 170));
    int w = inputGSMatrix1.GetWidth();

    QVector<QVector<double>> distMatrix;
    QVector<double> distVector;
    double minValue = std::numeric_limits<double>::max();
    double maxValue = std::numeric_limits<double>::min();
    double middleValue = 0;
    double modalValue = 0;

    //создаем матрицу расстояний между дескрипторов двух изображений
    for(int i = 0; i < descriptors1.size(); i++)
    {
        QVector<double> distRow;
        for(int j = 0; j < descriptors2.size(); j++)
        {
            double dist = GetDistanceBetweenDescriptors(descriptors1[i], descriptors2[j]);
            distRow.append(dist);
            distVector.append(dist);
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
    //qSort(distVector);
    std::sort(distVector.begin(),distVector.end());
    //modalValue = distVector[qMin(descriptors1.size(),descriptors2.size())];
    modalValue = distVector[2.5*qMin(descriptors1.size(),descriptors2.size())];
    middleValue /= descriptors1.size() * descriptors2.size();
    qDebug()<<"MinDist: " << minValue <<" MaxDist: " << maxValue << " middle: " << middleValue<< " modal: " << modalValue;



    qDebug() << "--Searching for pairs of near descriptors";
    QList<int> foundList;
    //для каждого дескриптора первой картинки ищем соответствующий из второго, т.е. до которого наименьшее расстояние
    for(int i = 0; i < descriptors1.size(); i++)
    {
        double firstMinValue = std::numeric_limits<double>::max();
        int firstMinIndex = 0;
        double secondMinValue = std::numeric_limits<double>::max();


        //для текущего дескриптора из первой картинки ищем ближайшие два дескриптора
        for(int j = 0; j < descriptors2.size(); j++)
        {
            double dist = distMatrix[i][j];
            if(dist < firstMinValue)
            {
                secondMinValue = firstMinValue;

                firstMinValue = dist;
                firstMinIndex = j;
            } else
            {
                if(dist < secondMinValue)
                {
                    secondMinValue = dist;
                }
            }
        }

        //для найденного ближайшего дескриптора из второй матрицы ищем ближайшие из первой,
        //чтобы проверить, нет ли ложных срабатываний
        double firstMinValue2 = std::numeric_limits<double>::max();
        double secondMinValue2 = std::numeric_limits<double>::max();
        for(int j = 0; j < descriptors1.size(); j++)
        {
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


        //проверяем, есть ли другой достаточно близкий дескриптор,
        //а также является ли найденное расстояние достаточно близким
        if(foundList.indexOf(firstMinIndex) == -1
                && firstMinValue <= modalValue)
        if(firstMinValue / secondMinValue < NNDR_NUM)
        if(firstMinValue2 / secondMinValue2 < NNDR_NUM)
        {
            foundList.append(firstMinIndex);
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
            qDebug()<<"Similar descriptors:[" << descriptors1[i].GetX()<< ";" << descriptors1[i].GetY() <<"] and [" << descriptors2[firstMinIndex].GetX()<< ";" << descriptors2[firstMinIndex].GetY() <<"]"<<" Distance = " << firstMinValue;
            qDebug()<<"firstMinValue: " << firstMinValue <<" secondMinValue: " << secondMinValue;
            qDebug()<<"firstMinValue2: " << firstMinValue2 <<" secondMinValue2: " << secondMinValue2;
            //рисуем прямую между двумя найденными точками
            painter.drawLine(descriptors1[i].GetX(), descriptors1[i].GetY(), descriptors2[firstMinIndex].GetX() + w, descriptors2[firstMinIndex].GetY());
        }
    }

    return resultImage;

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
    double gaussSum = 0;
    for (int u = -descriptorGridRadius; u <= descriptorGridRadius; u++)
    {
        QVector<double> gaussRow;
        for (int v = -descriptorGridRadius; v <= descriptorGridRadius; v++)
        {
            double value = coeff * exp(- (u * u + v * v) / delitel);
            gaussSum += value;
            gaussRow.append(value);
        }
         gaussKernel.append(gaussRow);
    }
//    for (int u = -descriptorGridRadius; u <= descriptorGridRadius; u++)
//        for (int v = -descriptorGridRadius; v <= descriptorGridRadius; v++)
//            gaussKernel[u+descriptorGridRadius][v+descriptorGridRadius] /= gaussSum;



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
                if(sqrt(i*i+j*j) < sqrt(2)*descriptorGridRadius){
                    //повернутые координаты для определения, в какую гистограмму писать значения
                    double angleCos = cos(angle * M_PI / 180.0);
                    double angleSin = sin(angle * M_PI / 180.0);

//                    double tempX = i * angleCos + j * angleSin;
//                    double tempY = j * angleCos - i * angleSin;
//                    int rotatedX = tempX;
//                    int rotatedY = tempY;
                    int rotatedX = i * angleCos + j * angleSin;
                    int rotatedY = j * angleCos - i * angleSin;

                    if(rotatedX > descriptorGridRadius)
                        rotatedX = descriptorGridRadius;
                    else if(rotatedX < -descriptorGridRadius)
                        rotatedX = -descriptorGridRadius;

                    if(rotatedY > descriptorGridRadius)
                        rotatedY = descriptorGridRadius;
                    else if(rotatedY < -descriptorGridRadius)
                        rotatedY = - descriptorGridRadius;

                    int     currX=x+i,//смещение координат в сетке
                            currY=y+j;

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

//                    b1Weight *= gaussKernel[(i+descriptorGridRadius)][(j+descriptorGridRadius)];
//                    b2Weight *= gaussKernel[(i+descriptorGridRadius)][(j+descriptorGridRadius)];

                    //определяем гистограмму, в которую записываем значение
                    int histCol = ((rotatedX+descriptorGridRadius) / (histogramGridSize) - 0.1);
                    int histRow = ((rotatedY+descriptorGridRadius) / (histogramGridSize) - 0.1);
                    int currHist =  descriptorSize * histRow + histCol;

                    pointDescriptor.addValueToBasket(currHist, basket1, currMagnitude * b1Weight);
                    pointDescriptor.addValueToBasket(currHist, basket2, currMagnitude * b2Weight);
                }
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
    double gaussSum = 0;
    for (int u = -descriptorGridRadius; u <= descriptorGridRadius; u++)
    {
        QVector<double> gaussRow;
        for (int v = -descriptorGridRadius; v <= descriptorGridRadius; v++)
        {
            double gaussVal = coeff * exp(-(u * u + v * v) / delitel);
            gaussRow.append(gaussVal);
            gaussSum = gaussVal;
        }
         gaussKernel.append(gaussRow);
    }

    for (int u = -descriptorGridRadius; u <= descriptorGridRadius; u++)
        for (int v = -descriptorGridRadius; v <= descriptorGridRadius; v++)
            gaussKernel[u+descriptorGridRadius][v+descriptorGridRadius] /= gaussSum;

    for(int index = 0; index < inputPoints.keyPoints.size(); index++) {

        double localBaskets[localBasketCount];
        for (int i = 0; i < localBasketCount; i++){
            localBaskets[i] = 0;
        }

        KeyFeatures::KeyPoint currPoint = inputPoints.keyPoints[index];

        for(int i = -descriptorGridRadius; i <= descriptorGridRadius; i++){
            for(int j = -descriptorGridRadius; j <= descriptorGridRadius; j++){

                //В пределах ?
                if(sqrt(i*i+j*j) < sqrt(2)*descriptorGridRadius){

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



QVector<Descriptor> DescriptorWorker::GetDescriptorsBlob(GrayScaleMatrix inputGSMatrix,
                                                                 int harrisRadius,
                                                                 int harrisPointsNum,
                                                                 int basketNum,
                                                                 int histogramGridSize,
                                                                 int descriptorSize)
{
    QVector<Descriptor> imageDescriptor;
    qDebug()<<"getting Gauss pyramid";
    ScaleOperation::Pyramid DoGPyramid;
    ScaleOperation::Pyramid inputPyramid = ScaleOperation::GetPyramidWithOverlay(inputGSMatrix,3,5,0,1.6);
    qDebug()<<"getting DoG pyramid";
    //для каждого уровня в каждой октаве находим значения Difference of Gaussian
    for(int i = 0; i < inputPyramid.octaves.size(); i++)
    {
        qDebug()<<"-octave "<< i;
        ScaleOperation::Octave currOctave;
        currOctave.level = i;
        //идем до предпоследнего элемента, т.к. из текущего слоя вычитаем значение следующего
        for(int j = 0; j < inputPyramid.octaves[i].layers.size()-1; j++)
        {
            qDebug()<<"--layer "<< j;
            ScaleOperation::Layer currLayer = inputPyramid.octaves[i].layers[j];
            ScaleOperation::Layer nextLayer = inputPyramid.octaves[i].layers[j+1];
//            currLayer.actualSigma = sqrt(currLayer.actualSigma * currLayer.actualSigma
//                                         + nextLayer.actualSigma*nextLayer.actualSigma);
//            currLayer.currentSigma = sqrt(currLayer.currentSigma * currLayer.currentSigma
//                                         + nextLayer.currentSigma * nextLayer.currentSigma);

            for(int ii = 0; ii < currLayer.matrix.GetWidth(); ii++)
            {
                for(int jj = 0; jj < currLayer.matrix.GetHeight(); jj++)
                {
                    currLayer.matrix.SetValue(ii,jj, abs(currLayer.matrix.GetValue(ii,jj)-nextLayer.matrix.GetValue(ii,jj)));
                }
            }
            currOctave.layers.append(currLayer);
        }
        DoGPyramid.octaves.append(currOctave);
    }

    qDebug()<<"find max and min response 3d ";
    //теперь ищем макс отклик в 3д
    QVector<QVector<KeyFeatures::KeyPointSet>> pyramidPoints;
    for(int i=0; i < DoGPyramid.octaves.size(); i++)
    {
        QVector<KeyFeatures::KeyPointSet> octavePoints;
        qDebug()<<"-octave "<< i;
        //отступаем от начала и конца, чтобы иметь слой сверху и снизу.
        for(int j=1; j < DoGPyramid.octaves[i].layers.size()-1; j++)
        {
            qDebug()<<"--layer "<< j;
            ScaleOperation::Layer currDoGLayer = DoGPyramid.octaves[i].layers[j];
            ScaleOperation::Layer currLayer = inputPyramid.octaves[i].layers[j-1];
            KeyFeatures::KeyPointSet currDoGLayerPoints;
            GrayScaleMatrix currHarrisMatrix = KeyFeatures::GetHarrisMatrix(currLayer.matrix, harrisRadius);
            currHarrisMatrix.NormalizeDouble();
            for(int ii = 0; ii < currDoGLayer.matrix.GetWidth(); ii++)
            {
                for(int jj = 0; jj < currDoGLayer.matrix.GetHeight(); jj++)
                {
                    double responseVal = currDoGLayer.matrix.GetValue(ii,jj);// значение проверяемого пикселя
                    if(responseVal > 0.01 )//отсекаем силшком маленькие отклики
                    {
                        bool isMinimal = true, isMaximal = true, canContinue = true;
                        //для трех слоев октавы
                        for(int k = -1; k <= 1 && canContinue; k++)
                        {
                            //для каждого слоя окрестность вокруг точки
                            for(int zi = -1; zi <= 1 && canContinue; zi++)
                            {
                                for(int zj = -1; zj <= 1 && canContinue; zj++)
                                {
                                    //если текущая точка окрестности не является проверяемой
                                    if((k != 0 || zi != 0 || zj != 0))
                                    {
                                        //текущее значение окрестности
                                        double tempVal = DoGPyramid.octaves[i].layers[j+k].matrix.GetValue(ii+zi,jj+zj);
                                        //если значение в окрестности больше проверяемого,
                                        //то преверяемое не может быть максимальным в области
                                        if(tempVal > responseVal)
                                            isMaximal = false;
                                        //если значение в окрестности меньше проверяемого,
                                        //то преверяемое не может быть минимальным в области
                                        if(tempVal < responseVal)
                                            isMinimal = false;

                                        //если проверяемое не максимальное и не минимальное,
                                        //то дальше проверять нет смысла
                                        if((!isMinimal && !isMaximal))
                                        {
                                            canContinue = false;
                                        }
                                    }
                                }
                            }
                        }
                        //если точка минимальная или максимальная в окрестности
                        if(isMaximal || isMinimal )
                        {
                            if(currHarrisMatrix.GetValue(ii,jj) > 0.01)//отсекаем слабый отклик
                            {
                                KeyFeatures::KeyPoint currPoint(ii,jj,currHarrisMatrix.GetValue(ii,jj));
                                currDoGLayerPoints.keyPoints.append(currPoint);
                            }
                        }
                    }
                }
            }
            qDebug()<<"---KeyPoints found "<< currDoGLayerPoints.keyPoints.size();
            octavePoints.append(currDoGLayerPoints);
        }
        pyramidPoints.append(octavePoints);
    }

    //считаем суммарное найденое число точек в пирамиде
    int totalPointsNum=0;
    QVector<int> octavePointsNum;
    for(int i=0; i < pyramidPoints.size(); i++)
    {
        int currOctavePointsNum = 0;
        for(int j=0; j < pyramidPoints[i].size(); j++)
            currOctavePointsNum += pyramidPoints[i][j].keyPoints.size();

        totalPointsNum += currOctavePointsNum;
        octavePointsNum.append(currOctavePointsNum);
    }
    qDebug()<<"totalPointsNum= "<< totalPointsNum;
    qDebug()<<"Compute descriptors for pyramid ";
    for(int i=0; i < pyramidPoints.size(); i++)
    {
        qDebug()<<"Octave "<<i;
        for(int j=0; j < pyramidPoints[i].size(); j++)
        {
            int currLevelReducedNum = harrisPointsNum * (double)pyramidPoints[i][j].keyPoints.size() / (double)totalPointsNum;
            qDebug()<<"---Level "<< j<<"currLevelPointsNum= "<< pyramidPoints[i][j].keyPoints.size()<<" currLevelReducedNum= "<< currLevelReducedNum;
            int maxRadius = qMin(inputGSMatrix.GetWidth(), inputGSMatrix.GetWidth()) / 2 / pow(2, i);
            GrayScaleMatrix dx = Convolution::DerivateX(inputPyramid.octaves[i].layers[j+1].matrix);
            GrayScaleMatrix dy = Convolution::DerivateY(inputPyramid.octaves[i].layers[j+1].matrix);
            GrayScaleMatrix gradientDirection = Convolution::GradientDirection(dx, dy);
            GrayScaleMatrix gradientMagnitude = Convolution::SobelOperator(dx, dy);
            pyramidPoints[i][j] = KeyFeatures::ReducePoints(pyramidPoints[i][j],
                                                                      currLevelReducedNum,maxRadius);

            int localHistogramGridSize = histogramGridSize * DoGPyramid.octaves[i].layers[j].currentSigma/DoGPyramid.octaves[i].layers[0].currentSigma;
            QVector<Descriptor> layerDescriptors = GetDescriptorsFromOneLayer(gradientDirection,
                                                                              gradientMagnitude,
                                                                              pyramidPoints[i][j],
                                                                              basketNum,
                                                                              localHistogramGridSize,
                                                                              descriptorSize);

            //преобразуем координаты изображения октавы в оригинальные координаты
            int coordMult = pow(2., i);
            foreach(Descriptor descr, layerDescriptors)
            {
                int xnew = static_cast<int>(descr.GetX() * coordMult);
                int ynew = static_cast<int>(descr.GetY() * coordMult);
                if(xnew >= inputGSMatrix.GetWidth()) xnew = inputGSMatrix.GetWidth()-1;
                if(ynew >= inputGSMatrix.GetHeight()) ynew =inputGSMatrix.GetHeight()-1;
                descr.SetX(xnew);
                descr.SetY(ynew);
                imageDescriptor.append(descr);
            }
            qDebug()<<"--Layer descriptor found ";
        }
    }
    qDebug()<<"Total descriptors found " << imageDescriptor.size();
    return imageDescriptor;
}



QVector<Descriptor> DescriptorWorker::GetDescriptorsFromOneLayer(GrayScaleMatrix gradientDirection,
                                                                 GrayScaleMatrix gradientMagnitude,
                                                                 KeyFeatures::KeyPointSet harrisPoints,
                                                                 int basketNum,
                                                                 int histogramGridSize,
                                                                 int descriptorSize)
{
    QVector<Descriptor> imageDescriptor;

//    GrayScaleMatrix dx = Convolution::DerivateX(inputGSMatrix);
//    GrayScaleMatrix dy = Convolution::DerivateY(inputGSMatrix);
//    GrayScaleMatrix gradientDirection = Convolution::GradientDirection(dx, dy);
//    GrayScaleMatrix gradientMagnitude = Convolution::SobelOperator(dx, dy);

    //GrayScaleMatrix harrisMatrix = KeyFeatures::GetHarrisMatrix(inputGSMatrix, harrisRadius);
    //KeyFeatures::KeyPointSet harrisPoints = KeyFeatures::GetPointsHarris(inputGSMatrix,harrisMatrix, harrisRadius, harrisPointsNum);
    //вычисляем уголы поворота интересных точек
    harrisPoints = DescriptorWorker::OrientPoints(harrisPoints, gradientDirection, gradientMagnitude, histogramGridSize, descriptorSize);


    double basketSize = 360. / basketNum;
    int descriptorRadius = descriptorSize / 2;
    int descriptorGridRadius = descriptorRadius * histogramGridSize;

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
                if(sqrt(i*i+j*j) < sqrt(2)*descriptorGridRadius){
                    //повернутые координаты для определения, в какую гистограмму писать значения
                    double angleCos = cos(angle * M_PI / 180.0);
                    double angleSin = sin(angle * M_PI / 180.0);

//                    double tempX = i * angleCos + j * angleSin;
//                    double tempY = j * angleCos - i * angleSin;
//                    int rotatedX = tempX;
//                    int rotatedY = tempY;
                    double rotatedX = i * angleCos + j * angleSin;
                    double rotatedY = j * angleCos - i * angleSin;

                    if(rotatedX > descriptorGridRadius)
                        rotatedX = descriptorGridRadius;
                    else if(rotatedX < -descriptorGridRadius)
                        rotatedX = -descriptorGridRadius;

                    if(rotatedY > descriptorGridRadius)
                        rotatedY = descriptorGridRadius;
                    else if(rotatedY < -descriptorGridRadius)
                        rotatedY = - descriptorGridRadius;

                    double     currX=x+i,//смещение координат в сетке
                            currY=y+j;

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

//                    b1Weight *= gaussKernel[(i+descriptorGridRadius)][(j+descriptorGridRadius)];
//                    b2Weight *= gaussKernel[(i+descriptorGridRadius)][(j+descriptorGridRadius)];

                    //определяем гистограмму, в которую записываем значение
                    double histCol = (rotatedX+descriptorGridRadius) / (histogramGridSize) - 0.1;
                    double histRow = (rotatedY+descriptorGridRadius) / (histogramGridSize) - 0.1;
                    int histCol1 = floor(histCol+0.5);
                    histCol1 = histCol1 >= descriptorSize ? descriptorSize-1 : (histCol1 < 0 ? 0 : histCol1);

                    int histCol2 = floor(histCol-0.5);
                    histCol2 = histCol2 >= descriptorSize ? descriptorSize-1 : (histCol2 < 0 ? 0 : histCol2);

                    int histRow1 = floor(histRow+0.5);
                    histRow1 = histRow1 >= descriptorSize ? descriptorSize-1 : (histRow1 < 0 ? 0 : histRow1);

                    int histRow2 = floor(histRow-0.5);
                    histRow2 = histRow2 >= descriptorSize ? descriptorSize-1 : (histRow2 < 0 ? 0 : histCol1);

                    int hist1 = descriptorSize * histRow1 + histCol1;
                    int hist2 = descriptorSize * histRow1 + histCol2;
                    int hist3 = descriptorSize * histRow2 + histCol1;
                    int hist4 = descriptorSize * histRow2 + histCol2;

                    double hist1w = sqrt(pow(histCol-histCol1,2) + pow(histRow-histRow1,2));
                    double hist2w = sqrt(pow(histCol-histCol1,2) + pow(histRow-histRow2,2));
                    double hist3w = sqrt(pow(histCol-histCol2,2) + pow(histRow-histRow1,2));
                    double hist4w = sqrt(pow(histCol-histCol2,2) + pow(histRow-histRow2,2));

                    double sumw = hist1w + hist2w + hist3w + hist4w;

//                    pointDescriptor.addValueToBasket(currHist, basket1, currMagnitude * b1Weight);
//                    pointDescriptor.addValueToBasket(currHist, basket2, currMagnitude * b2Weight);
                    pointDescriptor.addValueToBasket(hist1, basket1, currMagnitude * b1Weight * hist1w / sumw);
                    pointDescriptor.addValueToBasket(hist1, basket2, currMagnitude * b2Weight * hist1w / sumw);

                    pointDescriptor.addValueToBasket(hist2, basket1, currMagnitude * b1Weight * hist2w / sumw);
                    pointDescriptor.addValueToBasket(hist2, basket2, currMagnitude * b2Weight * hist2w / sumw);

                    pointDescriptor.addValueToBasket(hist3, basket1, currMagnitude * b1Weight * hist3w / sumw);
                    pointDescriptor.addValueToBasket(hist3, basket2, currMagnitude * b2Weight * hist3w / sumw);

                    pointDescriptor.addValueToBasket(hist4, basket1, currMagnitude * b1Weight * hist4w / sumw);
                    pointDescriptor.addValueToBasket(hist4, basket2, currMagnitude * b2Weight * hist4w / sumw);
                }
            }
        }
        pointDescriptor.NormalizeDescriptor();
        imageDescriptor.append(pointDescriptor);
    }
    qDebug() << "--histograms computed";
    return imageDescriptor;

}

