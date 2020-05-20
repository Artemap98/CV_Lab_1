#include "descriptorworker.h"
#include "imageaccessor.h"
#include <QPainter>

DescriptorWorker::DescriptorWorker()
{

}


double DescriptorWorker::GetDistance(Descriptor d1,Descriptor d2)
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
    
//    double  w = inputGSMatrix.GetWidth(),
//            h = inputGSMatrix.GetHeight();
    double basketSize = 360. / basketNum;
    int descriptorRadius = descriptorSize / 2;
    int histogramRadius = histogramGridSize / 2;


    //находим ядро Гаусса
    double sigma = static_cast<double>(histogramGridSize) / 6;
    QVector<QVector<double>> gaussKernel;

    double coeff = 1 / (2 * M_PI * sigma * sigma);
    double delitel = 2 * sigma * sigma;

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
    //для каждой интересной точки
    foreach(KeyFeatures::KeyPoint keyPoint, harrisPoints.keyPoints)
    {
        pointCount++;
        int     x = keyPoint.x,
                y = keyPoint.y;
        Descriptor pointDescriptor(basketNum,histogramGridSize,descriptorSize,x, y);
        int currHist = 0; // number of current histogram

        //qDebug() << "----Descriptor of " << pointCount << "point";
        //рассчитываем N*N гистограмм
        for(int ih=-descriptorRadius; ih<descriptorRadius; ih++)
        {
            for(int jh=-descriptorRadius; jh<descriptorRadius; jh++)
            {

                //для каждого пикселя в сетке
                for(int ig=-histogramRadius; ig<histogramRadius; ig++)
                {
                    for(int jg=-histogramRadius; jg<histogramRadius; jg++)
                    {
                        int     currX=x+jh*histogramGridSize+jg, //текущий х в сетке
                                currY=y+ih*histogramGridSize+ig; //текущий у в сетке

                        double currDirection = gradientDirection.GetValue(currX,currY);
                        currDirection = (currDirection < 0) ? currDirection + 360 : currDirection;
                        currDirection = (currDirection >= 360) ? currDirection - 360 : currDirection;
                        double currMagnitude = gradientMagnitude.GetValue(currX,currY);

                        double basketBetw = currDirection / basketSize;

                        int basket1 = floor(basketBetw);
                        double b1Weight = 1;

                        int basket2 = ceil(basketBetw);
                        double b2Weight = 0;

                        //если на границе с 0 или 360 градусов
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
        pointDescriptor.NormalizeDescriptor();
        imageDescriptor.append(pointDescriptor);
    }
    return imageDescriptor;

}



//GrayScaleMatrix DescriptorWorker::GetImageWithDescriptors(
//            GrayScaleMatrix inputGSMatrix, QVector<Descriptor> descriptors)
//{

//}


GrayScaleMatrix DescriptorWorker::MergeTwoMatrix(GrayScaleMatrix inputGSMatrix1, GrayScaleMatrix inputGSMatrix2)
{
    GrayScaleMatrix resultGSMatrix(inputGSMatrix1.GetWidth()+inputGSMatrix2.GetWidth(),
                                   std::max(inputGSMatrix1.GetHeight(), inputGSMatrix2.GetHeight()));

    for(int i=0; i<inputGSMatrix1.GetHeight(); i++)
    {
        for(int j=0; j< inputGSMatrix1.GetWidth(); j++)
        {
            resultGSMatrix.SetValue(j,i,inputGSMatrix1.GetValue(j,i));
        }
    }
    for(int i=0; i<inputGSMatrix2.GetHeight(); i++)
    {
        for(int j=0; j< inputGSMatrix2.GetWidth(); j++)
        {
            resultGSMatrix.SetValue(j+inputGSMatrix1.GetWidth(),i,inputGSMatrix2.GetValue(j,i));
        }
    }
    return resultGSMatrix;
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

    //Мин Макс и среднее
    for(int i = 0; i < descriptors1.size(); i++)
    {
        QVector<double> distRow;
        for(int j = 0; j < descriptors2.size(); j++)
        {
            double dist = GetDistance(descriptors1[i], descriptors2[j]);
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
    //Поиск соответствий
    for(int i = 0; i < descriptors1.size(); i++){
        double firstMinValue = std::numeric_limits<double>::max();
        int firstMinIndex = 0;
        double secondMinValue = std::numeric_limits<double>::max();
        //int secondMinIndex = 0;

        //ищем соотв. точку второй картинки
        for(int j = 0; j < descriptors2.size(); j++){
            double dist = distMatrix[i][j];
            if(dist < firstMinValue){
                secondMinValue = firstMinValue;
                //secondMinIndex = firstMinIndex;

                firstMinValue = dist;
                firstMinIndex = j;
            } else {
                if(dist < secondMinValue){
                    secondMinValue = dist;
                    //secondMinIndex = j;
                }
            }
        }

        //а теперь то же самое, но смотрим в обратную сторону. Чтобы для точки со второй картинки не было двух кандидатов из первой
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
        //берем точку если у нее NNDR < 0.6 (борьба с многозначностью). Также отбрасываем "ложные срабатывания"
        if(firstMinValue / secondMinValue < 0.6 && firstMinValue2 / secondMinValue2 < 0.6 && firstMinValue < middleValue * 0.1){

//            int r = rand() % 256;
//            int g = rand() % 256;
//            int b = rand() % 256;
//            painter.setPen(QColor(255,255,255, 170));

            for(int ii=-1; ii<=1; ii++)
            {
                for(int jj=-1; jj<=1; jj++)
                {
                    if(ii==0 || jj==0)
                    {
                        painter.drawPoint(descriptors1[i].GetX()+ii, descriptors1[i].GetY()+jj);
                        painter.drawPoint(descriptors2[firstMinIndex].GetX() + w +jj, descriptors2[firstMinIndex].GetY()+jj);

                    }
                }
            }


            painter.drawLine(descriptors1[i].GetX(), descriptors1[i].GetY(), descriptors2[firstMinIndex].GetX() + w, descriptors2[firstMinIndex].GetY());
        }
    }

//сохраняем
    return resultImage;

}
