#include "KeyFeatures.h"

KeyFeatures::KeyFeatures()
{

}


GrayScaleMatrix KeyFeatures::GetMoravecMatrix(GrayScaleMatrix inputGSMatrix, int windowRadius)
{
    qDebug() << "--GetMoravecMatrix(GrayScaleMatrix inputGSMatrix, int windowRadius)";

    int w = inputGSMatrix.GetWidth();
    int h = inputGSMatrix.GetHeight();

    GrayScaleMatrix moravecMatrix(w,h);
    inputGSMatrix = Convolution::GaussianFilter(inputGSMatrix,1.3);

    for(int i=windowRadius+1; i < h-windowRadius-1; i++)
    {
        for(int j=windowRadius+1; j < w-windowRadius-1; j++)
        {
            double localCoeff = std::numeric_limits <double>::max();
            for(int iy = -1; iy <= 1; iy++)
            {
                for(int ix = -1; ix <= 1; ix++)
                {
                    if(ix!=0 && iy!=0)
                    {
                        double sum = 0;
                        for (int u = -windowRadius; u <= windowRadius; u++) {
                            for (int v = -windowRadius; v <= windowRadius; v++) {
                                double tmp =  inputGSMatrix.GetValue(j,i) - inputGSMatrix.GetValue(j+ix+u,i+iy+v);
                                sum += tmp * tmp;
                            }
                        }
                        localCoeff = std::min(sum,localCoeff);
                    }
                }
                moravecMatrix.SetValue(j,i,localCoeff);
            }
        }
    }
    return moravecMatrix;
}

GrayScaleMatrix KeyFeatures::GetResultMoravec(GrayScaleMatrix inputGSMatrix, GrayScaleMatrix moravecMatrix, int windowRadius, int resultPointsNum)
{
    qDebug() << "--GetPointsMoravec(GrayScaleMatrix inputGSMatrix, GrayScaleMatrix moravecMatrix, int windowRadius, int resultPointsNum)";
    KeyPointSet interestingPoints = GetLocalMaximums(moravecMatrix, windowRadius, false);

    int w = inputGSMatrix.GetWidth();
    int h = inputGSMatrix.GetHeight();

    interestingPoints = ReducePoints(interestingPoints, resultPointsNum, std::min(w/2,h/2));


    foreach (KeyPoint point, interestingPoints.keyPoints) {
        for(int i=-1; i<=1; i++)
        {
            for(int j=-1; j<=1; j++)
            {
                if(i==0 || j==0)
                {
                    try
                    {
                        inputGSMatrix.SetValue(point.x+j,point.y+i,1);
                    } catch(_exception e){}
                }
            }
        }
    }
    return inputGSMatrix;
}

//получить локальные максимумы
KeyFeatures::KeyPointSet KeyFeatures::GetLocalMaximums(GrayScaleMatrix inputMatrix, int windowRadius, bool isHarris)
{
    qDebug() << "--GetLocalMaximums(GrayScaleMatrix moravecMatrix, int windowRadius, bool isHarris)";
    KeyPointSet points;

    int w = inputMatrix.GetWidth();
    int h = inputMatrix.GetHeight();

    //находим мин и макс для порога
    double  min = std::numeric_limits <double>::max(),
            max = std::numeric_limits <double>::min();
    for(int i = 0; i < h; i++)
    {
         for(int j = 0; j < w; j++)
         {
             double temp = inputMatrix.GetValue(j,i);
            if (max < temp) max = temp;
            if (min > temp) min = temp;
        }
    }

    //задаем порог
    double threshold = min + (max - min) * 0.005;
    if (isHarris)
        threshold = min + (max - min) * 0.004;


    //добавляем точки в список, только если они сильнейшие в своей окрестности
    for (int i = 0; i < h; i++)
    {
        for (int j = 0; j < w; j++)
        {
            bool is_correct = true;
            double sLocal = inputMatrix.GetValue(j,i);
            for (int px = -windowRadius; px <= windowRadius && is_correct; px++) {
                for (int py = -windowRadius; py <= windowRadius && is_correct; py++) {
                    if (px != 0 || py != 0) {
                        is_correct = sLocal > inputMatrix.GetValue(j+px,i+py);
                    }
                }
            }
            if (is_correct && sLocal > threshold) {
                points.keyPoints.append(KeyPoint(j,i,sLocal));
            }
        }
    }

    return points;
}

KeyFeatures::KeyPointSet KeyFeatures::ReducePoints(KeyFeatures::KeyPointSet points, int resultPointsNum, int maxRadius)
{
    qDebug() << "--ReducePoints(KeyFeatures::KeyPointSet points, int resultPointsNum, int maxRadius)";

    QVector<QVector<double>> distanceMatrix;
    QVector<QVector<bool>> sLocalGreaterMatrix;

    double minDist = std::numeric_limits <double>::max();

    //строим матрицу дистанций между всеми точками
    //и  матрицу отношений "больше" для локальных коэффициентов
    for(int i=0; i < points.keyPoints.size(); i++)
    {
        QVector<double> distanceRow;
        QVector<bool> sLocalGreaterRow;
        for(int j= 0; j < points.keyPoints.size(); j++)
        {
            double xd = points.keyPoints[i].x - points.keyPoints[j].x;
            double yd = points.keyPoints[i].y - points.keyPoints[j].y;
            double dist = sqrt(xd * xd + yd * yd);

            if(dist < minDist) minDist = dist;

            distanceRow.append(dist);
            sLocalGreaterRow.append(points.keyPoints[i].sLocal <= points.keyPoints[j].sLocal);
        }
        distanceMatrix.append(distanceRow);
        sLocalGreaterMatrix.append(sLocalGreaterRow);
    }

    int r = std::ceil(minDist);

    //пока точек слишком много и радиус в пределах допустимого
    while (points.keyPoints.size() > resultPointsNum && r < maxRadius)
    {
        for(int i=0; i< distanceMatrix.size(); i++)
        {
            for(int j=0; j<distanceMatrix[i].size(); j++)
            {
                if(distanceMatrix[i][j] <= r && i!=j)
                {
                    if(sLocalGreaterMatrix[i][j])
                    {
                        for(int ii=0; ii < i; ii++)
                        {
                            distanceMatrix[ii].remove(i);
                            sLocalGreaterMatrix[ii].remove(i);
                        }
                        distanceMatrix.remove(i);
                        sLocalGreaterMatrix.remove(i);
                        points.keyPoints.remove(i);
                        i--;
                        break;
                    }
                }
            }
        }
        r++;
        qDebug() << "--radius = " << r << "; num of key points = " << points.keyPoints.size();
    }

//    //строим треугольную матрицу дистанций между всеми точками
//    //и треугольную матрицу отношений "больше" для локальных коэффициентов
//    QVector<int> distCount = QVector<int>(2000);
//    for(int i=0; i < points.keyPoints.size(); i++)
//    {
//        QVector<double> distanceRow;
//        QVector<bool> sLocalGreaterRow;
//        for(int j= i + 1; j < points.keyPoints.size(); j++)
//        {
//            double xd = points.keyPoints[i].x - points.keyPoints[j].x;
//            double yd = points.keyPoints[i].y - points.keyPoints[j].y;
//            double dist = sqrt(xd * xd + yd * yd);

//            if(dist < minDist) minDist = dist;

//            distCount[floor(dist)]++;
//            distanceRow.append(dist);
//            sLocalGreaterRow.append(points.keyPoints[i].sLocal <= points.keyPoints[j].sLocal);
//        }
//        distanceMatrix.append(distanceRow);
//        sLocalGreaterMatrix.append(sLocalGreaterRow);
//    }

//    int r = std::ceil(minDist);

//    //пока точек слишком много и радиус в пределах допустимого
//    while (points.keyPoints.size() > resultPointsNum && r < maxRadius)
//    {
//        for(int i=0; i< distanceMatrix.size(); i++)
//        {
//            for(int j=0; j<distanceMatrix[i].size(); j++)
//            {
//                if(distanceMatrix[i][j] <= r)
//                {
//                    if(sLocalGreaterMatrix[i][j])
//                    {
//                        for(int ii=0; ii < i; ii++)
//                        {
//                            distanceMatrix[ii].remove(i-ii);
//                            sLocalGreaterMatrix[ii].remove(i-ii);
//                        }
//                        distanceMatrix.remove(i);
//                        sLocalGreaterMatrix.remove(i);
//                        points.keyPoints.remove(i);
//                        i--;
//                        break;
//                    }
////                    else
////                    {
////                        int jLocal = i+j+1;
////                        for(int ii=0; ii < jLocal; ii++)
////                        {
////                            distanceMatrix[ii].remove(jLocal-1-ii);
////                            sLocalGreaterMatrix[ii].remove(jLocal-1-ii);
////                        }
////                        distanceMatrix.remove(jLocal);
////                        sLocalGreaterMatrix.remove(jLocal);
////                        points.keyPoints.remove(jLocal);
////                        j--;
////                    }
//                }
//            }
//        }
//        r++;
//        qDebug() << "--radius = " << r << "; num of key points = " << points.keyPoints.size();
//    }
    return points;
}


GrayScaleMatrix KeyFeatures::GetHarrisMatrix(GrayScaleMatrix inputGSMatrix, int windowRadius)
{
    qDebug() << "--GetHarrisMatrix(GrayScaleMatrix inputGSMatrix, int windowRadius)";

    inputGSMatrix = Convolution::GaussianFilter(inputGSMatrix, 1.3);

    //находим производные
    GrayScaleMatrix derivateX = Convolution::DerivateX(inputGSMatrix);
    GrayScaleMatrix derivateY = Convolution::DerivateY(inputGSMatrix);

    int w = inputGSMatrix.GetWidth();
    int h = inputGSMatrix.GetHeight();

    QVector<QVector<double>>    a,
                                b,
                                c;

    //находим веса для окна - ядро Гаусса
    double sigma = static_cast<double>(windowRadius*2) / 6;
    QVector<QVector<double>>gaussKernel;

    double coeff = 1 / (2 * M_PI * sigma * sigma);
    double delitel = 2 * sigma * sigma;

    for (int u = -windowRadius; u <= windowRadius; u++)
    {
        QVector<double> gaussRow;
        for (int v = -windowRadius; v <= windowRadius; v++)
        {
            gaussRow.append(coeff * exp(- (u * u + v * v) / delitel));
        }
        gaussKernel.append(gaussRow);
    }

    //Вычисляем A, B, C для всех точек
    for (int i = 0; i < h; i++)
    {
        QVector<double> aRow,bRow,cRow;
         for (int j = 0; j < w; j++)
         {
            double sumA = 0, sumB = 0, sumC = 0;

            for (int u = -windowRadius; u <= windowRadius; u++)
            {
                for (int v = -windowRadius; v <= windowRadius; v++)
                {
                    double i_x = derivateX.GetValue(j + v, i + u);
                    double i_y = derivateY.GetValue(j + v, i + u);
                    sumA += i_x * i_x * gaussKernel[windowRadius+u][windowRadius+v];
                    sumB += i_x * i_y * gaussKernel[windowRadius+u][windowRadius+v];
                    sumC += i_y * i_y * gaussKernel[windowRadius+u][windowRadius+v];
                }
            }
            aRow.append(sumA);
            bRow.append(sumB);
            cRow.append(sumC);
        }
         a.append(aRow);
         b.append(bRow);
         c.append(cRow);
    }


    GrayScaleMatrix harrisMatrix(w,h);   //здесь будем хранить значения оператора

    for (int i =0; i < h; i++) {
        for (int j =0; j < w; j++) {
            double sc = a[i][j] + c[i][j];
            double d = a[i][j] * c[i][j] - b[i][j] * b[i][j];
            double det = sc * sc - 4 * d;
            double L1 = (sc + sqrt(det)) / 2;
            double L2 = (sc - sqrt(det)) / 2;
            double cHarris = std::min (L1, L2);
            harrisMatrix.SetValue(j, i, cHarris);
        }
    }

    return harrisMatrix;
}


KeyFeatures::KeyPointSet KeyFeatures::GetPointsHarris(GrayScaleMatrix inputGSMatrix, GrayScaleMatrix harrisMatrix, int windowRadius, int resultPointsNum)
{
    qDebug() << "--GetPointsHarris(GrayScaleMatrix inputGSMatrix, GrayScaleMatrix moravecMatrix, int windowRadius, int resultPointsNum)";
    KeyPointSet interestingPoints = GetLocalMaximums(harrisMatrix, 2, true);

    int w = inputGSMatrix.GetWidth();
    int h = inputGSMatrix.GetHeight();

    interestingPoints = ReducePoints(interestingPoints, resultPointsNum, std::min(w/2,h/2));

    return interestingPoints;
}


GrayScaleMatrix KeyFeatures::GetResultHarrisMatrix(GrayScaleMatrix inputGSMatrix, KeyPointSet interestingPoints)
 {
    foreach (KeyPoint point, interestingPoints.keyPoints)
    {
        for(int i=-1; i<=1; i++)
        {
            for(int j=-1; j<=1; j++)
            {
                if(i==0 || j==0)
                {
                    try
                    {
                        inputGSMatrix.SetValue(point.x+j,point.y+i,1);
                    } catch(_exception e){}
                }
            }
        }
    }
    return inputGSMatrix;
}

