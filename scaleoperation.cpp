#include "KeyFeatures.h"
#include "scaleoperation.h"
#include <QDebug>

ScaleOperation::ScaleOperation()
{

}

GrayScaleMatrix ScaleOperation::Downsample(GrayScaleMatrix inputMatrix)
{
    int     outW = inputMatrix.GetWidth()/2,
            outH = inputMatrix.GetHeight()/2;
    GrayScaleMatrix outputMatrix(outW,outH);

    for(int i=0; i < outH; i++)
    {
        for(int j=0; j < outW; j++)
        {
            outputMatrix.SetValue(j,i,inputMatrix.GetValue(j*2,i*2));
        }
    }
    return outputMatrix;
}


ScaleOperation::Pyramid ScaleOperation::GetPyramid(GrayScaleMatrix inputMatrix, int numOfOctaves, int numOfLayers, double sigmaA, double sigma0)
{
    Pyramid outputPyramid;
    double scaleInterval = pow(2.0, 1.0 / numOfLayers); //интервал между масштабами
    GrayScaleMatrix outputMatrix(inputMatrix);
    double sigmaB = sqrt(sigma0 * sigma0 - sigmaA * sigmaA);    //с каким сигма нужно сгладить, чтобы получить с требуемым sigma0
    outputMatrix = Convolution::GaussianFilter(inputMatrix,sigmaB);

    double actSigma = sigma0;//действительная сигма для октав
    GrayScaleMatrix currMatrix = outputMatrix;
    for(int i = 0; i < numOfOctaves; i++)
    {
        double currSigma = sigma0;//локальная сигма в октаве

        Octave currOctave;   //создаем новую октаву
        Layer currLayer(currMatrix,currSigma,actSigma);//нулевой слой, равный уменьшенному последнему из предыдущей октавы
        currOctave.layers.append(currLayer);
        for (int j = 1; j <= numOfLayers; j++)
        {
            double sigmaKoeff = currSigma;
            currSigma *= scaleInterval;
            actSigma *= scaleInterval;
            //находим коэфф., с которым нужно сгладить, чтобы получить требуемую сигму
            sigmaKoeff = sqrt(currSigma*currSigma - sigmaKoeff*sigmaKoeff);
            qDebug()<<"scaleInterval="<<scaleInterval<<"; currSigma="<< currSigma<<"; sigmaKoeff" << sigmaKoeff;
            currMatrix = Convolution::GaussianFilter(currMatrix, sigmaKoeff);
            currLayer.matrix = currMatrix;
            currLayer.currentSigma = currSigma;
            currLayer.actualSigma = actSigma;
            currOctave.layers.append(currLayer);
        }
        outputPyramid.octaves.append(currOctave);

        if (i < numOfOctaves - 1)
            currMatrix =  Downsample(currMatrix);
    }
    return outputPyramid;
}



ScaleOperation::Pyramid ScaleOperation::GetPyramidWithOverlay(GrayScaleMatrix inputMatrix,
                                                   int numOfOctaves, int numOfLayers,
                                                   double sigmaA, double sigma0,
                                                   int overlay// default = 3
                                                   )
{
    Pyramid outputPyramid;
    outputPyramid.overlay = overlay;
    double scaleInterval = pow(2.0, 1.0 / numOfLayers); //интервал между масштабами
    GrayScaleMatrix outputMatrix(inputMatrix);
    double sigmaB = sqrt(sigma0 * sigma0 - sigmaA * sigmaA);    //с каким сигма нужно сгладить, чтобы получить с требуемым sigma0
    outputMatrix = Convolution::GaussianFilter(inputMatrix,sigmaB);

    double actSigma = sigma0;//действительная сигма для октав
    GrayScaleMatrix currMatrix = outputMatrix;
    GrayScaleMatrix nextOctaveMatrix = currMatrix;
    for(int i = 0; i < numOfOctaves; i++)
    {
        double currSigma = sigma0;//локальная сигма в октаве

        Octave currOctave;   //создаем новую октаву
        Layer currLayer(currMatrix,currSigma,actSigma);//нулевой слой, равный уменьшенному последнему из предыдущей октавы
        currOctave.layers.append(currLayer);

        for (int j = 1; j <= numOfLayers + overlay; j++)
        {
            double sigmaKoeff = currSigma;
            currSigma *= scaleInterval;
            actSigma *= scaleInterval;
            //находим коэфф., с которым нужно сгладить, чтобы получить требуемую сигму
            sigmaKoeff = sqrt(currSigma*currSigma - sigmaKoeff*sigmaKoeff);
            qDebug()<<"scaleInterval="<<scaleInterval<<"; currSigma="<< currSigma<<"; sigmaKoeff" << sigmaKoeff;
            currMatrix = Convolution::GaussianFilter(currMatrix, sigmaKoeff);
            currLayer.matrix = currMatrix;
            currLayer.currentSigma = currSigma;
            currLayer.actualSigma = actSigma;
            currOctave.layers.append(currLayer);
            if(j == numOfLayers)
                nextOctaveMatrix = currMatrix;
        }
        outputPyramid.octaves.append(currOctave);

        if (i < numOfOctaves - 1)
            currMatrix =  Downsample(nextOctaveMatrix);
    }
    return outputPyramid;
}


//ScaleOperation::Pyramid ScaleOperation::GetDoGPyramid(ScaleOperation::Pyramid inputPyramid, int harrisRadius)
//{

//    Pyramid DoGPyramid;
//    //для каждого уровня в каждой октаве находим значения Difference of Gaussian
//    for(int i = 0; i < inputPyramid.octaves.size(); i++)
//    {
//        Octave currOctave;
//        currOctave.level = i;
//        //идем до предпоследнего элемента, т.к. из текущего слоя вычитаем значение следующего
//        for(int j = 0; j < inputPyramid.octaves[i].layers.size()-1; j++)
//        {
//            Layer currLayer = inputPyramid.octaves[i].layers[j];
//            Layer nextLayer = inputPyramid.octaves[i].layers[j];
////            currLayer.actualSigma = sqrt(currLayer.actualSigma * currLayer.actualSigma
////                                         + nextLayer.actualSigma*nextLayer.actualSigma);
////            currLayer.currentSigma = sqrt(currLayer.currentSigma * currLayer.currentSigma
////                                         + nextLayer.currentSigma * nextLayer.currentSigma);

//            for(int ii = 0; ii < currLayer.matrix.GetWidth(); ii++)
//            {
//                for(int jj = 0; jj < currLayer.matrix.GetHeight(); jj++)
//                {
//                    currLayer.matrix.SetValue(ii,jj, currLayer.matrix.GetValue(ii,jj)-currLayer.matrix.GetValue(ii,jj));
//                }
//            }
//            currOctave.layers.append(currLayer);
//        }
//        DoGPyramid.octaves.append(currOctave);
//    }
//    //теперь ищем макс отклик в 3д
//    for(int i=0; i < DoGPyramid.octaves.size(); i++)
//    {
//        //отступаем от начала и конца, чтобы иметь слой сверху и снизу.
//        for(int j=1; j < DoGPyramid.octaves[i].layers.size()-1; j++)
//        {
//            Layer currLayer = inputPyramid.octaves[i].layers[j];
//            for(int ii = 0; ii < currLayer.matrix.GetWidth(); ii++)
//            {
//                for(int jj = 0; jj < currLayer.matrix.GetHeight(); jj++)
//                {
//                    double responseVal = currLayer.matrix.GetValue(ii,jj);// значение проверяемого пикселя
//                    if(responseVal > 0.03)//отсекаем силшком маленькие отклики
//                    {
//                        bool isMinimal = true, isMaximal = true, canContinue = true;
//                        //для трех слоев октавы
//                        for(int k = -1; k <= 1 && canContinue; k++)
//                        {
//                            //для каждого слоя окрестность вокруг точки
//                            for(int zi = -1; zi <= 1 && canContinue; zi++)
//                            {
//                                for(int zj = -1; zj <= 1 && canContinue; zj++)
//                                {
//                                    //если текущая точка окрестности не является проверяемой
//                                    if(!(k != 0 && zi != 0 && zj != 0))
//                                    {
//                                        //текущее значение окрестности
//                                        double tempVal = DoGPyramid.octaves[i].layers[j+k].matrix.GetValue(ii+zi,jj+zj);
//                                        //если значение в окрестности больше проверяемого,
//                                        //то преверяемое не может быть максимальным в области
//                                        if(tempVal > responseVal)
//                                            isMaximal = false;
//                                        //если значение в окрестности меньше проверяемого,
//                                        //то преверяемое не может быть минимальным в области
//                                        if(tempVal < responseVal)
//                                            isMinimal = false;

//                                        //если проверяемое не максимальное и не минимальное,
//                                        //то дальше проверять нет смысла
//                                        if(!(isMinimal && isMaximal))
//                                        {
//                                            canContinue = false;
//                                        }
//                                    }
//                                }
//                            }
//                        }
//                        //если точка минимальная или максимальная в окрестности
//                        if(isMaximal || isMinimal)
//                        {
//                            //создаем матрицу окрестности вокруг точки
//                            GrayScaleMatrix pointLocality(harrisRadius*2+1,harrisRadius*2+1);
//                            for(int di = -harrisRadius; di <= harrisRadius; di++)
//                            {
//                                for(int dj = -harrisRadius; dj <= harrisRadius; dj++)
//                                {
//                                    pointLocality.SetValue(di,dj,DoGPyramid.octaves[i].layers[j].matrix.GetValue(ii+di,jj+dj));
//                                }

//                            }
//                            //получаем отклик Харриса для проверяемой точки
//                            pointLocality = KeyFeatures::GetHarrisMatrix(pointLocality, harrisRadius);

//                            if(pointLocality.GetValue(harrisRadius,harrisRadius) > 0.01)//отсекаем слабый отклик
//                            {

//                            }
//                        }
//                    }
//                }
//            }
//        }
//    }
//}


double ScaleOperation::GetL(ScaleOperation::Pyramid inputPyramid, int x, int y, double sigma)
{
    Layer targetLayer = inputPyramid.octaves[0].layers[0];
    int     octaveLevel = 0,
            octaveCount = 0;
    //находим октаву и уровень, на котором сигма будет наиболее близка к искомой
    foreach(Octave currOctave, inputPyramid.octaves)
    {
        foreach(Layer currLayer, currOctave.layers)
        {
            if(fabs(currLayer.actualSigma - sigma) < fabs(targetLayer.actualSigma - sigma))
            {
                targetLayer = currLayer;
                octaveLevel = octaveCount;
            }
        }
        octaveCount++;
    }

    //учитываем, что на следующих октавах изображения меньше по размеру
    int ynew = static_cast<int>(y / pow(2., octaveLevel));
    int xnew = static_cast<int>(x / pow(2., octaveLevel));
    if(xnew >= targetLayer.matrix.GetWidth()) xnew = targetLayer.matrix.GetWidth()-1;
    if(ynew >= targetLayer.matrix.GetHeight()) ynew = targetLayer.matrix.GetHeight()-1;
    return targetLayer.matrix.GetValue(xnew,ynew);
}
