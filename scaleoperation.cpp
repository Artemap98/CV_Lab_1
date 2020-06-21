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
