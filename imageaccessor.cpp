#include "imageaccessor.h"


GrayScaleMatrix ImageAccessor::GetMatrixFromImage(QString path)
{
    QImage inputImage1;
    QVector<QVector<unsigned char>> imageMatrix;
    if(inputImage1.load(path))
    {
        QImage inputImage = inputImage1.convertToFormat(QImage::Format_RGB32);
        {
            for(int i = 0; i < inputImage.height(); i++)
            {
                QVector<unsigned char> pixelLine;

                for(int j = 0; j < inputImage.width(); j++)
                {
                    float brightness;
                    QColor pColor = inputImage.pixelColor(j,i);
                    brightness = pColor.red()*0.299 + pColor.green()*0.587 + pColor.blue()*0.114;//0.213𝑅𝑅+0.715𝐺𝐺+0.072𝐵𝐵;

                    pixelLine.append(brightness);
                }
                imageMatrix.append(pixelLine);
            }
        }
    }
    GrayScaleMatrix gsMatrix(imageMatrix);
    return gsMatrix;
}

void ImageAccessor::DrawImageFromMatrix(GrayScaleMatrix gsMatrix, QString path)
{
    QVector<QVector<unsigned char>> imageMatrix = gsMatrix.GetMatrix255();
    QImage outputImage(imageMatrix.at(0).size(),imageMatrix.size(),QImage::Format_RGB32);


    for(int i=0; i<outputImage.height(); i++)
    {
        for(int j=0; j<outputImage.width(); j++)
        {
            unsigned char grayShade = imageMatrix[i][j];
            outputImage.setPixel(j,i,qRgb(grayShade,grayShade,grayShade));
        }
    }
    outputImage.save(path);
}

QImage ImageAccessor::GetImageFromMatrix(GrayScaleMatrix gsMatrix)
{
    QVector<QVector<unsigned char>> imageMatrix = gsMatrix.GetMatrix255();
    QImage outputImage(imageMatrix.at(0).size(),imageMatrix.size(),QImage::Format_RGB32);


    for(int i=0; i<outputImage.height(); i++)
    {
        for(int j=0; j<outputImage.width(); j++)
        {
            unsigned char grayShade = imageMatrix[i][j];
            outputImage.setPixel(j,i,qRgb(grayShade,grayShade,grayShade));
        }
    }
    return outputImage;
}
