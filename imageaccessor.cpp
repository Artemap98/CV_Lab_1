#include "imageaccessor.h"


ImageAccessor::ImageAccessor()
{   
}

QVector<QVector<unsigned char>> ImageAccessor::GetMatrixFromImage(QString path)
{
    QImage inputImage;
    QVector<QVector<unsigned char>> imageMatrix;
    if(!inputImage.load(path))
    {
        inputImage = inputImage.convertToFormat(QImage::Format_RGB32);
        {
            for(int i = 0; i < inputImage.width(); i++)
            {
                QVector<unsigned char> pixelLine;

                for(int j = 0; j < inputImage.height(); j++)
                {
                    float brightness;
                    QColor pColor = inputImage.pixelColor(i,j);
                    brightness = pColor.red()*0.299 + pColor.green()*0.587 + pColor.blue()*0.114;//0.213𝑅𝑅+0.715𝐺𝐺+0.072𝐵𝐵;

                    pixelLine.append(brightness);
                }
                imageMatrix.append(pixelLine);
            }
        }
    }

    return imageMatrix;
}

bool ImageAccessor::DrawImageFromMatrix(QVector<QVector<unsigned char>> imageMatrix, QString path)
{
    QImage outputImage(imageMatrix.at(0).size(),imageMatrix.size(),QImage::Format_RGB32);
    for(int i=0; i<outputImage.width(); i++)
    {
        for(int j=0; j<outputImage.height(); j++)
        {
            unsigned char grayShade = imageMatrix.at(i).at(j);
            outputImage.setPixel(i,j,qRgb(grayShade,grayShade,grayShade));
        }
    }
    outputImage.save(path);
    return true;
}
