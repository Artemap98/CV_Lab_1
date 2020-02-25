#ifndef IMAGEACCESSOR_H
#define IMAGEACCESSOR_H
#include <QImage>
#include <QVector>

class ImageAccessor
{
public:
    ImageAccessor();
    QVector<QVector<unsigned char>> GetMatrixFromImage(QString path);
    bool DrawImageFromMatrix(QVector<QVector<unsigned char>> imageMatrix, QString path);
};

#endif // IMAGEACCESSOR_H
