#ifndef DESCRIPTOR_H
#define DESCRIPTOR_H
#include <QVector>
#include "KeyFeatures.h"


class Descriptor
{
public:
    Descriptor(int basketNum, int histogramGridSize, int descriptorSize, int x, int y);
    double GetBasket(int histogramNum, int basketNum);
    QVector<double> GetHistogram(int histogramNum);
    void SetHistogram(int histogramNum, QVector<double> histogram);
    void SetBasket(int histogramNum, int basketNum, double value);
    void addValueToBasket(int histogramNum, int basketNum, double value); // прибавить значение к выбранной корзине
    int GetX();
    int GetY();
    int GetBasketNum(){return basketNum;}
    int GetHistogramGridSize(){return histogramGridSize;}
    int GetDescriptorSize(){return descriptorSize;}
    void NormalizeDescriptor();

private:
    int basketNum; //число корзин в гистограмме
    int histogramGridSize; //размерность сетки гистограммы M (сетка размером M*M)
    int descriptorSize; //размерность дескриптора N (количество гистограмм = N*N)
    int x,y; // координаты точки


    //Дескриптор<гистограмма<корзина>>
    //basketNum*descriptorSize*descriptorSize
    QVector<QVector<double>> pointDescriptor;
};

#endif // DESCRIPTOR_H
