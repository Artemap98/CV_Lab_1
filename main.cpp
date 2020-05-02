#include <QtCore>
#include <lab1.cpp>
/*
Лабораторная работа 1. Свертка
-Написать основу для представления изображений и их обработки свертками
-Реализовать вычисление частных производных и оператора Собеля
-Реализовать фильтр Гаусса
-Реализовать отображение полученных результатов
NB! Нормирование выходных данных
*/


int main()
{
    //Path to folder, image name, sigma for Gauss filter
    lab1("C:\\Users\\artem\\Desktop\\art\\pics\\","greetings-from-berlin.jpg",3);
    lab1("C:\\Users\\artem\\Desktop\\art\\pics\\","chicken.jpg",5);
    lab1("C:\\Users\\artem\\Desktop\\art\\pics\\","rabbit-on-a-train-detail.jpg",7);


    return 0;
}
