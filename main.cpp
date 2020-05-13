#include <QtCore>
#include <lab1.cpp>
#include <lab2.cpp>
#include <lab3.cpp>

int main()
{
//    Path to folder, image name, sigma for Gauss filter
//    lab1("C:\\Users\\artem\\Desktop\\art\\","VYH01gMwyHw.jpg",1);
//    lab1("C:\\Users\\artem\\Desktop\\art\\","MKWSqkFW3Z0.jpg",2);


//    lab2("C:\\Users\\artem\\Desktop\\art\\",    //Path to folder
//         "drive1",   //img name
//         ".jpg",    //extension
//         6,         //num of octaves
//         2,         //num of layers
//         0,         //sigmaA
//         1,         //sigma0
//         6.66       //sigmaL for L(x,y,sigma)
//         );

    lab3("C:\\Users\\artem\\Desktop\\art\\",    //Path to folder
        "Luffy",   //img name
        ".png",     //extension
        2,          //operator window radius
        1000         //num of key points
        );


    return 0;
}
