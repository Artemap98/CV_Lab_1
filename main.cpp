#include <QtCore>
#include <lab1.cpp>
#include <lab2.cpp>
#include <lab3.cpp>
#include <lab4.cpp>

int main()
{
//    Path to folder, image name, sigma for Gauss filter
//    lab1("C:\\Users\\artem\\Desktop\\art\\","VYH01gMwyHw.jpg",1);

//    lab1("C:\\Users\\artem\\Desktop\\art\\","MKWSqkFW3Z0.jpg",2);


    lab2("C:\\Users\\artem\\Desktop\\art\\",    //Path to folder
         "white",   //img name
         ".jpg",    //extension
         6,         //num of octaves
         2,         //num of layers
         0,         //sigmaA
         1,         //sigma0
         6.66       //sigmaL for L(x,y,sigma)
         );

//    lab3("C:\\Users\\artem\\Desktop\\art\\",    //Path to folder
//        "drive1",   //img name
//        ".jpg",     //extension
//        2,          //operator window radius
//        400         //num of key points
//        );

//    qDebug() << QString::fromUtf8("русская строка");
//    std::cout<<std::endl<<"test 1"<<std::endl<<std::endl;

//    lab4("C:\\Users\\artem\\Desktop\\art\\",    //Path to folder
//         "drive1",   //img1 name
//         ".jpg",     //extension1
//         "drive2",   //img2 name
//         ".jpg",     //extension2
//         2,
//         400,
//         8,
//         16,
//         4
//         );

//    std::cout<<std::endl<<"test 2"<<std::endl<<std::endl;

//    lab4("C:\\Users\\artem\\Desktop\\art\\",    //Path to folder
//         "drive1",   //img1 name
//         ".jpg",     //extension1
//         "drive3",   //img2 name
//         ".jpg",     //extension2
//         2,
//         400,
//         8,
//         16,
//         4
//         );

//    std::cout<<std::endl<<"test 3"<<std::endl<<std::endl;

//    lab4("C:\\Users\\artem\\Desktop\\art\\",    //Path to folder
//         "drive1",   //img1 name
//         ".jpg",     //extension1
//         "drive4",   //img2 name
//         ".jpg",     //extension2
//         2,
//         400,
//         8,
//         16,
//         4
//         );


    std::cout<<"FINISH"<<std::endl;
    return 0;
}
