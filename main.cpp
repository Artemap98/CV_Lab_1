#include <QtCore>
#include <lab1.cpp>
#include <lab2.cpp>
#include <lab3.cpp>
#include <lab4.cpp>
#include <lab5.cpp>

int main()
{
//    std::cout<<std::endl<<"lab 1"<<std::endl<<std::endl;
//    std::cout<<std::endl<<"test 1"<<std::endl<<std::endl;
//    lab1("C:\\Users\\artem\\Desktop\\art\\",
//         "lena.jpg",
//         1  //sigma
//         );

//    std::cout<<std::endl<<"test 2"<<std::endl<<std::endl;
//    lab1("C:\\Users\\artem\\Desktop\\art\\",
//         "luffy.png",
//         2  //sigma
//         );

//    std::cout<<std::endl<<"test 3"<<std::endl<<std::endl;
//    lab1("C:\\Users\\artem\\Desktop\\art\\",
//         "drive1.jpg",
//         3  //sigma
//         );

//    std::cout<<std::endl<<"lab 2"<<std::endl<<std::endl;
//    std::cout<<std::endl<<"test 1"<<std::endl<<std::endl;
//    lab2("C:\\Users\\artem\\Desktop\\art\\",    //Path to folder
//         "lena",   //img name
//         ".jpg",    //extension
//         6,         //num of octaves
//         2,         //num of layers
//         0,         //sigmaA
//         1,         //sigma0
//         6.66       //sigmaL for L(x,y,sigma)
//         );

//    std::cout<<std::endl<<"test 2"<<std::endl<<std::endl;
//    lab2("C:\\Users\\artem\\Desktop\\art\\",    //Path to folder
//         "luffy",   //img name
//         ".png",    //extension
//         4,         //num of octaves
//         4,         //num of layers
//         0,         //sigmaA
//         1,         //sigma0
//         4.20       //sigmaL for L(x,y,sigma)
//         );

//    std::cout<<std::endl<<"test 3"<<std::endl<<std::endl;
//    lab2("C:\\Users\\artem\\Desktop\\art\\",    //Path to folder
//         "drive1",   //img name
//         ".jpg",    //extension
//         3,         //num of octaves
//         3,         //num of layers
//         0,         //sigmaA
//         1,         //sigma0
//         2.28       //sigmaL for L(x,y,sigma)
//         );

//    std::cout<<std::endl<<"test 3"<<std::endl<<std::endl;
//    lab2("C:\\Users\\artem\\Desktop\\art\\",    //Path to folder
//         "drive11",   //img name
//         ".jpg",    //extension
//         3,         //num of octaves
//         10,         //num of layers
//         0,         //sigmaA
//         1,         //sigma0
//         2.28       //sigmaL for L(x,y,sigma)
//         );

    std::cout<<std::endl<<"lab 3"<<std::endl<<std::endl;
    std::cout<<std::endl<<"test 1"<<std::endl<<std::endl;
    lab3("C:\\Users\\artem\\Desktop\\art\\",    //Path to folder
        "lena",   //img name
        ".jpg",     //extension
        1,          //operator window radius
        150000        //num of key points
        );

    lab3("C:\\Users\\artem\\Desktop\\art\\",    //Path to folder
        "lenacopy",   //img name
        ".jpg",     //extension
        1,          //operator window radius
        1500       //num of key points
        );

    std::cout<<std::endl<<"test 2"<<std::endl<<std::endl;
    lab3("C:\\Users\\artem\\Desktop\\art\\",    //Path to folder
        "luffy",   //img name
        ".png",     //extension
        2,          //operator window radius
        100000        //num of key points
        );

    lab3("C:\\Users\\artem\\Desktop\\art\\",    //Path to folder
        "luffycopy",   //img name
        ".png",     //extension
        2,          //operator window radius
        1000        //num of key points
        );

    std::cout<<std::endl<<"test 3"<<std::endl<<std::endl;
    lab3("C:\\Users\\artem\\Desktop\\art\\",    //Path to folder
        "drive1",   //img name
        ".jpg",     //extension
        3,          //operator window radius
        50000        //num of key points
        );

    lab3("C:\\Users\\artem\\Desktop\\art\\",    //Path to folder
        "drive1copy",   //img name
        ".jpg",     //extension
        3,          //operator window radius
        500       //num of key points
        );

//    std::cout<<std::endl<<"lab 4"<<std::endl<<std::endl;
//    std::cout<<std::endl<<"test 1"<<std::endl<<std::endl;
//    lab4("C:\\Users\\artem\\Desktop\\art\\",    //Path to folder
//         "lena",   //img1 name
//         ".jpg",     //extension1
//         "lena2",   //img2 name
//         ".jpg",     //extension2
//         3,          //harris window size
//         800,        //num of points
//         8,          //num of baskets
//         4,          //grid size for one histogram
//         4           //num of histograms
//         );

//    std::cout<<std::endl<<"test 2"<<std::endl<<std::endl;
//    lab4("C:\\Users\\artem\\Desktop\\art\\",    //Path to folder
//         "luffy",   //img1 name
//         ".png",     //extension1
//         "luffy2",   //img2 name
//         ".png",     //extension2
//         2,          //harris window size
//         1600,        //num of points
//         8,          //num of baskets
//         4,          //grid size for one histogram
//         4           //num of histograms
//         );

//    std::cout<<std::endl<<"test 3"<<std::endl<<std::endl;
//    lab4("C:\\Users\\artem\\Desktop\\art\\",    //Path to folder
//         "drive1",   //img1 name
//         ".jpg",     //extension1
//         "drive2",   //img2 name
//         ".jpg",     //extension2
//         4,          //harris window size
//         400,        //num of points
//         8,          //num of baskets
//         4,          //grid size for one histogram
//         4           //num of histograms
//         );

//    std::cout<<std::endl<<"lab 5"<<std::endl<<std::endl;
//    std::cout<<std::endl<<"test 0"<<std::endl<<std::endl;
//    lab5("C:\\Users\\artem\\Desktop\\art\\",    //Path to folder
//         "lena",   //img1 name
//         ".jpg",     //extension1
//         "lena",   //img2 name
//         ".jpg",     //extension2
//         3,          //harris window size
//         800,        //num of points
//         8,          //num of baskets
//         4,          //grid size for one histogram
//         4           //num of histograms
//         );


//    std::cout<<std::endl<<"test 1"<<std::endl<<std::endl;
//    lab5("C:\\Users\\artem\\Desktop\\art\\",    //Path to folder
//         "drive1",   //img1 name
//         ".jpg",     //extension1
//         "drive5",   //img2 name
//         ".jpg",     //extension2
//         3,          //harris window size
//         400,        //num of points
//         8,          //num of baskets
//         16,          //grid size for one histogram
//         4           //num of histograms
//         );

//    std::cout<<std::endl<<"test 2"<<std::endl<<std::endl;
//    lab5("C:\\Users\\artem\\Desktop\\art\\",    //Path to folder
//         "lena",   //img1 name
//         ".jpg",     //extension1
//         "lena17",   //img2 name
//         ".jpg",     //extension2
//         3,          //harris window size
//         400,        //num of points
//         8,          //num of baskets
//         8,          //grid size for one histogram
//         4           //num of histograms
//         );

//    std::cout<<std::endl<<"test 3"<<std::endl<<std::endl;
//    lab5("C:\\Users\\artem\\Desktop\\art\\",    //Path to folder
//         "lena",   //img1 name
//         ".jpg",     //extension1
//         "lena-42",   //img2 name
//         ".jpg",     //extension2
//         3,          //harris window size
//         400,        //num of points
//         8,          //num of baskets
//         8,          //grid size for one histogram
//         4           //num of histograms
//         );

//    std::cout<<std::endl<<"test 4"<<std::endl<<std::endl;
//    lab5("C:\\Users\\artem\\Desktop\\art\\",    //Path to folder
//         "lena",   //img1 name
//         ".jpg",     //extension1
//         "lena-55",   //img2 name
//         ".jpg",     //extension2
//         3,          //harris window size
//         400,        //num of points
//         8,          //num of baskets
//         8,          //grid size for one histogram
//         4           //num of histograms
//         );

//    std::cout<<std::endl<<"test 5"<<std::endl<<std::endl;
//    lab5("C:\\Users\\artem\\Desktop\\art\\",    //Path to folder
//         "lena",   //img1 name
//         ".jpg",     //extension1
//         "lena-90",   //img2 name
//         ".jpg",     //extension2
//         3,          //harris window size
//         400,        //num of points
//         8,          //num of baskets
//         8,          //grid size for one histogram
//         4           //num of histograms
//         );

//    std::cout<<std::endl<<"test 5"<<std::endl<<std::endl;
//    lab5("C:\\Users\\artem\\Desktop\\art\\",    //Path to folder
//         "drive1",   //img1 name
//         ".jpg",     //extension1
//         "drive180",   //img2 name
//         ".jpg",     //extension2
//         3,          //harris window size
//         400,        //num of points
//         8,          //num of baskets
//         16,          //grid size for one histogram
//         4           //num of histograms
//         );




//    std::cout<<std::endl<<"test 2"<<std::endl<<std::endl;
//    lab5("C:\\Users\\artem\\Desktop\\art\\",    //Path to folder
//         "lena",   //img1 name
//         ".jpg",     //extension1
//         "lena5",   //img2 name
//         ".jpg",     //extension2
//         3,          //harris window size
//         800,        //num of points
//         8,          //num of baskets
//         16,          //grid size for one histogram
//         4           //num of histograms
//         );


//    std::cout<<std::endl<<"test 3"<<std::endl<<std::endl;
//    lab5("C:\\Users\\artem\\Desktop\\art\\",    //Path to folder
//         "luffy",   //img1 name
//         ".png",     //extension1
//         "luffy5",   //img2 name
//         ".png",     //extension2
//         3,          //harris window size
//         600,        //num of points
//         8,          //num of baskets
//         4,          //grid size for one histogram
//         4           //num of histograms
//         );



    std::cout<<"FINISH"<<std::endl;
    return 0;
}
