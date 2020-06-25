#include <QtCore>
#include <lab1.cpp>
#include <lab2.cpp>
#include <lab3.cpp>
#include <lab4.cpp>
#include <lab5.cpp>
#include <lab6.cpp>

int main()
{
    std::cout<<std::endl<<std::endl<<"TEST 1"<<std::endl<<std::endl<<std::endl;

    lab6("C:\\Users\\artem\\Desktop\\art\\",    //Path to folder
         "drive0",   //img1 name
         ".jpg",     //extension1
         "drive0zoomIn",   //img2 name
         ".jpg",     //extension2
         3,          //harris window size
         2000,        //num of points
         8,          //num of baskets
         8,          //grid size for one histogram
         4           //num of histograms
         );
    lab6("C:\\Users\\artem\\Desktop\\art\\",    //Path to folder
         "book1",   //img1 name
         ".jpg",     //extension1
         "book2",   //img2 name
         ".jpg",     //extension2
         3,          //harris window size
         2000,        //num of points
         8,          //num of baskets
         8,          //grid size for one histogram
         4           //num of histograms
         );
    lab6("C:\\Users\\artem\\Desktop\\art\\",    //Path to folder
         "1-1",   //img1 name
         ".jpg",     //extension1
         "2-2",   //img2 name
         ".jpg",     //extension2
         3,          //harris window size
         2000,        //num of points
         8,          //num of baskets
         8,          //grid size for one histogram
         4           //num of histograms
         );
    lab6("C:\\Users\\artem\\Desktop\\art\\",    //Path to folder
         "book1",   //img1 name
         ".jpg",     //extension1
         "book2",   //img2 name
         ".jpg",     //extension2
         3,          //harris window size
         2000,        //num of points
         8,          //num of baskets
         8,          //grid size for one histogram
         4           //num of histograms
         );
    lab6("C:\\Users\\artem\\Desktop\\art\\",    //Path to folder
         "mult1",   //img1 name
         ".jpg",     //extension1
         "mult2",   //img2 name
         ".jpg",     //extension2
         3,          //harris window size
         2000,        //num of points
         8,          //num of baskets
         8,          //grid size for one histogram
         4           //num of histograms
         );
    lab6("C:\\Users\\artem\\Desktop\\art\\",    //Path to folder
         "drive0",   //img1 name
         ".jpg",     //extension1
         "drive0zoomIn",   //img2 name
         ".jpg",     //extension2
         3,          //harris window size
         2000,        //num of points
         8,          //num of baskets
         8,          //grid size for one histogram
         4           //num of histograms
         );




    std::cout<<"FINISH"<<std::endl;
    return 0;
}
