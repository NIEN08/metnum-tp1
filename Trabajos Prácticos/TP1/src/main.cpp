#include "Matrix.h"
#include <fstream>

int main(int argc, char *argv[]) {

    unsigned int x,y,h,amount;
    double r,t,width,height;

    std::fstream outputFile;
    outputFile.open("entrada.txt");
    outputFile >> width >> height >> h >> amount;
    std::cout << width << " " << height << " " << h << " " << amount << std::endl;
    
    for(int i = 0; i < amount; i++){
            outputFile >> x >> y >> r >> t;
            std::cout << x << " " << y << " " << r << " " << t << " " << std::endl; 
    }

    return 0;
}
