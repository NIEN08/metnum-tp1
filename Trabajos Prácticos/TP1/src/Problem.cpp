//
// Created by julian on 3/21/15.
//

#include "Problem.h"
#include <fstream>

Problem::Problem(std::string input, std::string output, enum Method method) :
        input(input), output(output), method(method), width(0), height(0), h(0.0), amount(0) {
    //std::ifstream handle(input, std::ifstream::in);
    //handle >> width >> height >> h >> amount;

    
    std::fstream fs;
    fs.open ("entrada.txt", std::fstream::in | std::fstream::out | std::fstream::app);

    double width, height; 
    unsigned int h, amount;

    width = fs.peek(); 
    fs.peek(); 
    height = fs.peek(); 
    fs.peek(); 
    h = fs.peek(); 
    fs.peek(); 
    amount = fs.peek(); 

    unsigned int positionX[amount];
    unsigned int positionY[amount];
    double radio[amount];
    double temperature[amount];

    int i = 0;
    // Espero que esto tire exceptions...
    while (i < amount) {
        double r = 0, t = 0;
        unsigned int x = 0, y = 0;
        //handle >> x >> y >> r >> t;
        // Acá hay que meter los puntos en algún lado que sea útil

        positionX[i] = fs.peek();
        fs.peek();
        positionY[i] = fs.peek();
        fs.peek();
        radio[i] = fs.peek();
        fs.peek();
        temperature[i] = fs.peek();
        fs.peek();
        i++;
    }


    switch(method){
       case BAND_GAUSSIAN_ELEMINATION:
       case LU_FACTORIZATION:
       case SIMPLE_ALGORITHM:
       case SHERMAN_MORRISON: 
    ;}

      fs.close();
     
}

int Problem::run() {

}

int main() {

}