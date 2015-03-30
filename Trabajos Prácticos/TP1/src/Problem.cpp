//
// Created by julian on 3/21/15.
//

#include "Problem.h"
#include <fstream>
#include <iostream>
//using namespace std;

Problem::Problem(std::string input, std::string output, enum Method method) :
        input("entrada2.txt"), output("salida.txt"), method(method), width(0), height(0), h(0.0), amount(0) {
    //std::ifstream handle(input, std::ifstream::in);
    

    
    std::fstream fs;
    fs.open ("entrada.txt");

    double width, height; 
    unsigned int h, amount;

    fs >> width >> height >> h >> amount;

    unsigned int positionX[amount];
    unsigned int positionY[amount];
    double radio[amount];
    double temperature[amount];


    int i = 0;
    // Espero que esto tire exceptions...
    while (i < amount) {
        // Acá hay que meter los puntos en algún lado que sea útil
        fs >> positionX[i] >> positionY[i] >> radio[i] >> temperature[i]; 
        i++;
    }

    double matriz[(int)height+1][(int)width+1];    
    // la idea es guardar los resultados en esta matriz 
    // para despues imprimir el archivo de salida con los resultados

    switch(method){
       case BAND_GAUSSIAN_ELEMINATION:
       break;
       case LU_FACTORIZATION:
       break;
       case SIMPLE_ALGORITHM:
       break;
       case SHERMAN_MORRISON:
       break; 
    ;}

    
    // esta parte imprime la salida en un archivo
    // ANDA BIEN
    std::fstream outputFile;        
    outputFile.open("salida.txt");

    for(int i = 0; i < height; i++){
        for(int j = 0; j < width; j++){
            outputFile << i << " " << j << " " << matriz[i][j] << std::endl;
        }
    }
    
    fs.close();
    outputFile.close();
     
}

int Problem::run() {

}

int main() {
    //Problem p("entrada2.txt", "salida.txt", BAND_GAUSSIAN_ELEMINATION);
}