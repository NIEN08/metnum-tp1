#include "Problem.h"
#include <fstream>
#include "Matrix.h"

struct info { 
  BDouble posX;
  BDouble posY;
  BDouble radio;
  BDouble temperature;
  BDouble dif;
};

std::vector<info> vectorInfo;



int main(int argc, char *argv[]) {
    std::string input(argv[1]), output(argv[2]);
    enum Method solvingMethod = BAND_GAUSSIAN_ELIMINATION;


    if ((*argv[3]) == '1') {
        solvingMethod = LU_FACTORIZATION;
    } else if ((*argv[3]) == '2') {
        solvingMethod = SIMPLE_ALGORITHM;
    } else if ((*argv[3]) == '3') {
        solvingMethod = SHERMAN_MORRISON;
    } else {
        std::cout << "Unknown method to solve problem" << std::endl;
        return -1;
    }

    std::ifstream handle(input, std::ifstream::in);
    BDouble width, height, h;
    unsigned amount;

    handle >> width >> height >> h >> amount;

    int xCoordinates = std::round(width / h);
    int yCoordinates = std::round(height / h);

    // Creamos la matriz de temperaturas
    Matrix temperatures = Matrix(xCoordinates + 1, yCoordinates + 1);

    // Levantamos las posiciones de las sanguijuelas
    unsigned int i;
    for (i = 0; i < amount; i++) {
        BDouble x, y, r, t;
        handle >> x >> y >> r >> t;

        info aux;
        aux.posX = x;
        aux.posY = y;
        aux.radio = r;
        aux.temperature = t;
        // calculo la diferencia que hay entre el punto (x,y) y el punto critico.
        aux.dif = sqrt(pow((double)(width/2) - x, 2) + pow((double)(height/2) - y, 2));

        // Me fijo si el punto critico esta dentro del radio de la sanguijuela i.
        if(aux.dif < aux.radio){
            vectorInfo.push_back(aux);
        }


        // Distribuimos las temperaturas de las sanguijuelas
        int topX = std::floor(x + r);
        int bottomX = std::ceil(x - r);
        int topY = std::floor(y + r);
        int bottomY = std::ceil(y - r);

        // Ponemos las coordenadas en rango
        topX = std::min(std::max(topX, 0), xCoordinates);
        bottomX = std::min(std::max(bottomX, 0), xCoordinates);

        topY = std::min(std::max(topY, 0), yCoordinates);
        bottomY = std::min(std::max(bottomY, 0), yCoordinates);

        // Seteamos las temperaturas en la matriz.
        // Cabe destacar, la temperatura de cada sanguijuela es igual para todos los puntos que cubre.
        for (int l = bottomX; l <= topX; ++l) {
            for (int j = bottomY; j <= topY; ++j) {
                if (temperatures(l, j) < t) {
                    temperatures(l, j) = t;
                }
            }
        }
    }

    handle.close();

    // Ponemos los bordes en -100.0C
    for(unsigned int i = 0; i <= xCoordinates; ++i){
        if (i == 0 || i == xCoordinates) {
            for (std::size_t j = 0; j <= yCoordinates; ++j) {
                temperatures(i, j) = -100.0;
            }
        } else {
            temperatures(i, 0) = -100.0;
            temperatures(i, yCoordinates) = -100.0;
        }
    }

    Problem solver(solvingMethod, temperatures);
    temperatures = solver.run();

    std::ofstream out_handle(output, std::ofstream::out);

    for (std::size_t i = 0; i < temperatures.rows(); ++i) {
        for (std::size_t j = 0; j < temperatures.columns(); ++j) {
            out_handle << i << " " << j << " " << temperatures(i, j) << std::endl;
        }
    }

    out_handle.close();
    return 0;
}
/*
int main(int argc, char *argv[]) {
    Matrix A(3,3);

	A(0,0) = 2.0;
	A(0,1) = 4.0;
	A(0,2) = -2.0;
	A(1,0) = 4.0;
	A(1,1) = 9.0;
	A(1,2) = -3.0;
	A(2,0) = -2.0;
	A(2,1) = -3.0;
	A(2,2) = 7.0;
	
	BDouble b[3] = {2.0, 8.0, 10.0};
	
	std::pair<Matrix*, Matrix*> LU = A.LU_factorization(b);
	
	Matrix& L = *LU.first;
	Matrix& U = *LU.second;
	
	std::cout << "L: " << std::endl;
	std::cout << L << std::endl;
	std::cout << "U: " << std::endl;
	std::cout << U << std::endl;
	std::cout << "A: " << std::endl;
	std::cout << A << std::endl;

	
	return 0;
}

*/