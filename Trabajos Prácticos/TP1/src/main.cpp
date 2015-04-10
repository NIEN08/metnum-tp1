#include <sys/time.h>
#include <fstream>
#include <limits>
#include <list>
#include "Matrix.h"
#include "Problem.h"

int main(int argc, char *argv[]) {
	enum Method solvingMethod = BAND_GAUSSIAN_ELIMINATION;

	// Selección del método de resolución
	if ((*argv[3]) == '1') {
		solvingMethod = LU_FACTORIZATION;
	} else if ((*argv[3]) == '2') {
		solvingMethod = SIMPLE_ALGORITHM;
	} else if ((*argv[3]) == '3') {
		solvingMethod = SHERMAN_MORRISON;
	} else if ((*argv[3]) != '0') {
		std::cerr << "Invalid method parameter" << std::endl;
		return -1;
	}

	// Abrir el archivo de entrada
	std::ifstream handle(argv[1], std::ios::in);

	BDouble width, height, h;
	unsigned amount;

	handle >> width >> height >> h >> amount;

	// Levantamos las posiciones de las sanguijuelas
	std::list<Leech> leeches = std::list<Leech>();

	for (int i = 0; i < amount; i++) {
		Leech l = Leech();
		handle >> l.x >> l.y >> l.radio >> l.temperature;
		leeches.push_back(l);
	}

	handle.close();

	// Medimos tiempo de inicio
	timeval start, end;
	gettimeofday(&start, NULL);

	// Resolvemos el problema
	Problem solver(solvingMethod, width, height, h, leeches);
	Matrix temperatures = solver.run();

	// Medimos tiempo de finalización, e imprimimos el delta.
	gettimeofday(&end,NULL);
	std::cout.precision(std::numeric_limits<double>::digits10);
	std::cout << std::fixed << (1000000*(end.tv_sec-start.tv_sec)+(end.tv_usec-start.tv_usec))/1000000.0 << std::endl;

	// Imprimimos la salida del método
	std::ofstream out_handle(argv[2], std::ostream::out | std::ofstream::trunc);

	std::cout.precision(5);
	for (int i = 0; i < temperatures.rows(); ++i) {
		for (int j = 0; j < temperatures.columns(); ++j) {
			out_handle << i << '\t' << j << '\t' << std::fixed << temperatures(i, j) << std::endl;
		}
	}

	out_handle.close();

	return 0;
}