#include <fstream>
#include <list>
#include "Matrix.h"
#include "Problem.h"


int main(int argc, char *argv[]) {
    std::string input(argv[1]), output(argv[2]);
    enum Method solvingMethod = BAND_GAUSSIAN_ELIMINATION;

    // Selección del método de resolución
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

    // Abrir el archivo de entrada
    std::ifstream handle(input, std::ifstream::in);
    BDouble width, height, h;
    unsigned amount;

    handle >> width >> height >> h >> amount;

    // Levantamos las posiciones de las sanguijuelas
    std::list<Leech> leeches = std::list<Leech>();

    unsigned int i;
    for (i = 0; i < amount; i++) {
        Leech l = Leech();
        handle >> l.x >> l.y >> l.radio >> l.temperature;
        leeches.push_back(l);
    }

    handle.close();

    // Resolvemos el problema
    Problem solver(solvingMethod, width, height, h, leeches);
    Matrix temperatures = solver.run();

    // Imprimimos la salida del método
    std::ofstream out_handle(output, std::ostream::out);

    for (std::size_t i = 0; i < temperatures.rows(); ++i) {
        for (std::size_t j = 0; j < temperatures.columns(); ++j) {
            out_handle << i << " " << j << " " << temperatures(i, j) << std::endl;
        }
    }

    out_handle.close();

    return 0;
}

/* PROBAR SHERMAN MORRISON
 * 
int main(int argc, char *argv[]) {
    Matrix A(3,3);
	std::size_t N = std::min(A.rows(), A.columns());
	A(0,0) = 2.0;
	A(0,1) = 4.0;
	A(0,2) = -2.0;
	A(1,0) = 4.0;
	// Original 9.0
	A(1,1) = 18.0;
	A(1,2) = -3.0;
	A(2,0) = -2.0;
	A(2,1) = -3.0;
	A(2,2) = 7.0;
	
	BDouble b[3] = {2.0, 4.0, 8.0};
	
	std::pair<Matrix*, Matrix*> LU = A.LU_factorization(b);
	
	Matrix& L = *LU.first;
	Matrix& U = *LU.second;
	
	std::cout << "A: " << std::endl;
	std::cout << A << std::endl;
	std::cout << "L: " << std::endl;
	std::cout << L << std::endl;
	std::cout << "U: " << std::endl;
	std::cout << U << std::endl;
	
	BDouble* y;
	std::pair<BDouble *, enum Solutions> solution; 
	solution = L.forward_substitution(b);
	y = solution.first;
		
	//Finally we calculate x = y + z * k
	for (std::size_t h = 0; h < N; h++) {
		std::cout << y[h] << " ";
    }
    std::cout << std::endl; 
    
    solution = U.backward_substitution(y);
	BDouble* x = solution.first;
	
	//Finally we calculate x = y + z * k
	for (std::size_t h = 0; h < N; h++) {
		std::cout << x[h] << " ";
    }
    std::cout << std::endl; 
    
    std::cout << "A: " << std::endl;
	std::cout << A << std::endl;
	std::cout << "Mosified (" << 1 << ", " << 1 << ")" << std::endl;
    solution = A.sherman_morrison(&L, &U, 1, 1, 9.0, b);
    x = solution.first;
    for (std::size_t h = 0; h < N; h++) {
		std::cout << x[h] << " ";
    }
    std::cout << std::endl; 
	
	return 0;
}*/

/*PROBAR LU y forward sub
 * 
int main(int argc, char *argv[]) {
    Matrix A(3,3);
	std::size_t N = std::min(A.rows(), A.columns());
	A(0,0) = 2.0;
	A(0,1) = 4.0;
	A(0,2) = -2.0;
	A(1,0) = 4.0;
	// Original 9.0
	A(1,1) = 9.0;
	A(1,2) = -3.0;
	A(2,0) = -2.0;
	A(2,1) = -3.0;
	A(2,2) = 7.0;
	
	BDouble b[3] = {2.0, 4.0, 8.0};
	
	std::pair<Matrix*, Matrix*> LU = A.LU_factorization(b);
	
	Matrix& L = *LU.first;
	Matrix& U = *LU.second;
	
	std::cout << "A: " << std::endl;
	std::cout << A << std::endl;
	std::cout << "L: " << std::endl;
	std::cout << L << std::endl;
	std::cout << "U: " << std::endl;
	std::cout << U << std::endl;
	
	BDouble* y;
	std::pair<BDouble *, enum Solutions> solution; 
	solution = L.forward_substitution(b);
	y = solution.first;
		
	//Finally we calculate x = y + z * k
	for (std::size_t h = 0; h < N; h++) {
		std::cout << y[h] << " ";
    }
    std::cout << std::endl; 
    
    Matrix B(3,3, 1, 0);
    B(0,0) = 1.0;
    B(1,0) = 1.0;
	B(1,1) = 1.0;
	B(2,1) = 1.0;
	B(2,2) = 1.0;
	std::cout << "B: " << std::endl;
	std::cout << B << std::endl;
	BDouble b2[3] = {7.0, 10.0, -10.0};
	solution = B.forward_substitution(b2);
	y = solution.first;
	for (std::size_t h = 0; h < N; h++) {
		std::cout << y[h] << " ";
    }
    std::cout << std::endl; 
    
	
	
	Matrix C(3,3, 0, 0);
    C(0,0) = 1.0;
	C(1,1) = 1.0;
	C(2,2) = 1.0;
	std::cout << "C: " << std::endl;
	std::cout << C << std::endl;
	solution = C.forward_substitution(b2);
	y = solution.first;
	for (std::size_t h = 0; h < N; h++) {
		std::cout << y[h] << " ";
    }
    std::cout << std::endl; 
    
	
	
	
	
	return 0;
}
*/

/*PROBAR LU y backward sub
 * 
int main(int argc, char *argv[]) {
    Matrix A(3,3);
	std::size_t N = std::min(A.rows(), A.columns());
	A(0,0) = 2.0;
	A(0,1) = 4.0;
	A(0,2) = -2.0;
	A(1,0) = 4.0;
	// Original 9.0
	A(1,1) = 9.0;
	A(1,2) = -3.0;
	A(2,0) = -2.0;
	A(2,1) = -3.0;
	A(2,2) = 7.0;
	
	BDouble b[3] = {2.0, 4.0, 8.0};
	
	std::pair<Matrix*, Matrix*> LU = A.LU_factorization(b);
	
	Matrix& L = *LU.first;
	Matrix& U = *LU.second;
	
	std::cout << "A: " << std::endl;
	std::cout << A << std::endl;
	std::cout << "L: " << std::endl;
	std::cout << L << std::endl;
	std::cout << "U: " << std::endl;
	std::cout << U << std::endl;
	
	BDouble* y;
	std::pair<BDouble *, enum Solutions> solution; 
	solution = L.forward_substitution(b);
	y = solution.first;
		
	//Finally we calculate x = y + z * k
	for (std::size_t h = 0; h < N; h++) {
		std::cout << y[h] << " ";
    }
    std::cout << std::endl; 
    
	BDouble* x;
	solution = U.backward_substitution(y);
	x = solution.first;
	//Finally we calculate x = y + z * k
	for (std::size_t h = 0; h < N; h++) {
		std::cout << x[h] << " ";
    }
    std::cout << std::endl; 
    
    Matrix B(3,3, 0, 1);
    B(0,0) = 1.0;
    B(0,1) = 1.0;
	B(1,1) = 1.0;
	B(1,2) = 1.0;
	B(2,2) = 1.0;
	std::cout << "B: " << std::endl;
	std::cout << B << std::endl;
	BDouble b2[3] = {10.0, 29.0, 19.0};
	solution = B.backward_substitution(b2);
	y = solution.first;
	for (std::size_t h = 0; h < N; h++) {
		std::cout << y[h] << " ";
    }
    std::cout << std::endl;
	
	return 0;
}*/ 
