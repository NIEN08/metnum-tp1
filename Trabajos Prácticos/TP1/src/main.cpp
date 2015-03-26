#include "Matrix.h"
#include <iostream>
#include <istream>

int main(){
	BandMatrix<int, 0> m(1, 1, 3, 3);
	BandMatrix<int, 0> a(1, 1, 3, 3);
	//printf("%u\n", a==m);
	/*m(0,0) = 1;
	m(1,1) = 1;
	m(2,2) = 1;*/
}
