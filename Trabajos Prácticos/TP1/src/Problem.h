#ifndef _TP1_PROBLEM_H_
#define _TP1_PROBLEM_H_ 1

#include <list>
#include "Matrix.h"

enum Method {
    BAND_GAUSSIAN_ELIMINATION,
    LU_FACTORIZATION,
    SIMPLE_ALGORITHM,
    SHERMAN_MORRISON
};

typedef struct _Leech {
    BDouble x;
    BDouble y;
    BDouble radio;
    BDouble temperature;
    BDouble tempPC;
} Leech;

class Problem {
public:
    Problem(enum Method method,
            const BDouble &width,
            const BDouble &height,
            const BDouble &h,
            const std::list<Leech> &leeches)
            : width(width), height(height), h(h),
              xCoordinates(std::round(height / h)), yCoordinates(std::round(width / h)),
              leeches(leeches), method(method), temperatures(xCoordinates + 1, yCoordinates + 1) {

        for (std::list<Leech>::Iterator b = leeches.begin(); b != leeches.end(); ++leeches) {
            // Distribuimos las temperaturas de la sanguijuela
            int topX = std::floor(b.x + b.radio);
            int bottomX = std::ceil(b.x - b.radio);
            int topY = std::floor(b.y + b.radio);
            int bottomY = std::ceil(b.y - b.radio);

            // Ponemos las coordenadas en rango
            topX = std::min(std::max(topX, 0), xCoordinates);
            bottomX = std::min(std::max(bottomX, 0), xCoordinates);

            topY = std::min(std::max(topY, 0), yCoordinates);
            bottomY = std::min(std::max(bottomY, 0), yCoordinates);

            // Seteamos las temperaturas en la matriz.
            // Cabe destacar, la temperatura de cada sanguijuela es igual para todos los puntos que cubre.
            for (int l = bottomX; l <= topX; ++l) {
                for (int j = bottomY; j <= topY; ++j) {
                    // Sólo queda la temperatura más alta
                    if (temperatures(l, j) < t) {
                        temperatures(l, j) = t;
                    }
                }
            }
        }

        // Ponemos los bordes en -100.0C
        for(i = 0; i <= xCoordinates; ++i){
            if (i == 0 || i == xCoordinates) {
                for (std::size_t j = 0; j <= yCoordinates; ++j) {
                    temperatures(i, j) = -100.0;
                }
            } else {
                temperatures(i, 0) = -100.0;
                temperatures(i, yCoordinates) = -100.0;
            }
        }
    }

    Matrix run() {
        // Variamos nuestra solución según el método
        switch (method) {
            case BAND_GAUSSIAN_ELIMINATION:
                break;
            case LU_FACTORIZATION:
                break;
            case SIMPLE_ALGORITHM:
                break;
            case SHERMAN_MORRISON:
                break;
        }
    }
    
    
    std::tuple<bool, Leech> simple_algotithm(std::list<Leech> leeches){
   
      Leech result;
      std::tuple<bool, Leech> res(false, result);
      if(leeches.size() == 0){
	return res;
      }
      
      /*result.x = leeches.front().x;
      result.y = leeches.front().y;
      result.radio = leeches.front().radio;
      result.temperature = leeches.front().temperature;*/
      result.tempPC = 9999999999999;
     
      
      for (std::list<Leech>::Iterator a = leeches.begin(); a != leeches.end(); ++leeches) {
	
	    for (std::list<Leech>::Iterator b = leeches.begin(); b != leeches.end(); ++leeches) {
		  
		  if(a != b){
			// Distribuimos las temperaturas de la sanguijuela
			int topX = std::floor(b.x + b.radio);
			int bottomX = std::ceil(b.x - b.radio);
			int topY = std::floor(b.y + b.radio);
			int bottomY = std::ceil(b.y - b.radio);

			// Ponemos las coordenadas en rango
			topX = std::min(std::max(topX, 0), xCoordinates);
			bottomX = std::min(std::max(bottomX, 0), xCoordinates);

			topY = std::min(std::max(topY, 0), yCoordinates);
			bottomY = std::min(std::max(bottomY, 0), yCoordinates);

			// Seteamos las temperaturas en la matriz.
			// Cabe destacar, la temperatura de cada sanguijuela es igual para todos los puntos que cubre.
			for (int l = bottomX; l <= topX; ++l) {
			    for (int j = bottomY; j <= topY; ++j) {
				    // Sólo queda la temperatura más alta
				    if (temperatures(l, j) < t) {
				      temperatures(l, j) = t;
				    }
			    }
			}
		  }

	      }

	      
	      // Ponemos los bordes en -100.0C
	      for(i = 0; i <= xCoordinates; ++i){
		  if (i == 0 || i == xCoordinates) {
		      for (std::size_t j = 0; j <= yCoordinates; ++j) {
			  temperatures(i, j) = -100.0;
		      }
		  } else {
		      temperatures(i, 0) = -100.0;
		      temperatures(i, yCoordinates) = -100.0;
		  }
	      }
	      
	      if(temperatures(width/2, height/2) < 235){
		    
		    // guardo la temperatura del punto critico sin a, eso me sirve para saber cual es el que genera menor temperatura en el punto critico.
		    a.tempPC = temperatures(width/2, height/2);
		   
		    // me fijo si la temperatura del punto critico de la sanguijuela a es la que hace que el punto critico tenga la menor temperatura.
		    if(result.tempPC > a.tempPC){
			result.x = a.x;
			result.y = a.y;
			result.radio = a.radio;
			result.temperature = a.temperature;
			result.tempPC = a.tempPC;
		    }
	      }
      }     
      
      std::tuple<bool, Leech> res2(true, result);
      return res2;
      
    }
    
    
private:
    int xCoordinates, yCoordinates;
    BDouble width, height, h;
    std::list<Leech> leeches;
    Matrix temperatures;
    enum Method method;
};


#endif //_TP1_PROBLEM_H_
