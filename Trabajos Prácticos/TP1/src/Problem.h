#ifndef _TP1_PROBLEM_H_
#define _TP1_PROBLEM_H_ 1

#include "Matrix.h"
#include <vector>
#include <tuple>
#include <algorithm>

enum Method {
    BAND_GAUSSIAN_ELIMINATION,
    LU_FACTORIZATION,
    SIMPLE_ALGORITHM,
    SHERMAN_MORRISON
};

// extern vectorInfo;

class Problem {
public:
    Problem(enum Method method, const Matrix &input)
            : temperatures(input), method(method) { }

    Matrix run() {
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



struct info { 
  BDouble posX;
  BDouble posY;
  BDouble radio;
  BDouble temperature;
  BDouble dif;
};

bool myfunction (info i,info j) { return (i.dif < j.dif); }

std::tuple<bool, info> simple_algotithm(std::vector<info> myvector){
    
    // ordeno el vector de menor a mayor con respecto a la variable dif.
    // dif = sqrt[(x-x1)^2 + (y-y1)^2]   
    // dif es la distancia de un punto con respecto al punto critico.
    sort (myvector.begin(), myvector.end(), myfunction);

    info result;
    std::tuple<bool, info> res(false, result);

    if(myvector.size() == 0){
        return res;
    }

    result.posX = myvector[0].posX;
    result.posY = myvector[0].posY;
    result.radio = myvector[0].radio;
    result.temperature = myvector[0].temperature;
    result.dif = myvector[0].dif;

    for(int i = 1; i < myvector.size(); i++){
        
        if(myvector[i].temperature > result.temperature){
            result.posX = myvector[i].posX;
            result.posY = myvector[i].posY;
            result.radio = myvector[i].radio;
            result.temperature = myvector[i].temperature;
            result.dif = myvector[i].dif;
        }
    }

    std::tuple<bool, info> res2(true, result);
    return res2;
}



private:
    Matrix temperatures;
    enum Method method;
};


#endif //_TP1_PROBLEM_H_
