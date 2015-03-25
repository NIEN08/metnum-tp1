//
// Created by Julian Bayardo on 3/25/15.
//

#ifndef _TP1_BOXEDDOUBLE_H_
#define _TP1_BOXEDDOUBLE_H_


class BoxedDouble {
public:
    BoxedDouble(double, double);
    // TODO: implement operations with error margin.
private:
    double x;
    double tolerance;
};


#endif //_TP1_BOXEDDOUBLE_H_
