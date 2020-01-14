//
// Created by user on 12.11.2019.
//

#include "Vector.h"
// Сложение векторов
Vector Vector :: plus(Vector* b){
    if(b->getSize() != this->getSize())
        return Vector();
    Vector c(b->getSize());
    for(int i = 0; i < this->getSize(); i++)
        c.changeElem(this->getElem(i) + b->getElem(i), i);
    return c;
}
// Умножение вектора на число
Vector Vector :: mulDouble(double x){
    Vector c = *this;
    for(int i = 0; i < this->getSize(); i++)
        c.changeElem(c.getElem(i) * x, i);
    return c;
}