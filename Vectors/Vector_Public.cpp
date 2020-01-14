//
// Created by user on 29.10.2019.
//

#include "Vector.h"
#include <math.h>
#include <iostream>
// Печать матрицы
void Vector :: print(){
    for(auto elem : vect)
        cout << elem << " ";
}
// Заполнение вектора случайными числами
void Vector :: fillRand(){
    for(double & i : vect){
        i = rand();
    }
}
// Изменение вектора на процент
Vector Vector :: addErr(double percent){
    for(double & i : vect)
        i += i * ((double)rand()  * percent / RAND_MAX);
    return *this;
}
// Изменение элемента
void Vector :: changeElem(double val, int pos){
    vect[pos] = val;
}
// Возвращает элемент
double Vector :: getElem(int pos) const {
    return vect[pos];
}
// Возвращает размер
int Vector :: getSize() const {
    return vect.size();
}
// Смена элементов местами
void Vector :: swapElem(int first, int second){
    swap(vect[first],vect[second]);
}
// Оператор +=
Vector& Vector :: operator += (Vector& vec){
    vec = this->plus(&vec);
    return vec;
}

Vector Vector :: operator - (){
    Vector c = *this;
    return c.mulDouble(-1);

}

Vector& Vector :: operator -= (Vector& vec){
    Vector *min_Vec = new(Vector);
    *min_Vec = -vec;
    return *this += *min_Vec;
}

Vector Vector :: operator-(Vector b){
    return *this + -b;
}

Vector Vector :: operator+(Vector b){
    Vector c = *this;
    return c += b;
}
// Преобразование вектрора а массив
double* Vector :: toDouble(){
    auto* answer = new double[getSize()];
    for(int i = 0; i < getSize(); i++)
        answer[i] = getElem(i);
    return answer;
}

bool Vector :: operator == (Vector const &B){
    if(this->getSize() != B.getSize())
        return false;
    for(int i = 0; i < this->getSize(); i++)
        if(abs(this->getElem(i) - B.getElem(i)) > 0.0001)
            return false;
    return true;
}
// Норма вектора
double Vector :: norm(){
    double maximum = -INFINITY;
    for(auto elem : vect)
        if(abs(elem) > maximum)
            maximum = abs(elem);
    return maximum;
}
