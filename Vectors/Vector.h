//
// Created by user on 12.11.2019.
//

#ifndef VECTOR_DO_NOT_CHANGE_VECTOR_H
#define VECTOR_DO_NOT_CHANGE_VECTOR_H

#include <vector>
using namespace std;

class Vector{
private:
    vector <double> vect;

    Vector plus(Vector* b);     // Сложение векторов

    Vector mulDouble(double x); // Умножение вектора на число

public:
    Vector() = default;

    explicit Vector(int size){
        vect.resize(size, 0);
    }

    explicit Vector(vector <double> data){
        vect = std::move(data);
    }
    Vector(double* data, int size){
        for(int i = 0; i < size; i++)
            vect.push_back(data[i]);
    }
    void print(); // Печать матрицы

    void fillRand(); // Заполнение вектора случайными числами

    Vector addErr(double percent = 0.05); // Изменение вектора на процент

    void changeElem(double val, int pos); // Изменение элемента

    double getElem(int pos) const; // Возвращает элемент

    int getSize() const; // Возвращает размер

    void swapElem(int first, int second); // Смена элементов местами

    Vector& operator += (Vector& vec); // Оператор +=

    Vector operator - ();

    Vector& operator -= (Vector& vec);

    Vector operator-(Vector b);

    Vector operator+(Vector b);

    bool operator == (Vector const &b);

    double* toDouble(); // Преобразование вектрора а массив

    double norm(); // Норма вектора

};


#endif //VECTOR_DO_NOT_CHANGE_VECTOR_H
