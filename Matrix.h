//
// Created by user on 29.10.2019.
//

#ifndef MATRIXDLL_MATRIX_H
#define MATRIXDLL_MATRIX_H

#include <vector>
#include "Vectors/Vector.h"
#include <iostream>
#include <ctime>
#include <random>
using namespace std;

typedef struct{
    Vector res;
    Vector iterMax;
    int iterNum{};
}JacobiAns_t;


class Matrix{
protected:
    vector <vector<double>> matrix;

    bool isRight(); // Является ли матрица право-треугольной

    bool isLeft(); // Является ли матрица лево-треугольной

    bool swapNULL(int line, int column); // Если элемент нулевой, то ищет ненулевой элемент ниже него и меняет местами, если такого не нашлось, то возвращает false

    void lineDiv(int line, double divider); // Деление всей строки на число

    void columnDiv(int column, double divider); // Деление всей столбца на число

    void columnDiv(int column, double divider, int startLine); // Деление элементов столбца снизу от данного (начиная с него же) на число

    void lineDiv(int line, double divider, int startColumn); // Деление элементов строки справа от данного (начиная с него же) на число

    void minLines(int line1, int line2, double multiplier); // Вычитание из второй строки первой, умноженной на число

    void minLines(int line1, int line2, double multiplier, int startMin); // Вычитание из второй строки первой, умноженной на число, начиная с элемента (начиная с него же)

    int findMaxInColumn(int column); // Поиск макимума в столбце

    void swapLine(int i, int j); // Смена строк местами

    void swapColumn(int i, int j); // Смена столбцов местами

    double trace(); // След матрицы

    pair <Matrix, Matrix>  LU_Decomposition(Vector& b); // LU разложение

    void make_Right_Triangle(); // Сделать матрицу право-треугольной

    Matrix mul(Matrix* A); // Умножение матриц

    Vector matrixMulVector(Vector const & b); // Умножение матрицы на вектор

    double GaussMove(); // Прямой ход Гаусса

    long double GaussMove(Matrix* rightMatrix); // Прямой ход Гаусса с правой частью

    void GaussRevMove(Matrix* rightMatrix); // Обратный ход Гаусса с правой частью

    pair<Matrix, Vector> alpha_betta_Create(Vector* b); // Создание матрицы и вектора для метода простых итераций

    pair <int,int> findMax(); // Поиск максимального элемента матрицы вне диагонали

public:
    int addNum;
    int subNum;
    int multNum;
    int divNum;

    Matrix()= default;

    explicit Matrix(int size){ // Создание 0 матрицы с заданным размером
        matrix.resize(size, {});
        for(int i = 0; i < size; i ++)
            matrix[i].resize(size, 0);
        addNum=0;
        subNum=0;
        multNum=0;
        divNum=0;
    }

    Matrix(double* data, int size){ // создание матрицы с данными
        matrix.resize(size, vector<double>(size));
        for (int i = 0; i < size; i++)
            for (int j = 0; j < size; j++)
                matrix[i][j] = data[j * size + i];
            addNum=0;
            subNum=0;
            multNum=0;
            divNum=0;
    }

    explicit Matrix(vector<vector<double>> data){ // Создание матрицы по вектору
        matrix = std::move(data);
        addNum=0;
        subNum=0;
        multNum=0;
        divNum=0;
    }

    void hilbert(int size); // Генерация матрицы Гильберта

    Matrix inverse(); // Обратная матрица

    double determinate(); // Определитель

    void fill(int size); // Заполнение матрицы из консоли

    bool isIdentity(); // Проверка, является ли матрица единичной

    void print() const; // Печать матрицы

    void reFill(double** data, int size); // Заполенение матрицы из массива double

    void randFill(int size); // Заполнение матрицы случайными числами

    int getSize() const; // Возврашает размер матрицы

    void changeVal(int line, int column, double newVal); // Изменение элемента матрицы

    double getVal(int line, int column) const; // Возвращает элемент матрицы

    void setIdentity(); // Преобразование матрицы в единичную

    void transpose(); // Транспонирование

    double norm(); // Норма матрицы

    double cond(); // Число обусловленности

    Vector& operator *= (Vector const& vec); // *= оператор

    Vector operator*(Vector const& vec);

    Matrix& operator *= (Matrix A);

    Matrix operator*(const Matrix& B);

    bool operator == (Matrix const &B);

    Vector FindSolution_LU(double* b); // Поиск решения LU разложением

    Vector FindSolution_MSI(double * b, double epsilon); // Поиск решения методом простых итераций

    void GenerateDiagWithCond(int condDegree); // Генерирование диагональной матрицы с числом обусловленности

    void GenerateWithSmallDet(); //Генерирование матрицы с определителем, близким к нулю

    pair<Vector, vector<double>> FindSolution_MSI(double* b, double* x, double epsilon); // Поиск решения методом простых итераций с данными

    pair<JacobiAns_t,Matrix> JacobiRotateMethod(); // Поиск с.ч. и с.в. методом Якоби

    double* toDouble(); // Преобразование матрицы в массив

    bool isDiag(); // Проверка, является ли матрица диагональной

    void rotate(double phi, int i, int j);

    static double RND(){
        random_device rd;
        mt19937 generator(rd());
        return generator();
    }
    static double RND(double a, double b){
        random_device rd;
        mt19937 generator(rd());
        uniform_real_distribution<double> dis(a, b);
        return dis(generator);
    }

    pair<JacobiAns_t,Matrix> JacobiRotateMethod(double* Solution);

    pair<Vector, vector<double>>  FindSolution_MSI_withLimits(double* b_data, double* x, double epsilon, int operationsNum);
};


#endif //MATRIXDLL_MATRIX_H
