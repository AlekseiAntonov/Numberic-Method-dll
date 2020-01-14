//
// Created by user on 29.10.2019.
//
#include "Matrix.h"
// Является ли матрица право-треугольной
bool Matrix :: isRight(){
    for(int i = 1; i < this->getSize(); i++)
        for(int j = 0; j < i; j++)
            if(matrix[i][j])
                return false;
    return true;
}
// Является ли матрица лево-треугольной
bool Matrix :: isLeft(){
    for(int i = 0; i < getSize() - 1; i++)
        for(int j = i + 1; j < getSize(); j++)
            if(matrix[i][j])
                return false;
    return true;
}
// Если элемент нулевой, то ищет ненулевой элемент ниже него и меняет местами, если такого не нашлось, то возвращает false
bool Matrix :: swapNULL(int line, int column){
    if(matrix[line][column])
        return true;
    for(int i = line + 1; i < getSize(); i++)
        if(matrix[i][column]){
            swapLine(i, line);
            return true;
        }
    return false;
}
// Деление всей строки на число
void Matrix :: lineDiv(int line, double divider){
    for(int i = 0; i < getSize(); i++)
        matrix[line][i] /= divider;
}
// Деление всей столбца на число
void Matrix :: columnDiv(int column, double divider){
    for(int i = 0; i < getSize(); i++)
        matrix[i][column] /= divider;
}
// Деление элементов столбца снизу от данного (начиная с него же) на число
void Matrix :: columnDiv(int column, double divider, int startLine){
    for(int i = startLine; i < getSize(); i++)
        matrix[i][column] /= divider;
}
// Деление элементов строки справа от данного (начиная с него же) на число
void Matrix :: lineDiv(int line, double divider, int startColumn){
    for(int i = startColumn; i < getSize(); i++) {
        matrix[line][i] /= divider;
        divNum++;
    }
}
// Вычитание из второй строки первой, умноженной на число
void Matrix :: minLines(int line1, int line2, double multiplier){
    for(int i = 0; i < getSize(); i++) {
        matrix[line2][i] -= matrix[line1][i] * multiplier;
        subNum++;
        multNum++;
    }
}
// Вычитание из второй строки первой, умноженной на число, начиная с элемента (начиная с него же)
void Matrix :: minLines(int line1, int line2, double multiplier, int startMin){
    for(int i = startMin; i < getSize(); i++) {
        matrix[line2][i] -= matrix[line1][i] * multiplier;
        subNum++;
        multNum++;
    }
}
// Поиск макимума в столбце
int Matrix :: findMaxInColumn(int column) {
    int max = column;
    for (int i = column; i < this->getSize(); i++)
        if (this->getVal(max, max) < this->getVal(i, column))
            max = i;
    return max;
}
// Смена строк местами
void Matrix :: swapLine(int i, int j){
    swap(matrix[i], matrix[j]);
}
// Смена столбцов местами
void Matrix :: swapColumn(int i, int j){
    for(int k = 0; k < this->getSize(); k++)
        swap(matrix[k][i], matrix[k][j]);
}
// След матрицы
double Matrix :: trace(){
    double answer = 0;
    for(int i = 0; i < getSize(); i++)
        answer += matrix[i][i];
    return answer;
}
// LU разложение
pair <Matrix, Matrix>  Matrix :: LU_Decomposition(Vector& b) {
    int size = this->getSize();
    Matrix U = *this;
    Matrix L(size);
    for(int i = 0; i < getSize(); i++) {
        if(U.getVal(i,i) == 0)
            if (!swapNULL(i, i))
                continue;

        int max = U.findMaxInColumn(i);
        U.swapLine(i, max);
        L.swapLine(i,max);
        b.swapElem(i, max);

        double divider = U.getVal(i, i);

        L.changeVal(i, i, 1);
        for (int j = i + 1; j < getSize(); j++) {
            if (U.getVal(j,i) == 0)
                continue;
            double multiplier = U.getVal(j, i) / divider;
            divNum++;
            L.changeVal(j, i, multiplier);
            U.minLines(i, j, multiplier);
        }
    }
    addNum += L.addNum + U.addNum;
    subNum += L.subNum + U.subNum;
    divNum += L.divNum + U.divNum;
    multNum += L.multNum + U.multNum;
    return make_pair(L, U);
}
// Сделать матрицу право-треугольной
void Matrix :: make_Right_Triangle(){
    for(unsigned int i = 1; i < matrix.size(); i++)
        for(unsigned int j = 0; j < i ; j++)
            matrix[i][j] = 0;
}
// Умножение матриц
Matrix Matrix :: mul(Matrix* A){
    if(A->getSize() != this->getSize())
        return Matrix();
    Matrix C(A->getSize());
    for(int i = 0; i < this->getSize(); i++)
        for(int j = 0; j < A->getSize(); j++)
            for(int k = 0; k < this->getSize(); k++)
                C.changeVal(i, j, C.getVal(i, j) + this->getVal(i , k) * A->getVal(k, j));
    return C;
}
// Умножение матрицы на вектор
Vector Matrix :: matrixMulVector(Vector const &  b){
    if(getSize() != b.getSize())
        return {};
    Vector answer(getSize());
    for(int i = 0; i < getSize(); i++){
        for(int j = 0; j < getSize(); j++) {
            answer.changeElem(answer.getElem(i) + this->getVal(i, j) * b.getElem(j), i);
            addNum++;
            multNum++;
        }
    }
    return answer;
}
// Прямой ход Гаусса
double Matrix :: GaussMove(){
    double multiplierMatrix = 1;
    for(int i = 0; i < getSize(); i++) {
        if(!matrix[i][i]) {
            if (!swapNULL(i, i)) {
                multiplierMatrix = 0;
                continue;
            }
            multiplierMatrix *= -1;
        }
        double divider = getVal(i, i);
        multiplierMatrix *= divider;
        lineDiv(i, getVal(i, i), i);
        for (int j = i + 1; j < getSize(); j++) {
            if (!matrix[j][i])
                continue;
            double multiplier = matrix[j][i];
            minLines(i, j, multiplier);
        }
    }
    return multiplierMatrix;
}
// Прямой ход Гаусса с правой частью
long double Matrix :: GaussMove(Matrix* rightMatrix){
    long double multiplierMatrix = 1;
    for(int i = 0; i < getSize(); i++) {
        if(!matrix[i][i]) {
            if (!swapNULL(i, i)) {
                multiplierMatrix = 0;
                continue;
            }
            multiplierMatrix *= -1;
            multNum++;
        }
        double divider = getVal(i, i);
        multiplierMatrix *= divider;
        multNum++;
        lineDiv(i, divider, i);
        rightMatrix->lineDiv(i, divider);
        divNum+=rightMatrix->getSize();
        for (int j = i + 1; j < getSize(); j++) {
            if (!matrix[j][i])
                continue;
            double multiplier = matrix[j][i];
            minLines(i, j, multiplier);
            rightMatrix->minLines(i, j, multiplier);
        }
    }
    return multiplierMatrix;
}
// Обратный ход Гаусса с правой частью
void Matrix :: GaussRevMove(Matrix* rightMatrix){
    for(int i = getSize() - 1; i > 0; i--) {
        if (!matrix[i][i])
            continue;
        for (int j = i - 1; j >= 0; j--) {
            double multiplier = matrix[j][i];
            minLines(i, j, multiplier, i);
            rightMatrix->minLines(i, j, multiplier);
        }
    }
}
// Создание матрицы и вектора для метода простых итераций
pair<Matrix, Vector> Matrix :: alpha_betta_Create(Vector* b){
    Matrix Alpha(this->getSize());
    Vector Betta(this->getSize());
    for(int i  = 0; i < this->getSize(); i++) {
        for (int j = 0; j < this->getSize(); j++){
            if(i != j) {
                Alpha.changeVal(i, j, -(this->getVal(i, j) / this->getVal(i, i)));
                divNum++;
            }
            else
                Alpha.changeVal(i, j, 0);
        }
        divNum++;
        Betta.changeElem(b->getElem(i) / this->getVal(i, i), i);
    }
    return make_pair(Alpha, Betta);
}
// Поиск максимального элемента матрицы вне диагонали
pair <int, int> Matrix ::findMax() {
    double max = -INFINITY;
    pair <int, int> answer;
    for(int i = 0; i < getSize(); i++)
        for(int j = 0; j < getSize(); j++)
            if(i != j && abs(this->getVal(i, j)) > max){
                max = abs(this->getVal(i, j));
                answer = {i, j};
            }
    return answer;
}