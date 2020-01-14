//
// Created by user on 29.10.2019.
//
#include "Matrix.h"
#include <algorithm>

void Matrix :: hilbert(int size){ // Генерация матрицы Гильберта
    matrix.resize(size);
    for(int i = 0; i < size; i++)
        for(int j = 0; j < size; j++)
            matrix[i].push_back(1.0/(double)(i + j + 1));

}

Matrix Matrix :: inverse(){ // Обратная матрица
    if(this->determinate() == 0)
        return {};
    Matrix answer(getSize());
    Matrix copy = *this;
    answer.setIdentity();
    if(copy.GaussMove(&answer) == 0)
        return {};
    copy.GaussRevMove(&answer);
    return answer;
}

double Matrix :: determinate(){ // Определитель
    Matrix copy = *this;
    return copy.GaussMove();
}

void Matrix :: fill(int size){ // Заполнение матрицы из консоли
    matrix.resize(size, {});
    for(int i = 0; i < size; i++)
        for (int j = 0; j < size; j++) {
            matrix[i].push_back(0);
            cin >> matrix[i][j];
        }
}

bool Matrix :: isIdentity(){ // Проверка, является ли матрица единичной
    for(int i = 0; i < getSize(); i++)
        for(int j = 0; j < getSize(); j++) {
            if (i == j)
                if (matrix[i][j] != 1)
                    return false;
            if (i != j)
                if(matrix[i][j])
                    return false;
        }
    return true;
}
bool Matrix :: isDiag(){ // Проверка, является ли матрица диагональной
    for(int i = 0; i < getSize(); i++)
        for(int j = 0; j < getSize(); j++) {
            if (i != j)
                if(matrix[i][j])
                    return false;
            }
    return true;
}

void Matrix :: print() const { // Печать матрицы
    for(const auto& str : matrix) {
        for (auto elem : str)
            cout << elem << " ";
        cout << endl;
    }
}

void Matrix :: reFill(double** data, int size){ // Заполенение матрицы из массива double
    matrix.resize(size);
    for(int i = 0; i < size; i++)
        for(int j = 0; j < size; j++)
            matrix[i].push_back(data[i][j]);
}

void Matrix :: randFill(int size){ // Заполнение матрицы случайными числами
    matrix.resize(size, {});
    for(int i = 0; i < size; i++)
        for (int j = 0; j < size; j++) {
            matrix[i].push_back(RND());
        }
}

int Matrix :: getSize() const { // Возврашает размер матрицы
    return matrix.size();
}

void Matrix :: changeVal(int line, int column, double newVal){ // Изменение элемента матрицы
    if(line >= matrix.size() || column >= matrix.size())
        return;
    matrix[line][column] = newVal;
}

double Matrix :: getVal(int line, int column) const { // Возвращает элемент матрицы
    return matrix[line][column];
}

void Matrix :: setIdentity(){ // Преобразование матрицы в единичную
    for(unsigned int i = 0; i < matrix.size(); i++)
        for(unsigned int j = 0; j < matrix.size(); j++)
            matrix[i][j] = (i == j);
}

void Matrix :: transpose(){ // Транспонирование
    for(unsigned int i = 0; i < matrix.size(); i++)
        for(unsigned int j = i; j < matrix.size(); j++)
            swap(matrix[i][j], matrix[j][i]);
}

double Matrix :: norm(){ // Норма матрицы
    double max = -INFINITY;
    for(auto & j : matrix) {
        double sum = 0;
        for (auto & ai : j) {
            sum += abs(ai);
            addNum++;
        }
        if(sum > max)
            max = sum;

    }
    return max;

}

double Matrix :: cond(){ // Число обусловленности
    return this->norm() * (this->inverse()).norm();
}

Vector & Matrix :: operator *= (Vector const& vec){ // *= оператор
    auto* answer = new(Vector);
    *answer = this->matrixMulVector(vec);
    return *answer;
}

Vector Matrix :: operator*(Vector const&  vec){
    Vector answer = this->matrixMulVector(vec);
    return answer;
}

Matrix& Matrix :: operator *= (Matrix A) {
    *this = mul(&A);
    return *this;
}

Matrix Matrix :: operator*(Matrix  const& B){
    return *this *= B;
}

bool Matrix :: operator == (Matrix const &B){
    if(this->getSize() != B.getSize())
        return false;
    for(int i = 0; i < this->getSize(); i++)
        for(int j = 0; j < this->getSize(); j++)
            if(abs(this->getVal(i, j) - B.getVal(i, j)) > abs(0.0001)
               && abs(this->getVal(i, j) - B.getVal(i, j)) > abs(0.0001))
                return false;
    return true;
}

Vector Matrix :: FindSolution_LU(double* b_data){ // Поиск решения LU разложением
    Vector b(b_data, this->getSize());
    pair<Matrix, Matrix> LU = this->LU_Decomposition(b);
    Vector y = LU.first.inverse() * b;
    Vector x = LU.second.inverse() * y;
    addNum += LU.first.addNum + LU.second.addNum;
    subNum += LU.first.subNum + LU.second.subNum;
    divNum += LU.first.divNum + LU.second.divNum;
    multNum += LU.first.multNum + LU.second.multNum;
    multNum += 2 * getSize() * getSize();
    return x;
}

Vector Matrix :: FindSolution_MSI(double * b_data, double epsilon){ // Поиск решения методом простых итераций
    Vector b(b_data, getSize());
    pair <Matrix, Vector> Alpha_Betta = this->alpha_betta_Create(&b);
    Matrix Alpha = Alpha_Betta.first;
    Vector Betta = Alpha_Betta.second;
    double q = Alpha.norm();;
    Vector x_cur;
    Vector x_next = Betta;
    do{
        x_cur = x_next;
        x_next = Betta + Alpha * x_cur;
    }while ((x_cur - x_next).norm()>  epsilon);
    return x_cur;
}

pair<Vector, vector<double>> Matrix :: FindSolution_MSI(double* b_data, double* x, double epsilon){ // Поиск решения методом простых итераций с данными
    Vector b(b_data, getSize());
    pair <Matrix, Vector> Alpha_Betta = this->alpha_betta_Create(&b);
    Matrix Alpha = Alpha_Betta.first;
    Vector Betta = Alpha_Betta.second;
    double q = Alpha.norm();
    Vector x_cur;
    Vector x_next = Betta;
    vector <double> errs;
    Vector answer(x, getSize());
    cout << Alpha.cond();
    do{
        x_cur = x_next;
        errs.push_back((answer - x_cur).norm());

        x_next = Betta + Alpha * x_cur;
    }while ((x_cur - x_next).norm() > abs(1 - q) /q * epsilon);
    return make_pair(x_cur, errs);
}

void Matrix :: GenerateDiagWithCond(int condDegree){ // Генерирование матрицы с числом обусловленности
    vector<double> diag;
    double maxElem = -INFINITY;
    for(int i = 0; i < getSize(); i++){
        diag.push_back(RND());
        maxElem = max(maxElem, diag[i]);
    }
    diag[getSize() - 1] = maxElem / pow(10,condDegree);
    for(int i = 0; i < getSize(); i++) {
        changeVal(i, i, diag[i]);
        for (int j = 0; j < getSize(); ++j) {
            if (i != j)
                this->changeVal(i, j, RND(-abs(getVal(i, i)) / getSize(), abs(getVal(i, i)) / getSize()));
        }
    }
}

pair<JacobiAns_t,Matrix> Matrix :: JacobiRotateMethod(){ // Поиск с.ч. и с.в. методом Якоби
    Matrix Lambda(getSize());
    Lambda.setIdentity();
    Matrix Ak = *this;
    auto LineCol = findMax();
    double epsilon = 1e-10;
    double maxElem = getVal(LineCol.first, LineCol.second);
    vector <double> iterMax;
    int iter = 0;
    while (abs(maxElem) > epsilon){
        iter++;
        iterMax.push_back(abs(maxElem));
        Matrix B = Ak;
        int i = LineCol.first;
        int j = LineCol.second;
        //cout << '[' << i << ' ' << j << ']' << endl;
        double phi = 0.5 * atan(2.0 * maxElem / (Ak.getVal(i, i) - Ak.getVal(j, j)));
        /*double c = cos(phi);
        double s = sin(phi);
        for(int m = 0; m < getSize(); m++) {
            Lambda.changeVal(m, i, c * Lambda.getVal(m, i) - s * Lambda.getVal(m, j));
            Lambda.changeVal(m, j, s * Lambda.getVal(m, i) + c * Lambda.getVal(m, j));

            B.changeVal(i, m, c * Ak.getVal(m, i) + s * Ak.getVal(m, j));
            B.changeVal(m, i, Ak.getVal(i, m));
            B.changeVal(j, m, -s * Ak.getVal(m, i) + c * Ak.getVal(m, j));
            B.changeVal(m, j, Ak.getVal(j, m));
        }
        B.changeVal(i, i, c * c * Ak.getVal(i, i) + s * s * Ak.getVal(j, j) + 2 * s * c * Ak.getVal(i, j));
        B.changeVal(j, j, s * s * Ak.getVal(i, i) + c * c * Ak.getVal(j, j) - 2 * s * c * Ak.getVal(i, j));

        B.changeVal(i, j, 0);
        B.changeVal(j, i, 0);
        Ak = B;*/
        Ak.rotate(phi, i, j);
        Ak.changeVal(i,j, 0);
        Ak.changeVal(j,i,0);
        LineCol = Ak.findMax();
        maxElem = Ak.getVal(LineCol.first, LineCol.second);
    }
    vector <double> result;
    result.reserve(getSize());
    for(int i = 0; i < getSize(); i++)
        result.push_back(Ak.getVal(i,i));
    sort(result.begin(), result.end());
    Vector res(result);
    Vector iterations(iterMax);
    return {{res, iterations, iter}, Lambda};
}

double* Matrix :: toDouble(){ // Преобразование матрицы в массив
    auto* answer = new double[getSize() * getSize()];
    for(int i = 0; i < getSize(); i++)
        for(int j = 0; j < getSize(); j++)
            answer[i * getSize() + j] = getVal(j, i);
    return answer;
}

void Matrix :: GenerateWithSmallDet(){
    vector<double> diag;
    for(int i = 0; i < getSize(); i++) {
        diag.push_back(RND());
    }
    for(int i = 0; i < getSize() - 1; i++) {
        changeVal(i, i, diag[i]);
        for (int j = 0; j < getSize(); ++j) {
            if (i != j)
                this->changeVal(i, j, RND(-abs(getVal(i, i)) / getSize(), abs(getVal(i, i)) / getSize()));
        }
    }
    for(int i = 0; i < getSize() - 1;i++) this->changeVal(getSize() -1, i, 0);
    this->changeVal(getSize() - 1, getSize() - 1, 1e-50);
}

void Matrix ::rotate(double phi, int i, int j) {
    Matrix Jacobi(getSize());
    Jacobi.setIdentity();
    double c = cos(phi);
    double s = sin(phi);
    Jacobi.changeVal(i, i, c);
    Jacobi.changeVal(j, j, c);
    Jacobi.changeVal(i, j, -s);
    Jacobi.changeVal(j, i, s);
    *this *= Jacobi;
    Jacobi.transpose();
    Matrix mul = Jacobi * (*this);
    *this = mul;

}

pair<JacobiAns_t,Matrix> Matrix :: JacobiRotateMethod(double* Solution){ // Поиск с.ч. и с.в. методом Якоби
    Matrix Lambda(getSize());
    Lambda.setIdentity();
    Matrix Ak = *this;
    auto LineCol = findMax();
    double epsilon = 1e-15;
    double maxElem = getVal(LineCol.first, LineCol.second);
    vector <double> iterMax;
    int iter = 0;
    Vector answer(Solution, getSize());
    while (abs(maxElem) > epsilon){
        iter++;
        vector <double> kSol;
        for(int k = 0; k < getSize(); k++)
            kSol.push_back(Ak.getVal(k , k));
        sort(kSol.begin(),kSol.end());
        Vector kSolution(kSol);
        iterMax.push_back((kSolution - answer).norm());
        //iterMax.push_back(abs(maxElem));
        Matrix B = Ak;
        int i = LineCol.first;
        int j = LineCol.second;
        //cout << '[' << i << ' ' << j << ']' << endl;
        double phi = 0.5 * atan(2.0 * maxElem / (Ak.getVal(i, i) - Ak.getVal(j, j)));
        Ak.rotate(phi, i, j);
        Ak.changeVal(i,j, 0);
        Ak.changeVal(j,i,0);
        LineCol = Ak.findMax();
        maxElem = Ak.getVal(LineCol.first, LineCol.second);
    }
    vector <double> result;
    result.reserve(getSize());
    for(int i = 0; i < getSize(); i++)
        result.push_back(Ak.getVal(i,i));
    sort(result.begin(), result.end());
    Vector res(result);
    Vector iterations(iterMax);
    return {{res, iterations, iter}, Lambda};
}

pair<Vector, vector<double>> Matrix :: FindSolution_MSI_withLimits(double* b_data, double* x, double epsilon, int operationsNum){ // Поиск решения методом простых итераций с данными
    Vector b(b_data, getSize());
    pair <Matrix, Vector> Alpha_Betta = this->alpha_betta_Create(&b);
    Matrix Alpha = Alpha_Betta.first;
    Vector Betta = Alpha_Betta.second;
    double q = Alpha.norm();
    Vector x_cur;
    Vector x_next = Betta;
    vector <double> errs;
    Vector answer(x, getSize());
    int curOperNum  = 0;
    do{
        x_cur = x_next;
        errs.push_back((answer - x_cur).norm());
        x_next = Betta + Alpha * x_cur;
        multNum += getSize() * getSize();
        addNum += getSize()  + getSize() * getSize();
        subNum += getSize();
        curOperNum = (addNum + subNum + divNum + multNum);
    }while ((curOperNum < operationsNum || operationsNum == 0) && (x_cur - x_next).norm() > abs(1 - q) /q * epsilon);
    return make_pair(x_cur, errs);
}