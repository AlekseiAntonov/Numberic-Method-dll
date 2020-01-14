//
// Created by user on 12.11.2019.
//

#ifndef MATRIXDLL_MATRIX_WRAP_H
#define MATRIXDLL_MATRIX_WRAP_H

#ifdef __cplusplus
extern "C" {
#endif

struct mather;

typedef struct mather Matrix_t;

typedef struct{
    double* solution;
    double* err;
    int iterNum;
}answer_t;

typedef struct {
    answer_t* answer;
    double * operations;
}Kurs_answer_t;

typedef struct {
    double* eigenvalue;
    double* eigenvector;
    double* IterMax;
    int iterNum;
}APSZ_t;


Matrix_t * matrix_create(double* data, int size); // создание матрицы с данными

void matrix_destroy(Matrix_t *matrix); // Удаление матрицы

void matrix_print(Matrix_t* matrix); // Печать матрицы

Matrix_t* matrix_inverse(Matrix_t *matrix); // Обратная матрица

int matrix_getSize(Matrix_t* matrix); // Вернуть размер матрицы

double matrix_determinate(Matrix_t *matrix); // Найти определитель

void matrix_fill(Matrix_t *matrix, int size); // Заполнение матрицы из консоли

//bool matrix_isIdentity(Matrix_t *matrix);

void matrix_refill(Matrix_t *matrix , double** data, int size); // Заполенение матрицы из массива double

void matrix_randFill(Matrix_t *matrix, int size); // Заполнение матрицы случайными числами

void matrix_changeVal(Matrix_t *matrix, int line, int column, double newVal); // Изменение элемента матрицы

double matrix_getVal(Matrix_t *matrix, int line, int column); // Возврашает элемент матрицы

void matrix_setIdentity(Matrix_t *matrix); // Преобразование матрицы в единичную

void matrix_transpose(Matrix_t *matrix); // Транспонирование матрицы

double matrix_norm(Matrix_t *matrix); // Норма матрицы

double matrix_cond(Matrix_t *matrix); // Число обусловленности

double* matrix_FindSolution_LU(Matrix_t* matrix, double* b);  // Поиск решения LU разложением

double* matrix_FindSolution_MSI(Matrix_t* matrix, double* b, double epsilon); // Поиск решения методом простых итераций

answer_t* matrix_FindSolution_MSI_Convergence(Matrix_t* matrix, double* b, double* x, double epsilon); // Поиск решения методом простых итераций с данными

double* matrix_getAnswer(answer_t* answer); // Ответ метода простых итераций

double* matrix_getErrs(answer_t* answer); // Погрешность метода простых итераций

void matrix_GenerateDiagWithCond(Matrix_t* matrix, int degree); // Генерирование диагональной матрицы с числом обусловленности

double* matrix_mulVector(Matrix_t* matrix, double* vector); // Умножение матрицы на вектор

int matrix_getIterNum(answer_t* answer); // Число итераций

APSZ_t* matrix_JacobiRotateMethod(Matrix_t* matrix); // Поиск решения методом Якоби

double* matrix_GetEigenvalue(APSZ_t* source); // Собственные числа

double* matrix_GetEigenvector(APSZ_t* source); // Собственные вектора

double* matrix_GetIterMax(APSZ_t* source); // Максимальный элемент на каждой итерации

int matrix_iterNum(APSZ_t* source); // Число итераций

double* matrix_return(Matrix_t* matrix); // Преобразование матрицы в массив

void matrix_GenerateWithSmallDet(Matrix_t* matrix);

APSZ_t* matrix_JacobiRotateMethod_WithConvergence(Matrix_t* matrix, double* Solution);

Kurs_answer_t*  matrix_FindSolution_MSI_Convergence_With_Op(Matrix_t* matrix, double* b, double* x, double epsilon, int lim);

double* matrix_getOperations(Kurs_answer_t* answer); // Число операций

answer_t* matrix_getAnswerStruct(Kurs_answer_t* answer);

Kurs_answer_t* matrix_FindSolution_LU_With_Op(Matrix_t* matrix, double* b);
#ifdef __cplusplus
}
#endif

#endif //MATRIXDLL_MATRIX_WRAP_H
