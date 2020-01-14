//
// Created by user on 12.11.2019.
//

//
// Created by user on 30.10.2019.
//

#include <cstdlib>
#include "Matrix_Wrap.h"
#include "Matrix.h"

struct mather {
    void *obj;
};

Matrix_t * matrix_create(double* data, int size){
    Matrix_t *m;
    Matrix *obj;
    m = (Matrix_t *)malloc(sizeof(*m));
    obj = new Matrix(data, size);
    m->obj = obj;
    return m;
}

void matrix_destroy(Matrix_t *m){
    if (m == nullptr)
        return;
    delete static_cast<Matrix *>(m->obj);
    free(m);
}

void matrix_print(Matrix_t *matrix){
    Matrix *obj;
    if (matrix == nullptr)
        return;
    obj = static_cast<Matrix *>(matrix->obj);
    obj->print();
}

Matrix_t* matrix_inverse(Matrix_t *matrix){
    Matrix *obj;
    if (matrix == nullptr)
        return nullptr;
    obj = static_cast<Matrix *>(matrix->obj);
    Matrix_t *m;
    m = (Matrix_t *)malloc(sizeof(*m));
    obj->inverse();
    m->obj = obj;
    return m;
}

double matrix_determinate(Matrix_t *matrix){
    Matrix *obj;
    if (matrix == nullptr)
        return 0;
    obj = static_cast<Matrix *>(matrix->obj);
    return  obj->determinate();
}

void matrix_fill(Matrix_t *matrix, int size){
    Matrix *obj;
    if (matrix == nullptr)
        return;
    obj = static_cast<Matrix *>(matrix->obj);
    obj->fill(size);
}

bool matrix_isIdentity(Matrix_t *matrix){
    Matrix *obj;
    if (matrix == nullptr)
        return false;
    obj = static_cast<Matrix *>(matrix->obj);
    return obj->isIdentity();
}

void matrix_refill(Matrix_t *matrix , double** data, int size){
    Matrix *obj;
    if (matrix == nullptr)
        return;
    obj = static_cast<Matrix *>(matrix->obj);
    obj->reFill(data, size);
}

void matrix_randFill(Matrix_t *matrix, int size){
    Matrix *obj;
    if (matrix == nullptr)
        return;
    obj = static_cast<Matrix *>(matrix->obj);
    obj->randFill(size);
}

int matrix_getSize(Matrix_t *matrix){
    Matrix *obj;
    if (matrix == nullptr)
        return 0;
    obj = static_cast<Matrix *>(matrix->obj);
    return obj->getSize();
}

void matrix_changeVal(Matrix_t *matrix, int line, int column, double newVal){
    Matrix *obj;
    if (matrix == nullptr)
        return;
    obj = static_cast<Matrix *>(matrix->obj);
    obj->changeVal(line, column, newVal);
}

double matrix_getVal(Matrix_t *matrix, int line, int column){
    Matrix *obj;
    if (matrix == nullptr)
        return 0;
    obj = static_cast<Matrix *>(matrix->obj);
    return obj->getVal(line, column);
}

void matrix_setIdentity(Matrix_t *matrix){
    Matrix *obj;
    if (matrix == nullptr)
        return;
    obj = static_cast<Matrix *>(matrix->obj);
    obj->setIdentity();
}

void matrix_transpose(Matrix_t *matrix){
    Matrix *obj;
    if (matrix == nullptr)
        return;
    obj = static_cast<Matrix *>(matrix->obj);
    obj->transpose();
}

double matrix_norm(Matrix_t *matrix){
    Matrix *obj;
    if (matrix == nullptr)
        return 0 ;
    obj = static_cast<Matrix *>(matrix->obj);
    return obj->norm();
}

double matrix_cond(Matrix_t *matrix){
    Matrix *obj;
    if (matrix == nullptr)
        return 0;
    obj = static_cast<Matrix *>(matrix->obj);
    return obj->cond();
}

double* matrix_FindSolution_LU(Matrix_t* matrix, double* b){
    Matrix *obj;
    if (matrix == nullptr)
        return nullptr;
    obj = static_cast<Matrix *>(matrix->obj);
    return obj->FindSolution_LU(b).toDouble();
}
double* matrix_FindSolution_MSI(Matrix_t* matrix, double* b, double epsilon){
    Matrix *obj;
    if (matrix == nullptr)
        return nullptr;
    obj = static_cast<Matrix *>(matrix->obj);
    return obj->FindSolution_MSI(b, epsilon).toDouble();
}
void matrix_GenerateDiagWithCond(Matrix_t* matrix, int degree){
    Matrix *obj;
    if (matrix == nullptr)
        return;
    obj = static_cast<Matrix *>(matrix->obj);
    obj->GenerateDiagWithCond(degree);
}

double* matrix_mulVector(Matrix_t* matrix, double* vec){
    Matrix *obj;
    obj = static_cast<Matrix *>(matrix->obj);
    Vector b(vec, obj->getSize());
    return (*obj * b).toDouble();
}

answer_t* matrix_FindSolution_MSI_Convergence(Matrix_t* matrix, double* b, double* x, double epsilon){
    Matrix *obj;
    if (matrix == nullptr)
        return nullptr;
    obj = static_cast<Matrix *>(matrix->obj);
    pair<Vector, vector<double>> answer = obj->FindSolution_MSI(b, x, epsilon);
    auto* answer1 = new answer_t;
    Vector err(answer.second);
    answer1->err = err.toDouble();
    answer1->iterNum = err.getSize();
    answer1->solution = answer.first.toDouble();
    return answer1;
}


double* matrix_getAnswer(answer_t* answer){
    return answer->solution;
}

double* matrix_getErrs(answer_t* answer){
    return answer->err;
}

int matrix_getIterNum(answer_t* answer){
    return answer->iterNum;
}

APSZ_t* matrix_JacobiRotateMethod(Matrix_t* matrix){
    Matrix *obj;
    if (matrix == nullptr)
        return nullptr;
    obj = static_cast<Matrix *>(matrix->obj);
    auto* answer1 = new APSZ_t;
    auto  answer = obj->JacobiRotateMethod();
    answer1->eigenvalue = answer.first.res.toDouble();
    answer1->eigenvector = answer.second.toDouble();
    answer1->IterMax = answer.first.iterMax.toDouble();
    answer1->iterNum =answer.first.iterNum;
    return answer1;
}

double* matrix_GetEigenvalue(APSZ_t* source){
    return source->eigenvalue;
}

double* matrix_GetEigenvector(APSZ_t* source){
    return source->eigenvector;
}

double* matrix_GetIterMax(APSZ_t* source){
    return source->IterMax;
}

int matrix_iterNum(APSZ_t* source){
    return source->iterNum;
}

double* matrix_return(Matrix_t* matrix){
    Matrix *obj;
    if (matrix == nullptr)
        return nullptr;
    obj = static_cast<Matrix *>(matrix->obj);
    return obj->toDouble();
}

void matrix_GenerateWithSmallDet(Matrix_t* matrix){
    Matrix *obj;
    if (matrix == nullptr)
        return;
    obj = static_cast<Matrix *>(matrix->obj);
    obj->GenerateWithSmallDet();
}

APSZ_t* matrix_JacobiRotateMethod_WithConvergence(Matrix_t* matrix, double* Solution){
    Matrix *obj;
    if (matrix == nullptr)
        return nullptr;
    obj = static_cast<Matrix *>(matrix->obj);
    auto* answer1 = new APSZ_t;
    auto  answer = obj->JacobiRotateMethod(Solution);
    answer1->eigenvalue = answer.first.res.toDouble();
    answer1->eigenvector = answer.second.toDouble();
    answer1->IterMax = answer.first.iterMax.toDouble();
    answer1->iterNum =answer.first.iterNum;
    return answer1;
}

Kurs_answer_t*  matrix_FindSolution_MSI_Convergence_With_Op(Matrix_t* matrix, double* b, double* x, double epsilon, int lim){
    Matrix *obj;
    if (matrix == nullptr)
        return nullptr;
    obj = static_cast<Matrix *>(matrix->obj);
    pair<Vector, vector<double>> answer = obj->FindSolution_MSI_withLimits(b, x, epsilon, lim);
    auto* answer1 = new answer_t;
    Vector err(answer.second);
    answer1->err = err.toDouble();
    answer1->iterNum = err.getSize();
    answer1->solution = answer.first.toDouble();

    auto* answer2 = new Kurs_answer_t;
    answer2->answer = answer1;
    auto operations = new double[5];
    operations[0] = (double)obj->addNum;
    operations[1] = (double)obj->subNum;
    operations[2] = (double)obj->multNum;
    operations[3] = (double)obj->divNum;
    operations[4] = operations[0] + operations[1] + operations[2] + operations[3];

    answer2->operations = operations;

    obj->addNum = 0;
    obj->subNum = 0;
    obj->multNum = 0;
    obj->divNum = 0;

    return answer2;
}

double* matrix_getOperations(Kurs_answer_t* answer){
    return answer->operations;
}

answer_t* matrix_getAnswerStruct(Kurs_answer_t* answer){
    return answer->answer;
}

Kurs_answer_t* matrix_FindSolution_LU_With_Op(Matrix_t* matrix, double* b){
    Matrix *obj;
    if (matrix == nullptr)
        return nullptr;
    obj = static_cast<Matrix *>(matrix->obj);
    auto* answer1 = new answer_t;
    answer1->solution = obj->FindSolution_LU(b).toDouble();

    auto* answer2 = new Kurs_answer_t;
    answer2->answer = answer1;
    auto operations = new double[5];
    operations[0] = (double)obj->addNum;
    operations[1] = (double)obj->subNum;
    operations[2] = (double)obj->multNum;
    operations[3] = (double)obj->divNum;
    operations[4] = operations[0] + operations[1] + operations[2] + operations[3];

    obj->addNum = 0;
    obj->subNum = 0;
    obj->multNum = 0;
    obj->divNum = 0;


    answer2->operations = operations;

    return answer2;
}