#include <iostream>

// 参数：一个order阶矩阵，和矩阵的阶数
void lowerUpperFactor(double *matrix, int order) {
    printf("--------原矩阵：--------\n");
    printMatrix(matrix, order,order);
    
    // 结果变量  L矩阵和U矩阵都是order阶矩阵
    double *L = new double[order*order];
    double *U = new double[order*order];

    // 初始化全为0
    for (int i = 0; i < order; i++) {
        // 初始化U下三角为0
        for (int j = 0; j < i; j++) {
            *(U + i * order + j) = 0;
        }
        //初始化L对角线为1，上三角为0
        *(L + i * order + i) = 1;
        for (int j = i + 1; j < order; j++) {
            *(L + i * order + j) = 0;
        }
    }
    // 计算U的第一行和L的第一列
    int i = 0;
    for (i = 0; i < order; i++) {
        *(U + i) = *(matrix + i);
    }
    for (i = 1; i < order; i++) {
        *(L + i * order) = *(matrix + i * order) / *U;
    }
    // 计算其余行列
    int temp;
    for (int i = 1; i < order; i++) {
        // 计算矩阵U
        for (int j = i; j < order; j++) {
            temp = 0;
            for (int k = 0; k < i; k++) {
                temp+= (*(U + k * order + j) * (*(L + i * order + k)));
            }
            *(U + i * order + j) = *(matrix + i * order + j) - temp;
        }
        // 计算矩阵L
        for (int j = i+1; j < order; j++) {
            temp = 0;
            for (int k = 0; k < i; k++) {
                temp += *(U + k * order + i) * (*(L + j * order + k));
            }
            *(L + j * order + i) = (*(matrix +j * order + i) - temp) / (*(U+i* order + i));
        }
    }

    printf("------矩阵U------\n");
    printMatrix(U, order,order);
    printf("------矩阵L------\n");
    printMatrix(L, order, order);

    if (L) {
        delete[] L;
    }
    if (U) {
        delete[] U;
    }
}
void printMatrix(double *matrix, int row, int column) {
    for(int i = 0; i < row; i++) {
        for (int j = 0; j < column; j++) {
            printf("%6.2lf ", *(matrix + i * column + j));
        }
        printf("\n");
    }

}

int main() {

    double matrix[] = { 1,2,3,0,1,5,6,0 };
    int order = 3;
    // double matrix1[] = { 1,1,-1,2,1,2,0,2,-1,-1,2,0,0,0,-1,1 };
    // int order1 = 4;

    lowerUpperFactor(matrix, order);
}
