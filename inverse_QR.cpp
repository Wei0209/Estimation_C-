#include <iostream>
#include <vector>
#include <cmath>
#include "fMatrix.h"
using namespace std;

int main() /* 矩阵A的QR分解*/
{
    vector<vector<double>> a = { {1,2,3},{0,1,4},{5,6,0} };
    int n = a.size();
    vector<vector<double>> q(n, vector<double>(n));
    vector<vector<double>> r(n, vector<double>(n));

    cout << "A:" << endl; //输出矩阵A
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            printf("%.4f ", a[i][j]);
        }
        cout << endl;
    }

    for (int k = 0; k < n; k++)
    {
        double MOD = 0;
        for (int i = 0; i < n; i++)
        {
            MOD += a[i][k] * a[i][k]; 
        }
        r[k][k] = sqrt(MOD); // 计算A第k列的模长，由公式(4)等于R的对角线元素||A:k||
        for (int i = 0; i < n; i++)
        {
            q[i][k] = a[i][k] / r[k][k]; // 由公式(2)，A第k列标准化之后成为Q的第k列
        }

        for (int i = k + 1; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                r[k][i] += a[j][i] * q[j][k]; // 由公式(4)，计算R的上三角部分
            }
            for (int j = 0; j < n; j++)
            {
                a[j][i] -= r[k][i] * q[j][k]; // 由公式(1)，计算更新A的每一列
            }
        }
    }

    cout << endl;
    cout << "Q:" << endl; //输出矩阵Q
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            printf("%.4f ", q[i][j]);
        }
        cout << endl;
    }

    cout << endl;
    cout << "R:" << endl; //输出矩阵R
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            printf("%.4f ", r[i][j]);
        }
        cout << endl;
    }

    return 0;
}