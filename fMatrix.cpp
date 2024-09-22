#include <iostream>
#include "fMatrix.h"
#include "fVector.h"
#include <cstring> 
#include <math.h>
#include <stdlib.h>
#include <cmath>
using namespace std;

fMatrix::fMatrix( int n_rows, int n_cols)
{   
    this->rows = n_rows;
    this->cols = n_cols;
    elem = new Float[rows*cols];
}
fMatrix::fMatrix( const fMatrix & a)
{
    rows = a.rows;
    cols = a.cols;
    elem = new Float[rows*cols];
    memcpy(elem,&a.elem[0],sizeof(Float)*rows*cols);
} 

fMatrix::fMatrix( Float *Array, int n_rows , int n_cols  )
{
    elem = NULL;
    rows = n_rows;
    cols = n_cols;
    elem = new Float[rows*cols];
    memcpy(elem,Array,sizeof(Float)*rows*cols);
}

fMatrix::fMatrix( int n_rows , int n_cols , Float *Array )
{
    elem = NULL;
    rows = n_rows;
    cols = n_cols;
    elem = new Float[rows*cols];
    memcpy(elem,Array,sizeof(Float)*rows*cols);
}

fMatrix::~fMatrix(){
    if (elem) {
        delete[] elem; // 釋放動態分配的矩陣資源
    }
    elem = nullptr; // 將elem指向空指針，避免出現野指針
    rows = 0; // 將rows設為0，表示矩陣已被刪除
    cols = 0; // 將cols設為0，表示矩陣已被刪除
}

fMatrix  &fMatrix::operator=( const fMatrix &M )
{
    if(this != &M)
    {
        if(this -> cols != M.cols && this -> rows != M.rows)
        {
            if(this -> elem)
            {
                delete[] this -> elem;
            }
            if(M.rows > 0 && M.cols > 0)
            {
                this -> rows =M.rows;
                this -> cols =M.cols;
                this -> elem = new Float[this -> rows * this -> cols ];
                memcpy(this->elem, M.elem, this -> rows *  this -> cols * sizeof(Float));
            }
            else
            {
                this->rows = 0;
                this->cols = 0; 
                this->elem = NULL;
            }
        }
        else
        {
            memcpy(this->elem, M.elem, this -> rows *  this -> cols * sizeof(Float));
        }
    }
    return *this;
}
fMatrix  &fMatrix::operator=( Float s )
{
    if(this->rows != 0 && this->cols != 0)
    {
        int size = this->rows * this->cols;
        for(int i = 0; i < size; ++i)
        {
            this -> elem[i] = s;

        }
    }
    return *this;
}
fMatrix &fMatrix::SwapRows(int i1, int i2)
{
    if (this->rows > i1 && this -> rows > i2 && this->cols != 0)
    {
        Float *row1 = this->elem + this->cols * i1;
        Float *row2 = this->elem + this->cols * i2;
        Float *temp = new Float[this->cols];
        int Crows = sizeof(Float)* this->cols;
        memcpy(temp, row1, Crows);
        memcpy(row1, row2, Crows);
        memcpy(row2, temp, Crows);
        delete[] temp;
    }
    return *this;
}
fMatrix &fMatrix::SwapCols(int j1, int j2)
{
    if (this->cols > j1 && this->cols > j2 && this->rows != 0)
    {
    //元素交換
        Float *col1 = this->elem + j1;
        Float *col2 = this->elem + j2;
        for (int i = 0; i < this->rows; ++i)
        {
            Float temp = *col1;
            *col1 = *col2;
            *col2 = temp;
            col1 += this->cols;
            col2 += this->cols;
        }
    }
    return *this;
}
void fMatrix::SetCol(int col, const fVector &ColVec)
{
    if (this->cols > col && this->rows == ColVec.Size())
    {
        Float *coln = this->elem + col;
        for (int i = 0; i < this->rows; ++i)
        {
            *coln = ColVec(i);
            coln += this->cols;
        }
    }
}
void fMatrix::SetRow(int row, const fVector &RowVec)

{
    if (this->rows > row && this->cols == RowVec.Size())
    {
        Float *rown = this->elem + this->cols * row;
        for (int i = 0; i < this->rows; ++i)
            {
                *(rown + i) = RowVec(i);
            }
    }
}
void    fMatrix::SetBlock( int imin, int imax, int jmin, int jmax, const fMatrix &Mat )
{
    if(imax >= imin && jmax >= jmin)
    {
        if(this -> rows > imax && this -> cols > jmax)
        {
            int copysize = (jmax - jmin + 1) * sizeof(Float);
            for (int i = imin;i <= imax; ++i)
            {
                Float *stemp = Mat.elem + Mat.cols * (i-imin);
                Float *dtemp = this -> elem + this -> cols * i + jmin;
                memcpy(dtemp, stemp, copysize);
            }
        }
    }
}
void    fMatrix::SetSize( int rows, int cols )
{
    if( rows != 0 && cols != 0)
    {
        if(this->rows == rows && this-> cols == cols)
        {
            memset(this->elem, 0, this->rows * this->cols * sizeof(Float));
        }
        
    }
    else
    {
        this->rows = rows;
        this->cols = cols;
        if (this->elem)
        {
            delete[] this->elem;
        }
        this->elem = new Float[this->rows * this->cols];
    }
}
fVector fMatrix::GetCol(int col) const
{
    fVector Vec(this->rows);
    if (this -> cols > col)
    {
        Float *temp = this->elem + col;
        for (int i = 0; i < this->rows; ++i)
        {
            Vec(i) = *(temp + this->cols * i);
        }
    }
    return Vec;
}
fVector fMatrix::GetRow(int row) const
{
    fVector Vec(this->cols);
    if (this -> rows > row)
    {
        Float *temp = this -> elem + this->cols * row;
        for (int i = 0; i < this->cols; ++i)
        {
            Vec(i) = *(temp + i);
        }
    }
    return Vec;
}
fMatrix  fMatrix::GetBlock( int imin, int imax, int jmin, int jmax ) const
{
    fMatrix Mat(imax - imin + 1,jmax -jmin +1);
    if(imax >= imin && jmax >= jmin)
    {
        if(this -> rows > imax || this -> cols > jmax )
        {
            for (int i = 0; i < Mat.rows; i++) 
            {
               Float *temp = this->elem + this->cols * (i + imin) + jmin;
               for (int j = 0; j < Mat.cols; ++j)
               {
                    Mat.elem[i,j] = *(temp + j);
               }
            }
           
        }
    }   
    return Mat;
}

fMatrix  Transp      ( const fMatrix &A )
{   
    fMatrix c(A.cols,A.rows);
    for(int i=0;i<A.cols;i++)
    {
        for(int j=0;j<A.rows;j++)
        {
            c.elem[i*A.rows+j] = A.elem[j*A.cols+i];
        }
    }
    return c;
}
fMatrix  AATransp    ( const fMatrix &a )
{
    return a * Transp(a);
}
fMatrix  ATranspA    ( const fMatrix &a )
{
    return Transp(a) * a;
}

fMatrix  Outer       ( const fVector &A, const fVector &B ) 
{
	fMatrix c(A.Size(),B.Size());
    for(int i=0;i<A.Size()*B.Size();i++)
    {
        c.elem[i] = A.Array()[i/c.rows]*B.Array()[i%c.rows];
    }
    return c;
}

fMatrix  Identity	 ( int nSize ) 
{
	fMatrix a(nSize, nSize);
	for (int i = 0; i < nSize; i++)
	{
		for (int j = 0; j < nSize; j++)
		{
			if(i == j)
			{
				a.elem[i*a.cols+j] = 1;
			}
			else
			{
				a.elem[i*a.cols+j] = 0;
			}
		}
	}
	return a;
}


fMatrix  Diag        ( const fVector &A)
{
    fMatrix c(A.Size(),A.Size());
    for(int i=0;i<A.Size()*A.Size();i++)
    {
        if(i/A.Size()==i%A.Size())
        {
            c.elem[i] = A.Array()[i/A.Size()];
        }
        else
        {
            c.elem[i] = 0;
        }    
    }
    return c;
}

fVector  Diag        ( const fMatrix &A )
{
    if(A.rows != A.cols)
    {
        cout << "shape error!" <<endl;
        return 0;
    }

    fVector c(A.rows);
    for(int i=0;i<A.rows*A.cols;i++)
    {
        if(i/A.rows==i%A.cols)
        {
            c.Array()[i/A.rows] = A.elem[i];
        }    
    }
    return c;
}


fMatrix  Diag        ( Float a, Float b, Float c ) 
{
	fMatrix fm(3, 3);
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			if(i == j)
			{
				if(i == 0)
				{
					fm.elem[i*fm.cols+j] = a;
				}
				else if(i == 1)
				{
					fm.elem[i*fm.cols+j] = b;
				}
				else if(i == 2)
				{
					fm.elem[i*fm.cols+j] = c;
				}			
			}
			else
			{
				fm.elem[i*fm.cols+j] = 0;
			}
		}
	}
	return fm;
}
double   Determinant ( const fMatrix &A )
{
    if(A.cols!=A.rows)
    {
        cout << "shape error!" <<endl;
        return 0;
    }
    double mul;
    double sum = 0;

    if (A.rows == 2)
    {
        sum = A.elem[0]*A.elem[3]-A.elem[1]*A.elem[2];
        return sum;
    }
    
    for(int i=0;i<A.cols;i++)
    {
        int n=0;
        fMatrix c(A.cols-1,A.cols-1);
        for(int j=A.cols;j<A.cols*A.rows;j++)
        {
            if(j%A.cols != i)
            {
                c.elem[n] = A.elem[j];
                n++;
                // cout << j << " ";
            }
        }
        // cout << endl;
        
        if(i%2 == 0)
            sum += A.elem[i]*Determinant(c);
            // cout << Determinant(c) << " ";
        if(i%2 == 1)
            sum -= A.elem[i]*Determinant(c);
            // cout << Determinant(c) << " ";
    }
    return sum;
}
double   Trace       ( const fMatrix &A )
{
    if(A.cols!=A.rows)
    {
        cout << "shape error!" <<endl;
        return 0;
    }
    double sum = 0;
    for(int i=0;i<A.rows*A.cols;i++)
    {
        if(i/A.cols==i%A.rows)
        {
            sum += A.elem[i];
        }
    }
    return sum;
}

double   OneNorm     ( const fMatrix &a )
{
    double sum = 0;
    for (int i = 0; i < a.rows; ++i) {
        double row_sum = 0.0;
        for (int j = 0; j < a.cols; ++j) {
            row_sum += abs(a.elem[i,j]);
        }
        sum = max(sum, row_sum);
    }
    return sum;
}
double   InfNorm     ( const fMatrix &a )
{
    double max_row_sum = 0.0;
    for (int i = 0; i < a.rows; ++i) {
        double row_sum = 0.0;
        for (int j = 0; j < a.cols; ++j) {
            row_sum += abs(a.elem[i,j]);
        }
        max_row_sum = max(max_row_sum, row_sum);
    }
    return max_row_sum;
}
fMatrix  Inverse  ( const fMatrix &A)
{
    if(A.cols!=A.rows)
    {
        cout << "shape error!" <<endl;
        return 0;
    }
    if(A.cols==2)
    {
        fMatrix c(2,2);
        c.elem[0] = A.elem[3];
        c.elem[1] = -A.elem[1];
        c.elem[2] = -A.elem[2];
        c.elem[3] = A.elem[0];
        return c/Determinant(A);
    }
    fMatrix c(A.rows,A.cols);
    fMatrix d(A.rows-1,A.cols-1);
    for(int i=0;i<A.cols*A.rows;i++)
    {
        // fMatrix d(A.rows-1,A.cols-1);
        int n = 0;
        for(int j=0;j<A.cols*A.rows;j++)
        {
            if(j/A.rows != i/A.rows && j%A.cols != i%A.cols)
            {
                d.elem[n] = A.elem[j];
                // cout<<j;
                n++;
            }
        }
        if (i%2 == 0)
        {
            c.elem[i] = Determinant(d)/Determinant(A);
        }
        else
        {
            c.elem[i] = -Determinant(d)/Determinant(A);
        }
    }
    return Transp(c);
}
fMatrix  Cholesky	 ( const fMatrix &fm1 ) // Computes Cholesky decomposition of a square matrix.	
{ 
	fMatrix fm(fm1.rows, fm1.cols);
	fm = 0;
	for (int i = 0; i < fm.cols; i++)
	{
		for (int j = 0; j < (i+1); j++) 
		{
			double val1 = 0;
			for (int k = 0; k < j; k++)
			{
				val1 += fm.elem[i*fm.cols+k] * fm.elem[j*fm.cols+k];
			}
			if (i == j)
			{
				fm.elem[i*fm.cols+j] = sqrt(fm1.elem[i*fm1.cols+i] - val1);
			}
			else
			{
				fm.elem[i*fm.cols+j] = (1.0 / fm.elem[j*fm.cols+j] * (fm1.elem[i*fm.cols+j] - val1));
			}			
		}
	}
	return fm;
}
void fMatrix::fill(double data)
{   
    cout << this->elem[2,2] << "\n";
    for(int i=0;i<rows;i++){
        for (int j=0;j<cols;j++)
        {
            this->elem[i * cols + j] = data;
        }
    }
}
void fMatrix::zeros(){
    fill(0);
}
void fMatrix::ones(){
    fill(1);
}
fMatrix  Inverse_LU(const fMatrix &a) {
    int n = a.rows;
    int m = a.cols;

    fMatrix p(a.rows, a.cols);
    fMatrix l(a.rows, a.cols);
    fMatrix u(a.rows, a.cols);
    fMatrix l_inv(l), u_inv(u);
    l_inv.zeros();
    u_inv.zeros();
    // Forward substitution to solve for invL
    for (int j = 0; j < n; j++) {
        l_inv.elem[j,j] = 1 / l.elem[j,j];
        for (int i = j + 1; i < n; i++) {
            double sum = 0;
            for (int k = j; k < i; k++) {
                sum += l.elem[i,k] * l_inv.elem[k,j];
            }
            l_inv.elem[i,j] = -sum / l.elem[i,i];
        }
    }
    // Backward substitution to solve for invU
    for (int j = n - 1; j >= 0; j--) {
        u_inv.elem[j,j] = 1 / u.elem[j,j];
        for (int i = j - 1; i >= 0; i--) {
            double sum = 0;
            for (int k = j; k > i; k--) {
                sum += u.elem[i,k] * u_inv.elem[k,j];
            }
            u_inv.elem[i,j] = -sum / u.elem[i,i];
        }
    }
    
    return  u_inv * l_inv;
}
// fMatrix Inverse_QR(const fMatrix &a)
// {
//     fMatrix result(a.rows, a.cols);
//     int n = a.rows;
//     double b[n][n]={0};
//     double t[n][n]={0};
//     int i;//数组  行
//     int j;//数组  列
//     int k;//代表B的角标
//     int l;//数组  列
//     double dev;
//     double numb;//计算的中间变量
//     double numerator,denominator;
//     double ratio;
//     double Q[n][n]={0};//正交矩阵
//     double QH[n][n]={0};//正交矩阵的转置共轭
//     double R[n][n]={0};//
//     double invR[n][n]={0};//R的逆矩阵
//     double invA[n][n]={0};//A的逆矩阵，最终的结果
// //={0};//
//     double matrixR1[n][n]={0};
//     double matrixR2[n][n]={0};
//     fMatrix q(a.rows, a.cols); //Q是标准正交矩阵
//     fMatrix r(a.rows, a.cols); //R是一个上三角矩阵

//     cout << "A:" << endl; //输出矩阵A
//     for (int i = 0; i < a.rows; i++)
//     {
//         for (int j = 0; j < a.cols; j++)
//         {
//             cout <<a.elem[i * a.cols + j] << " ";
//         }
//         cout << endl;
//     }
//     for(int i=0;i<n;++i)
//     {
//         for(int j=0;j<n;++j)
//         {
//             b[j][i]=a.elem[j,i];
//         }
//         for(int k=0;k<i;++k)
//         {
//             if(i)
//             {
//                 numerator=0.0;
//                 denominator=0.0;
//                 for(l=0;l<n;++l)
//                 {
//                     numerator+=a.elem[l,i]*b[l][k];
//                     denominator+=b[l][k]*b[l][k];
//                 }
//                 dev=numerator/denominator;
//                 t[k][i]=dev;
//                 for(j=0;j<n;++j)
//                 {
//                     b[j][i]-=t[k][i]*b[j][k];//t  init  =0  !!!
//                 }
//             }
//         }
//     }
//     for(i=0;i<n;++i)
//     {
//         numb=0.0;
//         for(j=0;j<n;++j)
//         {
//             numb+=(b[j][i]*b[j][i]);
//         }
//         dev=sqrt(numb);
//         for(j=0;j<n;++j)
//         {
//             Q[j][i]=b[j][i]/dev;
//         }
//         matrixR1[i][i]=dev;
//     }
//     for(i=0;i<n;++i)
//     {
//         for(j=0;j<n;++j)
//         {
//             if(j<i)
//             {
//                 matrixR2[j][i]=t[j][i];
//             }
//             else if(j==i)   
//             {
//                 matrixR2[j][i]=1;
//             }
//             else
//             {
//                 matrixR2[j][i]=0;
//             }
//         }
//     }
// ///////////////////////QR分解完毕//////////////////////////
//     printf("QR分解：\n");
//     printf("Q=\n");
//     for(i=0;i<n;++i)
//     {
//         for(j=0;j<n;++j)
//         {
//             printf("%2.4f    ",Q[i][j]);
//         //  
//         }
//         printf("\n");
//     }
//     printf("R=\n");
//     for(i=0;i<n;++i)
//     {
//         for(j=0;j<n;++j)
//         {
//             printf("%2.4f    ",R[i][j]);
//         //  
//         }
//         printf("\n");
//     }
//     return 0;
// }
// float DotProduct(const fMatrix &a, int rowA, const fMatrix &b, int rowB) {
//     float result = 0;
//     for (int i = 0; i < a.numCols(); i++) {
//         result += a(rowA, i) * b(rowB, i);
//     }
//     return result;
// }

// float Magnitude(const fMatrix& a, int row) {
//     float result = 0;
//     for (int i = 0; i < a.numCols(); i++) {
//         result += a(row, i) * a(row, i);
//     }
//     return sqrt(result);
// }
// fMatrix HouseholderReflect(const fMatrix &x) {
//     fMatrix result(x.rows, x.cols);
//     float mag = Magnitude(x, 0);
//     result.elem[0,0] = 1 - 2 * x.elem[0,0] * x.elem[0,0] / (mag * mag);
//     for (int i = 1; i < x.rows; i++) {
//         result.elem[0,i] = -2 * x.elem[0,i] * x.elem[0,0] / (mag * mag);
//     }
//     for (int i = 1; i < x.cols; i++) {
//         for (int j = 0; j < x.rows; j++) {
//             if (i == j) {
//                 result.elem[i,j] = 1;
//             } else {
//                 result.elem[i,j] = 0;
//             }
//         }
//     }
// }
// fMatrix Inverse_QR(const fMatrix &a) {
//     int j = 0;
//     fMatrix q(a.rows, a.cols);
//     for (int i = 0; i < a.rows; i++) {
//         for (int j = 0; j < a.cols; j++) {
//             if (i == j) {
//                 q.elem[i, j] = 1;
//             } else {
//                 q.elem[i, j] = 0;
//             }
//         }
//     }
//     for (int k = 0; k < a.cols && k < a.rows - 1; k++) {
//         fMatrix x(a.rows - k, 1);
//         for (int i = k; i < a.rows; i++) {
//             x.elem[i - k, 0] = a.elem[i, k];
//         }
//         fMatrix H = HouseholderReflect(x);
//         fMatrix temp = H * a;
//         for (int i = k; i < a.rows; i++) {
//             for (int j = k; j < a.cols; j++) {
//                 a.elem[i,j] = temp.elem[i - k, j - k];
//             }
//         }
//         for (int i = a.cols - 1; i >= k; i--) {
//             fMatrix x(q.rows - k, 1);
//             for (int j = k; j < q.rows; j++) {
//                 x.elem[j-k,0]= q.elem[j,i];
//             }
//             fMatrix H = HouseholderReflect(x);
//             fMatrix temp = H * q;
//             for (int j = k; j < q.rows; j++) {
//                 q.elem[j,i] = temp.elem[j-k,i];
//             }
//         }
//     }
//     fMatrix r = a;
//     for (int i = 0; i < a.rows; i++) {
//         for (int j = 0; j < i && j < a.cols; j++) {
//             r.elem[i,j] = 0;
//         }
//     }
//     fMatrix result(a.rows, a.cols);
//     for (int i = a.cols - 1; i >= 0; i--) {
//         fMatrix b(a.rows, 1);
//         b.elem[i,0] = 1;
//         for (int j = i + 1; j < a.cols; j++) {
//             float dot = DotProduct(r, i, result, j);
//             for (int k = 0; k < a.rows; k++) {
//                 b.elem[k,0] -= dot * result.elem[k,j];
//             }
//         }
//         float mag = Magnitude(r, i);
//         for (int k = 0; k < a.rows; k++) {
//             result.elem[k,i] = b.elem[k,0] / mag;
//         }
//     }
//     return result * Transp(q);
// }


fVector  Mean		( const fMatrix &A )
{
    fVector c(A.rows);
    Float n;
    for(int i=0;i<A.rows;i++)
    {
        n = 0;
        for(int j=0;j<A.cols;j++)
        {
            n += A.elem[i+j*A.cols];
            // cout << A.elem[i*A.cols+j] << " ";
        }
        // cout << endl;
        c.Array()[i] = n / A.cols;
    }
    return c;
}
fMatrix  Cov    ( const fMatrix &A )
{
    fVector Am = Mean(A);
    fMatrix c(A);
    for(int i=0;i<A.cols;i++)
    {
        for(int j=0;j<A.cols;j++)
        {
            c.elem[i*A.cols+j] -= Am.Array()[j]; 
        }
    }
    return (Transp(c)*c) / (A.rows-1);
}

fMatrix  Cov	( const fVector &A )
{
    return Cov(Outer(A,A)); 
}

void fMatrix::Show() const
{
    for(int i = 0; i < rows; i++)
    {
        for(int j = 0; j < cols; j++)
        {
            cout << elem[ i * cols + j] << " ";
        }
        cout << endl;
    } 
}

fMatrix  operator +  ( const fMatrix &a, const fMatrix &b)
{
    fMatrix c(a.rows,a.cols);
    for(int i=0; i < a.rows * a.cols; i++)
    {
        c.elem[i] = a.elem[i] + b.elem[i];
    }
    return c;
}
fMatrix  operator -  ( const fMatrix &a                 )
{
    fMatrix c(a.rows,a.cols);
    for(int i = 0; i < a.rows * a.cols; i++)
    {
        c.elem[i] = -a.elem[i];
    }
    return c;
}
fMatrix  operator -  ( const fMatrix &a, const fMatrix &b)
{
    fMatrix c(a.rows,a.cols);
    for(int i=0; i < a.rows * a.cols; i++)
    {
        c.elem[i] = a.elem[i] - b.elem[i];
    }
    return c;
}
fMatrix  operator *  ( const fMatrix &a,       Float n   )
{
    fMatrix c(a.rows,a.cols);
    for(int i=0; i < a.rows * a.cols; i++)
    {
        c.elem[i] = n *  a.elem[i];
    }
    return c;
}
fMatrix  operator *  (       Float n ,  const fMatrix &a )
{
    fMatrix c(a.rows,a.cols);
    for(int i = 0; i < a.rows * a.cols; i++)
    {
        c.elem[i] = n * a.elem[i];
    }
    return c;
}
fMatrix  operator /  (  const fMatrix &a,   Float   n )
{
    fMatrix c(a.rows,a.cols);
    for(int i = 0; i < a.rows * a.cols; i++)
    {
        c.elem[i] = a.elem[i]/n;
    }
    return c;
}
fMatrix  operator *  ( const fMatrix &A, const fMatrix &B )
{
    fMatrix c(A.rows,B.cols);
    if(A.cols != B.rows)
    {
        cout<<"shape error !!!"<<endl;
        return 0;
    }
    for(int i=0;i<A.rows*B.cols;i++)
    {
        Float sum = 0;
        int d = i/B.cols;
        int n = i%B.cols;

        for(int j=0;j<A.cols;j++)
        {
            sum += A.elem[d*A.cols+j]*B.elem[n+j*B.cols];
        }
        
        c.elem[i] = sum;
    }
    // cout << A.rows*B.cols<<endl;
    return c;
}
fVector  operator *  ( const fMatrix &a, const fVector &b )
{
    fVector c(a.rows);
    for(int i=0;i<a.rows;i++)
    {
        double sum = 0;
        for(int j=0;j<b.Size();j++)
        {
            sum += a.elem[i*b.Size()+j]*b.Array()[j];
        }
        c.Array()[i] = sum;
    }
        // cout<< B.Array()[i]<<endl;
    return c;
}
fVector  operator *  ( const fVector &a, const fMatrix &b )
{
    fVector c(b.cols);
    for(int i = 0;i < b.rows;i++)
    {
        double sum = 0;
        for(int j = 0; j < a.Size(); j++)
        {
            sum += a.Array()[j] * b.elem[i+j * b.cols];
        }
        c.Array()[i] = sum;
    }
    return c;
}
fMatrix& operator += (fMatrix &a, const fMatrix &b )
{
    for(int i = 0;i < a.rows * a.cols; i++)
    {
        a.elem[i] += b.elem[i];
    }
    return a;
}

fMatrix& operator -= (fMatrix &a, const fMatrix &b )
{
    for(int i = 0;i < a.rows * a.cols; i++)
    {
        a.elem[i] -= b.elem[i];
    }
    return a;
}

fMatrix& operator *= (fMatrix &A,  Float n )
{
    for(int i=0;i<A.rows*A.cols;i++)
    {
        A.elem[i] *= n;
    }
    return A;
}

fMatrix& operator /= (fMatrix &a,  Float n )
{
    for(int i = 0; i < a.rows * a.cols; i++)
    {
        a.elem[i] /= n;
    }
    return a;
}

fMatrix& operator *= (fMatrix &a, const fMatrix &b )
{   
    if(a.cols != b.rows)
    {
        cout<<"不能乘"<<endl;
        exit(1);
    }
    fMatrix c(a.rows, b.cols);
    for(int i = 0; i < a.rows * b.cols; i++)
    {
        Float sum = 0;
        int d = i / b.cols;
        int n = i % b.cols;

        for(int j=0;j<a.cols;j++)
        {
            sum += a.elem[ d * a.cols + j] * b.elem[ n + j * b.cols];
        }
        c.elem[i] = sum; 
    }
    a = c;
    
    return a;
}

fVector& operator *= (fVector &a, const fMatrix &b )
{
    fVector c(a.Size());
    for(int i = 0; i < b.cols; i++)
    {
        Float sum = 0;
        for(int j = 0;j < a.Size(); j++)
        {   
            sum += a.Array()[j] * b.elem[j*3+i];
        }
        c.Array()[i] = sum;
    }
    a=c;
    return a;
}