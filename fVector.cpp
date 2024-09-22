#include <iostream>
#include "fMatrix.h"
#include "fVector.h"
#include <cstring>
#include <math.h>
#include <stdlib.h>
using namespace std;

// Initinalize constructor.
fVector::fVector(int size)
{
    if (size != 0)
    {
        this -> size = size;//將Vector物件元素數量做設定
        this-> elem = new Float[this-> size];//將元素所需空間作記憶體配置
        memset(this -> elem, 0, this -> size * sizeof(Float));
    }
    else
    {
        this -> size = 0;
        this -> elem =NULL;
    }
    nVecCount;
    // this->size = size;
    // this->elem = new Float[size];
}
// Copy constructor.
fVector::fVector( const fVector &a )
{
    if(a.size != 0)
    {
        this -> size = a.size;
        this -> elem = new Float[this -> size];
        memcpy(this -> elem, a.elem, this -> size * sizeof(Float));
    }
    else
    {
        this -> size = 0;
        this -> elem =NULL;
    }
    nVecCount;
}
// Assign constructor.
fVector::fVector( const Float *x, int n )
{
    if (n != 0 && x)
    {
        this->size = n;
        this->elem = new Float[this->size];
        memcpy(this->elem, x, this->size * sizeof(Float));
    }
    else
    {
        this->size = 0; 
        this->elem = NULL; 
    } 
    nVecCount;
}
fVector::fVector( int n, const Float *x )
{
    if (n != 0 && x)
    {
        this -> size = n;
        this->elem = new Float[this->size];
        memcpy(this->elem, x, this->size * sizeof(Float));
    }
    else
    {
        this->size = 0; 
        this->elem = NULL;
    }
    nVecCount;
}
fVector::fVector( Float x, Float y)
{
    this->size = 2;
    this->elem = new Float[this->size];
    this->elem[0] = x; 
    this->elem[1] = y;
    nVecCount;

}
fVector::fVector( Float x, Float y, Float z)
{
    this->size = 3;
    this->elem = new Float[this->size];
    this->elem[0] = x;
    this->elem[1] = y; 
    this->elem[2] = z;  
    nVecCount;
}
fVector &fVector::operator=(const fVector &Vec)
{
    if (this != &Vec)
    {
        if (this->size != Vec.size)
        {
            if (this->elem)
                delete [] this->elem;
            if (Vec.size != 0)
            {
                this->size = Vec.size;
                this->elem = new Float[this->size];
                memcpy(this->elem, Vec.elem, this->size * sizeof(Float));
            }
            else
            {
                this->size = 0;
                this->elem = NULL;
            }
        }
        else 
        memcpy(this->elem, Vec.elem, this->size * sizeof(Float));
    }
    return *this;
}
void fVector::operator=(Float Val)
{
    if (this->size != 0)
    {
        for (int i = 0; i < this->size; ++i)
        {
            this->elem[i] = Val;
        }
    }
}
void fVector::SetSize(int size)
{
    if (this->elem)
        delete[] this->elem; 
    if (size != 0)
    {
        this->size = size; 
        this->elem = new Float[this->size]; 
    }
    else
    {
        this->size = 0; 
        this->elem = NULL; 
    }
}
fVector &fVector::Swap(int i , int j)
{
    if(this -> size > i && this -> size > j)
    {
        Float temp = this->elem[i];
        this -> elem[i] = this->elem[j];
        this -> elem[j] = temp;
    }
    return *this;
}

fVector fVector::GetBlock(int i, int j) const
{
    fVector oVec;
    if (j >= i && this->size > j)
    {
        int BlockSize = j - i + 1; //算出要提取的元素數量
        oVec.SetSize(BlockSize); //製作一個新的Vector
        
        memcpy(oVec.elem, this->elem + i, sizeof(float)*BlockSize);
    }
    return oVec;
}
void fVector::SetBlock(int i, int j, const fVector &Vec)
{
    if (j >= i && Vec.Size() == (j - i + 1)) //確定在合理範圍
    {
        if (this->size > j) //判斷所要交換的值是否超過元素數量
        {
            for (int addr = 0; addr < Vec.Size(); ++addr)

            {
                *(this->elem + addr + i) = Vec(addr); //將區塊元素值置入
            }
        }
    }
}

fVector::~fVector()
{
    if (this->elem)
    {
        delete[] this->elem; 
    }
    nVecCount;
}
fVector operator+( const fVector &a, const fVector &b )
{
    fVector c(a.Size());
    for(int i = 0; i < a.Size(); i++)
    {
        c(i) = a(i) + b(i);
        
    }
    return c;
}
fVector  operator -  ( const fVector &a, const fVector &b )// Binary minus.
{
    fVector c(a.Size());
	for(int i = 0; i < a.Size(); i++)
    {
        c(i) = a(i) - b(i);
    }
	return c;
}
fVector  operator -  ( const fVector &a         ) // Unary minus.
{
    fVector c(a.Size());
	for(int i = 0; i < a.Size(); i++) 
    {
        c(i) = 0 - a(i);
    }
	return c;
}
fVector  operator -  ( const fVector &a, Float b)
{
    fVector c(a.Size());
	
	for(int i = 0; i < a.Size(); i++)
    {
        c(i) = a(i) - b;
    }
	return c;
}
fVector  operator -  (		Float b	, const fVector &a )
{
    fVector c(a.Size());
	
	for(int i = 0; i < a.Size(); i++)
    {
        c(i) = b - a(i);
    }
	return c;
}
fVector  operator *  ( const fVector &a,        Float b )
{
    fVector c(a.Size());
    for(int i = 0; i < a.Size(); i++)
    {
        c(i) = a(i) * b;
    }
	return c;
}

fVector  operator *  (       Float b   , const fVector &a )
{
    fVector c(a.Size());
    for(int i = 0; i < a.Size(); i++)
    {
        c(i) = b * a(i);
    }
	return c;
}


fVector  operator /  ( const fVector &a,        Float b)
{
    fVector c(a.Size());
    for(int i = 0; i < a.Size(); i++) 
    {
        c(i) = a(i) / b;
    }
	return c;
}
fVector  operator /  ( const fVector &a, const fVector &b ) // Element-wise division
{
    fVector c(a.Size());
    if(a.Size() == b.Size())
    {
        for(int i = 0; i < a.Size(); i++)
        {
            c(i) = a(i) / b(i);
        }
    }
    return c;
}
double operator*(const fVector &a, const fVector &b) // Inner-product between two vectors
{
    double Val = 0.0;
    if (a.Size() == b.Size()) //判斷兩Vector元素數量是否相同
    {
        for (int i = 0; i < a.Size(); ++i)
        {
            Val += (double)a(i) * (double)b(i); //兩Vector元素相乘後相加(內積)
        }
    }
    return Val;
}
fVector  operator ^  ( const fVector &A, const fVector &B)
{
    if(A.size != 3 && B.size !=3)
    {
        cout << "vector shape error";
        return 0; 
    }
    fVector c(A.size);
    c.elem[0] = A.elem[1]*B.elem[2] - A.elem[2]*B.elem[1];
    c.elem[1] = -A.elem[0]*B.elem[2] + A.elem[2]*B.elem[0];
    c.elem[2] = A.elem[0]*B.elem[1] - A.elem[1]*B.elem[0];
    
    return c;
}
fVector& operator += (       fVector &a, const fVector &b )
{
    for(int i = 0; i < a.Size(); i++){
        a(i) = a(i) + b(i);
    }
    return a;
}
fVector& operator -= (       fVector &a, const fVector &b )
{
    for(int i = 0; i < a.Size(); i++){
        a(i) = a(i) - b(i);
    }
    return a;
}
fVector& operator *= (       fVector &a,        Float   b)
{
    for(int i = 0;i < a.Size(); i++){
        a(i) = a(i) * b;
    }
    return a;
}
fVector& operator /= (       fVector &a,        Float  b)
{
    for(int i = 0; i < a.Size(); i++){
        a(i) = a(i) / b;
    }
    return a;
}
fVector Min ( const fVector &A, const fVector &B )
{
    fVector c(A.size);
    for(int i=0;i<A.size;i++)
    {
        if(A.elem[i]>=B.elem[i])
            c.elem[i] = B.elem[i];
        if(A.elem[i]<B.elem[i])
            c.elem[i] = A.elem[i];
    }
    return c;
}

fVector Max ( const fVector &A, const fVector &B )
{
    fVector c(A.size);
    for(int i=0;i<A.size;i++)
    {
        if(A.elem[i]>=B.elem[i])
            c.elem[i] = A.elem[i];
        if(A.elem[i]<B.elem[i])
            c.elem[i] = B.elem[i];
    }
    return c;
}
double   Dist        ( const fVector &a, const fVector &b )// Returns two norm distance between two vectors
{
    double dis =0;
    if(a.Size() == b.Size())
    {
        for(int i = 0; i < a.Size(); i++) 
        {
            dis = dis + pow((a(i) - b(i)) , 2);
        }
    }
	return sqrt(dis);
} 

fVector  Normalize   ( const fVector &a ) // Normalizes a vector into an unit vector
{ 
	double count = 0;
	for(int i = 0 ; i < a.Size() ; i++) 
		count = count + pow(a(i) , 2);

	return sqrt(count);

}
double   OneNorm     ( const fVector &a ) // Returns one norm value of a vector
{
    double ElemSum = 0.0;
    if (a.Size() != 0)
    {
        for (int i = 0; i < a.Size(); ++i)
        {
            double ElemT = fabs((double)a(i));
            if (ElemSum < ElemT){
                ElemSum = ElemT;
            }
        }
    }
    return ElemSum;
}
double   TwoNormSqr     ( const fVector &a ) // Returns square of the two norm value of a vector 
{
    double ElemSum = 0.0;
    if (a.Size() != 0) //判斷元素數量是否為0
    {
        //Two Norm運算法則，但不開根號
        for (int i = 0; i < a.Size(); ++i)
        {
            ElemSum += (double)a(i) * (double)a(i);
        }
    }
    return ElemSum;
}
double   TwoNorm  ( const fVector &a )// Returns two norm value of a vector
{
    double ElemSum = 0.0;
    if (a.Size() != 0) //判斷元素數量是否為0
    {
    //Two Norm運算法則
        for (int i = 0; i < a.Size(); ++i)
        {
            ElemSum += (double)a(i) * (double)a(i);
        }
        ElemSum = sqrt(ElemSum);
    }
    return ElemSum;
}

fVector  Sqrt   ( const fVector &A )
{
    fVector c(A.size);
    for(int i=0;i<A.size;i++)
    {
        c.elem[i] = sqrt(A.elem[i]);
    }
    return c;
}
double   Mean		( const fVector &a )// Mean value of a vector.
{
    Float ElemSum = 0.0f;
    if (a.Size() != 0) //判斷元素數量是否為0
    {
        for (int i = 0; i < a.Size(); ++i)
        {
            ElemSum += a(i); //所有元素相加
        }
        ElemSum /= a.Size();//取平均
    }
    return ElemSum;
}   
double   Var			( const fVector &a )// Variance of a vector.
{
    double avg = 0 , Variance = 0;
	for(int i = 0 ; i < a.Size() ; i++) 
		avg = avg + a(i)/a.Size();
	for(int i = 0 ; i < a.Size() ; i++) 
		Variance = Variance + pow((a(i) - avg),2);

	return Variance;
}  
double   Std	( const fVector &a )// Standard derivation of a vector.  
{
    return sqrt(Var(a));
}  	
void fVector::Show(VecType Type) const {
    if (Type == ColVec) {
        // 列向量顯示
        for (int i = 0; i < size; i++) {
            cout << elem[i] << endl;
        }
    } else {
        // 行向量顯示
        for (int i = 0; i < size; i++) {
            cout << elem[i] << " ";
        }
        cout << endl;
    }
}


// void ShowVector(const fVector &iVec, VecType Type = ColVec )
// {
//     if (iVec.Size() != 0) //判斷元素數量是否為0
//     {
//         cout << "[ ";
//         for (int i = 0; i < iVec.Size(); ++i)
//         {
//             cout << iVec(i) << " ";
//         }
//         cout << "]\n";
//     }
// }