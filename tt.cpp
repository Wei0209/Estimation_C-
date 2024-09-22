#include<iostream>
#include <stdlib.h>
#include"fVector.h"
#include <vector>
#include"fMatrix.h"
#include "ParamEstimator.h" 
using namespace std;

void testParamEstimator()
{
	Float A[5] = {  2.0, 2.1, 2.1, 2.03, 2.04 };
	Float C[10] = { 1, 2.6,
					1, 2.72,
					1, 2.75,
					1, 2.67,
					1, 2.68 };

	Float d[9] = {	0.9649, 0.9572, 0.1419,
					0.1576, 0.4854, 0.4218,
					0.9706, 0.8003, 0.9157 };

	fVector VecA(A, 5);
	fMatrix matA(C,5,2);
	fMatrix matd(d,3,3);
	fVector VecZ(2);
	Float x[50] = {
		1,50,
		1,53,
		1,54,
		1,55,
		1,56,
		1,59,
		1,62,
		1,65,
		1,67,
		1,71,
		1,72,
		1,74,
		1,75,
		1,76,
		1,79,
		1,80,
		1,82,
		1,85,
		1,87,
		1,90,
		1,93,
		1,94,
		1,95,
		1,97,
		1,100,
	};
	Float y[25] = {
		122,
		118,
		128,
		121,
		125,
		136,
		144,
		142,
		149,
		161,
		167,
		168,
		162,
		171,
		175,
		182,
		180,
		183,
		188,
		200,
		194,
		206,
		207,
		210,
		219,
	};
	fMatrix matX(x,25,2);
	fMatrix matW(25,25);
	fVector vecY(y,25);
	fVector vecZ(25);
	CParamEstimator c;

	c.SetParamEstiMethod(LS);
	cout << c.GetParamEstiMethod() << endl;
	LS_Param param = {&matX,&vecY,NULL};
	c.SetMethodParameters(LS,&param);
	c.SolveOptParam(&vecZ);

	// c.SetParamEstiMethod(WLS);
	// cout << c.GetParamEstiMethod() << endl;
	// // c.SetMethodParameters(LS,&{matA,VecA,matA});
	// WLS_Param param = {&matX,&matW,&vecY,NULL};
	// c.SetMethodParameters(WLS,&param);
	// c.SolveOptParam(&vecZ);

	// cout << "\n" ;
	// fMatrix A(3,4);
	// A.Show();

	// cout << "\n";
	// Float B[4] = {1,1,1,1};
	// fVector B_(B,4);
	// A.SetRow(1,B_);
	// A.SetRow(2,B_*3);
	// A.Show();
}

int main()
{
	// testVectorFuns();
	// testMatrixFuns();
	testParamEstimator();
	return 0;
}