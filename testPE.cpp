#include<iostream>
#include "fVector.h"
#include "fMatrix.h"
#include "ParamEstimator.h" 
using namespace std;

void testVectorFuns()
{
    Float A[3] = {1.1, 2.2, 3.3};
	Float B[3] = {4.4, 5.5, 6.6};

	fVector VecA(A, 3);
	fVector VecB(B, 3);

	cout << "\nVecA = " << endl;
	VecA.Show();
	
	cout << "\nVecB = " << endl;
	VecB.Show();
	
	cout << "\nStarts to test vector operators..." << endl;
	// 1. A+B
	cout << "\n1. A+B" << endl;
	(VecA+VecB).Show();

	// 2. A-B
	cout << "\n2. A-B" << endl;
	(VecA-VecB).Show();
    
    // 3. -A
	cout << "\n3. -A" << endl;
	(-VecA).Show();

    // 4. A-2
	cout << "\n4. A-2" << endl;
	(VecA-2).Show();

    // 5. 2-A
	cout << "\n5. 2-A" << endl;
	(2-VecA).Show();

    // 6. A*2
	cout << "\n6. A*2" << endl;
	(VecA*2).Show();

    // 7. 0.5*B
	cout << "\n7. 0.5*B" << endl;
	(0.5*VecB).Show();
    
    // 8. A/4
	cout << "\n8. A/4" << endl;
	(VecA/4).Show();

    // 9. A/B
	cout << "\n9. A/B" << endl;
	(VecA/VecB).Show();

    // 10. A*B
	cout << "\n10. A*B" << endl;
	cout << "A*B = " << VecA*VecB << endl;

	// 11. A^B
	cout << "\n11. A^B" << endl;
    (VecA^VecB).Show();

    // 12. A+=B
	cout << "\n12. A+=B" << endl;
	(VecA+=VecB).Show();

    // 13. A-=B
	cout << "\n13. A-=B" << endl;
	(VecA-=VecB).Show();

    // 14. A*=5
	cout << "\n14. A*=5" << endl;
	(VecA*=5).Show();

    // 15. A/=5
	cout << "\n15. A/=5" << endl;
	(VecA/=5).Show();

    //Copy Vector
    cout << "\nCopy Vector A -> C" <<endl;
    fVector VecC(VecA);
    VecC.Show();

	cout << "\nStarts to test vector functions..." << endl;

	// 16. Min(A,B)
	cout << "\n16. Min(A,B)" << endl;
	(Min(VecA,VecB)).Show();

	// 17. Max(A,B)
	cout << "\n17. Max(A,B)" << endl;
	(Max(VecA,VecB)).Show();

	// VecA.Show();
	// VecB.Show();

	// 18. Dist(A,B)
	cout << "\n18. Dist(A,B)" << endl;
	cout << "Dist(A,B) = " << Dist(VecA,VecB) << endl;

	// 19. Normalize(A,B)
	cout << "\n19. Normalize(A)" << endl;
	(Normalize(VecA)).Show();

	// 20. OneNorm(A)
	cout << "\n20. OneNorm(A)" << endl;
	cout << "OneNorm(A) = " << OneNorm(VecA) << endl;

	// 21. TwoNorm(A)
	cout << "\n21. TwoNorm(A)" << endl;
	cout << "TwoNorm(A) = " << TwoNorm(VecA) << endl;

	// 22. TwoNormSqr(A)
	cout << "\n22. TwoNormSqr(A)" << endl;
	cout << "TwoNormSqr(A) = " << TwoNormSqr(VecA) << endl;

	// 23. Sqrt(A,B)
	cout << "\n23. Sqrt(A)" << endl;
	(Sqrt(VecA)).Show();

	// 24. Mean(A)
	cout << "\n24. Mean(A)" << endl;
	cout << "Mean(A) = " << Mean(VecA) << endl;

	// 25. Var(A)
	cout << "\n25. Var(A)" << endl;
	cout << "Var(A) = " << Var(VecA) << endl;

	// 26. Std(A)
	cout << "\n26. Std(A)" << endl;
	cout << "Std(A) = " << Std(VecA) << endl;

}

void testMatrixFuns()
{
	Float A[3] = {1.1, 2.2, 3.3};
	Float B[3] = {4.4, 5.5, 6.6};
	Float C[9] = {0.9649, 0.9572, 0.1419, 
				  0.1576, 0.4854, 0.4218, 
				  0.9706, 0.8003, 0.9157};
	Float D[9] = {0.8147, 0.9134, 0.2785, 
				  0.9058, 0.6324, 0.5469, 
				  0.1270, 0.0975, 0.9575};
	Float X[15] = {0.7922, 0.9340, 0.6555, 
				   0.9595, 0.6787, 0.1712, 
				   0.6557, 0.7577, 0.7060, 
				   0.0357, 0.7431, 0.0318, 
				   0.8491, 0.3922, 0.2769};
	// Float s[4] = {1,2,3,4};			   
	fVector VecA(A, 3);
	fVector VecB(B, 3);
	fVector VecC(VecA);
	fMatrix MatA(C, 3, 3);
	fMatrix MatB(D, 3, 3);
	fMatrix MatC(MatA);
	fMatrix MatX(X, 5, 3);
	// fMatrix Mats(s, 2, 2);
	cout << "\nMatA = " << endl;
	MatA.Show();
	
	cout << "\nMatB = " << endl;
	MatB.Show();

	cout << "\nMatX = " << endl;
	MatX.Show();

	cout << "\nStarts to test matrix operators..." << endl;
	// 1. A+B
	cout << "\n1. A+B" << endl;
	(MatA+MatB).Show();

	// 2. A-B
	cout << "\n2. A-B" << endl;
	(MatA-MatB).Show();

	// 3. -B
	cout << "\n3. -B" << endl;
	(-MatB).Show();

	// 4. 2*A
	cout << "\n4. 2*A" << endl;
	(2*MatA).Show();

	// 5. A*0.5
	cout << "\n5. A*0.5" << endl;
	(MatA*0.5).Show();

	// 6. A/2
	cout << "\n6. A/2" << endl;
	(MatA/2).Show();

	// cout << "\n MatX" << endl;
	// MatX.Show();

	// 7. A*B
	cout << "\n7. A*B" << endl;
	(MatA*MatB).Show();

	// 8. MatA*VecA
	cout << "\n8. MatA*VecA" << endl;
	(MatA*VecA).Show();

	// 9. VecA*MatA
	cout << "\n9. VecA*MatA" << endl;
	(VecA*MatA).Show(RowVec);

	// 10. A+=B
	cout << "\n10. A+=B" << endl;
	(MatA+=MatB).Show();

	// 11. A-=B
	cout << "\n11. A-=B" << endl;
	(MatA-=MatB).Show();

	// 12. A*=2
	cout << "\n12. A*=2" << endl;
	(MatA*=2).Show();

	// 13. A/=2
	cout << "\n13. A/=2" << endl;
	(MatA/=2).Show();
	cout << "\nStarts to test matrix functions..." << endl;
	
	// 14. A*=B
	cout << "\n14. A*=B" << endl;
	(MatA*=MatB).Show();
	MatA = MatC;
	
	// cout << "\n MatA=" << endl;
	// MatA.Show();
	// 15. VecA*=MatA
	cout << "\n15. VecA*=MatA" << endl;
	(VecA*=MatA).Show(RowVec);
	VecA = VecC;

	// MatA.Show();
	// VecC.Show();

	// MatA.Show();
	// 16. Xt
	cout << "\n16. Xt" << endl;
	Transp(MatX).Show();

	// 17. X*Xt
	cout << "\n17. X*Xt" << endl;
	AATransp(MatX).Show();

	// 18. Xt*X
	cout << "\n18. Xt*X" << endl;
	ATranspA(MatX).Show();

	// 22, Diag(MatA)
	cout << "\n22. Diag(MatA)" << endl;
	(Diag(MatA)).Show();

	// 24. Determinant(MatA)
	cout << "\n24. Determinant(MatA) = " << Determinant(MatA) << endl;

	// 25. Trace(MatA)
	cout << "\n25. Trace(MatA) = " << Trace(MatA) << endl;

	// 28. Inverse(MatA)
	cout << "\n28. Inverse(MatA) = " << endl; 
	(Inverse(MatA)).Show();

	// cout << "\nMatA"<<endl;
	// MatA.Show();

	// 30. Mean(MatA)
	cout << "\n30. Mean(MatA) = " << endl; 
	(Mean(MatA)).Show(RowVec);

	// 31. Cov(MatA)
	cout << "\n31. Cov(MatA) = " << endl; 
	(Cov(MatA)).Show();


	
	// Float s[25] = {1,8,1,6,7,5,4,3,1,4,9,6,4,1,2,3,5,5,6,8,9,7,5,12,5};
	// fMatrix Mats(s, 5, 5);

	// cout << "\n5x5 det = " << endl << Determinant(Mats) <<endl;

	// cout << "\n28. Inverse(Mats) = " << endl; 
	// (Inverse(Mats)).Show();

	// cout << MatA(2,0) << endl;
}

void testParamEstimator()
{
	// Float A[5] = {  2.0, 2.1, 2.1, 2.03, 2.04 };
	// Float C[10] = { 1, 2.6,
	// 				1, 2.72,
	// 				1, 2.75,
	// 				1, 2.67,
	// 				1, 2.68 };

	// Float d[9] = {	0.9649, 0.9572, 0.1419,
	// 				0.1576, 0.4854, 0.4218,
	// 				0.9706, 0.8003, 0.9157 };

	// fVector VecA(A, 5);
	// fMatrix matA(C,5,2);
	// fMatrix matd(d,3,3);
	// fVector VecZ(2);
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

	// c.SetParamEstiMethod(LS);
	// cout << c.GetParamEstiMethod() << endl;
	// LS_Param param = {&matX,&vecY,NULL};
	// c.SetMethodParameters(LS,&param);
	// c.SolveOptParam(&vecZ);

	c.SetParamEstiMethod(WLS);
	cout << c.GetParamEstiMethod() << endl;
	// c.SetMethodParameters(LS,&{matA,VecA,matA});
	WLS_Param param = {&matX,&matW,&vecY, NULL};	
	c.SetMethodParameters(WLS,&param);
	c.SolveOptParam(&vecZ);

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