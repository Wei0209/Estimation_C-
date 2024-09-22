#include<iostream>
#include "ParamEstimator.h"
#include "fMatrix.h"
// #include "fVector.h"
// #include "fMatrix.cpp"
using namespace std;

CParamEstimator::CParamEstimator()
{
    LSP = new LS_Param;
    WLSP = new WLS_Param;
    MLP = new ML_Param;

    LSP -> pMat_H = new fMatrix;
    LSP -> pVec_Z = new fVector;
    LSP -> pMat_Vz = new fMatrix;

    WLSP -> pMat_H = new fMatrix;
    WLSP -> pVec_Z = new fVector;
    WLSP -> pMat_Vz = new fMatrix;
    WLSP -> pMat_W = new fMatrix;

    MLP -> pMat_H = new fMatrix;
    MLP -> pVec_Z = new fVector;
    MLP -> pMat_Vz = new fMatrix;
    MLP -> pMat_W = new fMatrix;
}
CParamEstimator::~CParamEstimator()
{
    delete LSP;
    delete WLSP;
    delete MLP;

    delete LSP->pMat_H;
    delete LSP->pVec_Z;
    delete LSP->pMat_Vz;

    delete WLSP->pMat_H;
    delete WLSP->pVec_Z;
    delete WLSP->pMat_Vz;
    delete WLSP->pMat_W;

    delete MLP->pMat_H;
    delete MLP->pVec_Z;
    delete MLP->pMat_Vz;


}

fMatrix*	CParamEstimator::SolveOptParam(fVector*	pfVecOptParam)
{
    if(this -> EstiMethod == LS) //Z(m*1)=H(m*n)*theta(n*1)
    {
        fMatrix H((*LSP->pVec_Z).Size(),(*LSP->pVec_Z).Size());
        fVector v((*LSP->pVec_Z).Size());
        *pfVecOptParam = (Inverse(ATranspA(*(LSP->pMat_H)))*Transp(*(LSP->pMat_H))*(*(LSP->pVec_Z)));
        // ((Ht * H)inverse) * Ht * Z
        pfVecOptParam->Show();
        v = (*LSP->pVec_Z) - (*LSP->pMat_H)*(*pfVecOptParam);
    }
    else if(this -> EstiMethod == WLS) 
    {
        fVector V((*WLSP -> pVec_Z).Size());
        fMatrix Vv((*WLSP  -> pVec_Z).Size(),(*WLSP -> pVec_Z).Size());
        fMatrix Inv_ATAH(Inverse(ATranspA((*WLSP -> pMat_H))));

        *pfVecOptParam = ( Inv_ATAH * Transp((*WLSP -> pMat_H)) * ((*WLSP -> pVec_Z)) );
        // pfVecOptParam->Show();
        V = (*WLSP -> pVec_Z) - ((*WLSP -> pMat_H)*(*pfVecOptParam));
        // V.Show();
        // Vv = Cov(V);
        // Vv.Show();
        // Vv = Diag(Diag(Vv));
        // Vv.Show();
        *WLSP -> pMat_W = Inverse(Vv);
        fMatrix Inv_ATAVH(Inverse((Transp(*WLSP -> pMat_H))*(*WLSP -> pMat_W)*(*WLSP -> pMat_H)));
        // ((Transp(*this->WLSP.pMat_H))).Show();
        // Inv_ATAVH.Show();
        *pfVecOptParam = Inv_ATAVH * Transp(*WLSP->pMat_H) * (*WLSP->pMat_W) * (*(WLSP -> pVec_Z));
        // pfVecOptParam->Show();
    }
    else if(this -> EstiMethod == ML)
    {
        fMatrix H((*LSP->pVec_Z).Size(),(*LSP->pVec_Z).Size());
        fVector v((*LSP->pVec_Z).Size());
        *pfVecOptParam = (Inverse(ATranspA(*(MLP->pMat_H))))  *  Transp(*(MLP->pMat_H)) * (*(MLP->pVec_Z));
        //((Ht * w * H)inverse) * Ht * w * Z
        pfVecOptParam->Show();
        v = (*MLP->pVec_Z) - ((*MLP->pMat_H)*(*pfVecOptParam));
        // H = Cov(v);
        // H = Diag(Diag(H));
        *MLP->pMat_W = Inverse(H);
        fMatrix Inv_ATAVH(Inverse((Transp(*MLP->pMat_H))*(*MLP->pMat_W)*(*MLP->pMat_H)));
        ((Transp(*MLP->pMat_H))).Show();
        Inv_ATAVH.Show();
        *pfVecOptParam = Inv_ATAVH * Transp(*MLP->pMat_H) * (*MLP->pMat_W) * (*(MLP->pVec_Z));
        pfVecOptParam -> Show();
    }
}
void CParamEstimator::SetParamEstiMethod(ParamEstiMethod Method)
{
    EstiMethod = Method;
}
void CParamEstimator::SetMethodParameters(ParamEstiMethod Method, void*	pParam)
{
    if(Method == LS)
    {
        LSP = ((LS_Param *)pParam);
    }
    if(Method == WLS)
    {
        WLSP =((WLS_Param *)pParam);
    }
    if(Method == ML)
    {
        MLP = ((ML_Param *)pParam);
    }
}
ParamEstiMethod	CParamEstimator::GetParamEstiMethod(void) const
{
    return ParamEstiMethod(EstiMethod);
}