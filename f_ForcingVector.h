#ifndef F_FORCINGVECTOR_H
#define F_FORCINGVECTOR_H

#include <math.h>
#include <iostream>
#include <tuple>

#include "C_FEM_BasisFunction_1D.h"
#include "C_FEM_GaussPoint_1D.h" 
#include "C_Mesh_1D.h"
#include "C_Mesh_Frame.h"
#include "C_Matrix_Dense.h"

#include "f_MiscellaneousFunctions.h"

//! Construct/Assemble 1D Axial Stiffness Matrix
/*!
This function computes the element matrices and positions them in the global matrix, kGlobal, for
the linear elastic stiffness.

    \param flags object containing material flags
    \param material object containing material data
    \param kGlobal Pass-by-reference global sparse matrix. Initialized in COO form.
    \param fGlobal Pass-by-reference forcing vector.
    \param itEl Current element index which is being initialized.
    \param elNodes Vector containing start and end locations (in x) of element nodes.
    \param elDOF Degree-of-Freedom element matrix.
    \param GP_Data Object containing Gauss-Points.

Author: Dominic Jarecki
Date: 4-14-2022
*/
class C_Material{
    public:
        // Elastic Constants
        double EI = 1;
        double EA = 1;
        
        // Elastic-Plastic Constants
        double my = 15.66; // Plastic Yield Moment
        double a  = 10.10; // Sharpness Parameter

    C_Material(){}
};

void calculateForceVector_AxialBar(int itEl, C_Mesh_1D& mesh, C_Material& mat, 
    C_Matrix_Dense& K_norm, C_Matrix_Dense& KuGlobal, C_Matrix_Dense& uGlobal, C_LagrangeBasis& feL, C_GaussPoint_1D& GP_Data){
    
    std::vector<int>      elDoF = mesh.elements[itEl];
    std::vector<double> elNodes = {mesh.nodes[elDoF[0]], mesh.nodes[elDoF[1]]};

    double le = elNodes[1] - elNodes[0];

    C_Matrix_Dense dsp_L;
    for (int itGp = 0; itGp < GP_Data.num_GP; itGp++) {

        int ns_L = feL.num_sp;
        
        // Extract Rows at each Gauss Point
        dsp_L = (2/le)*feL.dsp(itGp, intspace(0,ns_L));
        double JxW = 0.5*GP_Data.wt[itGp]*le;
        C_Matrix_Dense S11 = JxW*mat.EA*(dsp_L.T()*dsp_L);
        K_norm(itEl,0) = S11.norm();

        // Add to Global Force Vector
        C_Matrix_Dense tmp=uGlobal(elDoF,0);
        KuGlobal.add_matr(S11*tmp, elDoF, {0});
    }
}


void calculateForceVector_Truss(int itEl, C_Mesh_Frame& mesh, C_Material& mat, 
    C_Matrix_Dense& K_norm, C_Matrix_Dense& KuGlobal, C_Matrix_Dense& uGlobal, C_LagrangeBasis& feL, C_GaussPoint_1D& GP_Data){
    std::vector<int> elCon=mesh.elements[itEl];
    std::vector<int>      elDoF = {mesh.elements[itEl][0]*2, mesh.elements[itEl][0]*2+1, mesh.elements[itEl][1]*2, mesh.elements[itEl][1]*2+1};
    std::vector< std::vector<double> > elNodes = {mesh.nodes[elCon[0]], mesh.nodes[elCon[1]]};

    double le = sqrt(pow(elNodes[1][0] - elNodes[0][0], 2) + pow(elNodes[1][1] - elNodes[0][1], 2) + pow(elNodes[1][2] - elNodes[0][2], 2));
    double theta;
    double eps=1e-9;
    if ((elNodes[1][0] - elNodes[0][0])<=eps)
    {
        if (elNodes[0][1] < elNodes[1][1])
        {
            theta = M_PI/2.0;
        }
        else if (elNodes[0][1] > elNodes[1][1])
        {
            theta = -M_PI/2.0;
        }
        
    }
    else
    {
        theta = atan2((elNodes[1][1] - elNodes[0][1]), (elNodes[1][0] - elNodes[0][0]));
    }
    double cos = cosf64(theta);
    double sin = sinf64(theta);
    C_Matrix_Dense T(4,4);
    T(0,0) = cos;
    T(0,1) = sin;
    T(1,0) = -sin;
    T(1,1) = cos;
    T(2,2) = cos;
    T(2,3) = sin;
    T(3,2) = -sin;
    T(3,3) = cos;
    

    C_Matrix_Dense dsp_L;
    for (int itGp = 0; itGp < GP_Data.num_GP; itGp++) {

        int ns_L = feL.num_sp;
        
        // Extract Rows at each Gauss Point
        dsp_L = (2/le)*feL.dsp(itGp, intspace(0,ns_L));
        double JxW = 0.5*GP_Data.wt[itGp]*le;
        C_Matrix_Dense S11_local = JxW*mat.EA*(dsp_L.T()*dsp_L);
        C_Matrix_Dense S11_local_app(4,4);
        std::vector<int> loc = {0,2};
        S11_local_app.add_matr(S11_local, loc, loc);
        C_Matrix_Dense S11 = T.T()*S11_local_app*T;
        K_norm(itEl,0) = S11.norm();

        // Add to Global Force Vector
        C_Matrix_Dense tmp=uGlobal(elDoF,0);
        KuGlobal.add_matr(S11*tmp, elDoF, {0});
    }
}


#endif