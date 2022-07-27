#include <math.h>
#include <vector>

#include "C_FEM_BasisFunction_1D.h"
#include "C_FEM_GaussPoint_1D.h" 
#include "C_Mesh_1D.h"
#include "C_Matrix_Sparse.h"
#include "C_Matrix_Dense.h"

#include "f_ForcingVector.h"

int main()
{
    C_Material mat;
    C_GaussPoint_1D GP_Data(2);
    C_Mesh_1D       mesh(10.0, 11);
    C_LagrangeBasis feL(1, GP_Data);
    int totDof=mesh.num_Nd;
    C_Matrix_Dense KuGlobal(totDof,1);
    C_Matrix_Dense uGlobal(totDof,1);
    C_Matrix_Dense fGlobal(totDof,1);

    for (int itEl = 0; itEl < mesh.num_El; itEl++) {
        calculateForceVector_AxialBar(itEl, mesh, mat, KuGlobal, uGlobal, feL, GP_Data);
    }

    // i. Neumann Boundary Conditions
    int fDof=totDof-1;
    fGlobal(fDof,0)=10.0;

    // ii. Dirichlet Boundary Conditions
    int bcDof=0;
    double bcValue = 0.;
    uGlobal(bcDof,0) = bcValue; 

    
    
    std::cout << uGlobal;
    std::cout << fGlobal;

    double m=1, c=500., dt=0.01, t=0.;
    C_Matrix_Dense v(totDof,1), a(totDof,1);

    while (t<5)
    {
        a = (1/m)*((fGlobal - KuGlobal) - c*v);
        v = v + (a*dt);
        uGlobal = uGlobal + (v*dt);
        // ii. Dirichlet Boundary Conditions
        int bcDof=0;
        double bcValue = 0.;
        uGlobal(bcDof,0) = bcValue;
        for (int itEl = 0; itEl < mesh.num_El; itEl++) {
            calculateForceVector_AxialBar(itEl, mesh, mat, KuGlobal, uGlobal, feL, GP_Data);
        }
        t = t + dt;
    }
    std::cout<<"U global: ";
    std::cout<<uGlobal<<"\n";
    std::cout<<"V global: ";
    std::cout<<v<<"\n";
    std::cout<<"A global: ";
    std::cout<<a;
    
    
    
    return 0;
}