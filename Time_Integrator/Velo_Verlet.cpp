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
    C_Matrix_Dense K_norm(mesh.num_El,1);
    C_Matrix_Dense KuGlobal(totDof,1);
    C_Matrix_Dense uGlobal(totDof,1);
    C_Matrix_Dense fGlobal(totDof,1);

    for (int itEl = 0; itEl < mesh.num_El; itEl++) {
        calculateForceVector_AxialBar(itEl, mesh, mat, K_norm, KuGlobal, uGlobal, feL, GP_Data);
    }

    // i. Neumann Boundary Conditions
    int fDof=totDof-1;
    fGlobal(fDof,0)=10.0;

    // ii. Dirichlet Boundary Conditions
    int bcDof=0;
    double bcValue = 0.;
    uGlobal(bcDof,0) = bcValue; 

    // Parameters for time integration scheme
    int sum = 0.0;
    for (int i = 0; i < K_norm.row_size; i++)
    {
        sum = sum + K_norm(i,0);
    }
    double k_avg = sum/K_norm.row_size;             // Request Dominic if he can introduce avg() method for Dense matrix class
    double c_ratio=0.95, time_cnst=1, w_n=1/(time_cnst*c_ratio);         // At t=time_cnst, the amplitude would drop by 1/e.
    double m=k_avg/(w_n*w_n), stability_factor=0.2, t=0.;
    double dt=stability_factor/w_n;
    double c_crit = 2*sqrt(m*k_avg);
    double c=c_ratio*c_crit;
    C_Matrix_Dense v(totDof,1), a(totDof,1), a_n(totDof,1);

    while (t<5)
    {
        uGlobal = uGlobal + (v*dt) + (a*dt*dt*0.5);
        // ii. Dirichlet Boundary Conditions
        int bcDof=0;
        double bcValue = 0.;
        uGlobal(bcDof,0) = bcValue;
        // Calculate new acceleration
        for (int itEl = 0; itEl < mesh.num_El; itEl++) {
            calculateForceVector_AxialBar(itEl, mesh, mat, K_norm, KuGlobal, uGlobal, feL, GP_Data);
        }
        a_n = (1/m)*((fGlobal - KuGlobal) - c*v);
        // Update velocity
        v = v + (a+a_n)*(dt*0.5);
        // Update acceleration
        a = a_n;
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