#include <math.h>
#include <vector>
#include<fstream>

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
    C_Matrix_Dense uGlobal_prev(totDof,1);
    C_Matrix_Dense uGlobal_new(totDof,1);
    C_Matrix_Dense fGlobal(totDof,1);

    // ii. Dirichlet Boundary Conditions
    int bcDof=0;
    double bcValue = 0.;
    uGlobal(bcDof,0) = bcValue; 

    for (int itEl = 0; itEl < mesh.num_El; itEl++) {
        calculateForceVector_AxialBar(itEl, mesh, mat, K_norm, KuGlobal, uGlobal, feL, GP_Data);
    }

    // i. Neumann Boundary Conditions
    int fDof=totDof-1;
    fGlobal(fDof,0)=10.0;

    // Parameters for time integration scheme
    int sum = 0.0;
    for (int i = 0; i < K_norm.row_size; i++)
    {
        sum = sum + K_norm(i,0);
    }
    double k_avg = sum/K_norm.row_size;             // Request Dominic if he can introduce avg() method for Dense matrix class
    double c_ratio=1., time_cnst=0.001, w_n=1/(time_cnst*c_ratio);         // At t=time_cnst, the amplitude would drop by 1/e.
    double m=k_avg/(w_n*w_n), stability_factor=0.2, t=0., t_f=5.0;
    double dt=stability_factor/w_n;
    double c_crit = 2*sqrt(m*k_avg);
    double c=c_ratio*c_crit;
    double KE;
    double TOL=1e-3, v_norm;
    C_Matrix_Dense v_init(totDof,1), a_init(totDof,1), v(totDof,1);
    a_init = (1/m)*(-KuGlobal-c*v_init);
    uGlobal_prev = uGlobal - dt*v_init + (dt*dt/2)*a_init;

    std::ofstream file_handle;
    file_handle.open("KE.txt");

    while (t<t_f)
    {
        uGlobal_new = (1/(m/(dt*dt)+c/(2*dt)))*(fGlobal + 2/(dt*dt)*m*uGlobal_prev - KuGlobal -1/(dt*dt)*m*uGlobal_prev + 1/(2*dt)*c*uGlobal_prev);
        // ii. Dirichlet Boundary Conditions
        int bcDof=0;
        double bcValue = 0.;
        uGlobal_new(bcDof,0) = bcValue;
        // Calculate new acceleration
        KuGlobal.setZero();
        for (int itEl = 0; itEl < mesh.num_El; itEl++) {
            calculateForceVector_AxialBar(itEl, mesh, mat, K_norm, KuGlobal, uGlobal_new, feL, GP_Data);
        }

        uGlobal = uGlobal_new;
        // Calculate velocity
        v = 1/(2*dt)*(uGlobal_new-uGlobal_prev);
        KE = 0.5*m*v.inner_product(v);
        file_handle<<t<<" ,"<<KE<<"\n";
        t = t + dt;

        v_norm = v.norm();
        if (v_norm<TOL)
        {
            break;
        }
    }
    file_handle.close();


    std::cout<<"U global: ";
    std::cout<<uGlobal<<"\n";
    std::cout<<"V global: ";
    std::cout<<v<<"\n";
    std::cout<<(v_norm<TOL)<<"\n";
    
    
    return 0;
}