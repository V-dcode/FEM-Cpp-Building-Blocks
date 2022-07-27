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
    C_Mesh_1D       mesh(10.0, 10);
    C_LagrangeBasis feL(1, GP_Data);
    int totDof=mesh.num_Nd;
    C_Matrix_Dense K_norm(mesh.num_El,1);
    C_Matrix_Dense KuGlobal(totDof,1);
    C_Matrix_Dense uGlobal(totDof,1);
    C_Matrix_Dense uGlobal_prev(totDof,1);
    C_Matrix_Dense uGlobal_new(totDof,1);
    C_Matrix_Dense fGlobal(totDof,1);

    int bcDof=0;
    double bcValue = 0.;
    uGlobal(bcDof,0) = bcValue; 
    int fDof=totDof-1;
    double m=1e-7;
    double c=0.05;
    double dt=0.002;
    double kEner;
    double t, tf=5.0;
    C_Matrix_Dense v_init(totDof,1), a_init(totDof,1), v(totDof,1);
    a_init = (1/m)*(-KuGlobal-c*v_init);
    uGlobal_prev = uGlobal - dt*v_init + (dt*dt/2)*a_init;

    std::ofstream file_handle, fHandleDisp;
    file_handle.open("KE.txt");
    t=dt;
    while (t<tf)
    {
        std::cout<<t<<"\n";
        if (t<1.0)
        {
            fGlobal(fDof,0)=t;
        }
        int bcDof=0;
        KuGlobal.setZero();
        for (int itEl = 0; itEl < mesh.num_El; itEl++) {
            calculateForceVector_AxialBar(itEl, mesh, mat, K_norm, KuGlobal, uGlobal, feL, GP_Data);
        }
        uGlobal_new = (1/(m/(dt*dt)+c/(2*dt)))*(fGlobal + 2/(dt*dt)*m*uGlobal_prev - KuGlobal -1/(dt*dt)*m*uGlobal_prev + 1/(2*dt)*c*uGlobal_prev);
        double bcValue = 0.;
        uGlobal_new(bcDof,0) = bcValue;
        uGlobal = uGlobal_new;
        // Calculate velocity
        v = 1/(2*dt)*(uGlobal-uGlobal_prev);
        uGlobal_prev=uGlobal;
        kEner = 0.5*v.inner_product(v);
        file_handle<<t<<","<<kEner<<"\n";
        t = t + dt;
    }
    file_handle.close();
    fHandleDisp.open("disp.txt");
    for(int i=0; i<mesh.num_Nd; i++)
    {
        fHandleDisp<<mesh.nodes[i]<<","<<uGlobal(i,0)<<"\n";
    }
    fHandleDisp.close();
    return 0;
}