#include <math.h>
#include <vector>
#include<fstream>

#include "C_FEM_BasisFunction_1D.h"
#include "C_FEM_GaussPoint_1D.h" 
#include "C_Mesh_Frame.h"
#include "C_Matrix_Sparse.h"
#include "C_Matrix_Dense.h"

#include "f_ForcingVector.h"

int main()
{
    std::vector< std::vector<double> > x_mesh = {{0.,0.,0.}, {2.,0.,0.}, {1.,1.,0.}};
    std::vector< std::vector<double> > B = {{0,1,1}, {1,2,1}, {2,0,1}};
    C_Material mat;
    C_GaussPoint_1D GP_Data(2);
    C_Mesh_Frame       mesh(2);
    mesh.construct_frame(x_mesh, B);
    C_LagrangeBasis feL(1, GP_Data);
    int totDof=2*mesh.num_Nd;
    C_Matrix_Dense K_norm(mesh.num_El,1);
    C_Matrix_Dense KuGlobal(totDof,1);
    C_Matrix_Dense uGlobal(totDof,1);
    C_Matrix_Dense uGlobal_prev(totDof,1);
    C_Matrix_Dense uGlobal_new(totDof,1);
    C_Matrix_Dense fGlobal(totDof,1);

    int DOF_per_Node = 2;
    int NSPV=3;
    C_Matrix_Dense ISPV(NSPV,2);
    C_Matrix_Dense VSPV(NSPV,1);
    ISPV(0,0) = 0;
    ISPV(0,1) = 0;
    ISPV(1,0) = 0;
    ISPV(1,1) = 1;
    ISPV(2,0) = 1;
    ISPV(2,1) = 1;
    VSPV(0,0) = 0.;
    VSPV(1,0) = 0.;
    VSPV(2,0) = 0.;
    
    int NSSV=1;
    C_Matrix_Dense ISSV(NSSV,2);
    C_Matrix_Dense VSSV(NSSV,1);
    ISSV(0,0) = 2;
    ISSV(0,1) = 0;
    VSSV(0,0) = 1.0;
    
    double m=1e-2;
    double c=5e-1;
    double dt=0.002;
    double kEner;
    double t, tf=5.0;
    double TOL = 1e-3, Error;
    C_Matrix_Dense v_init(totDof,1), a_init(totDof,1), v(totDof,1);
    a_init = (1/m)*(-KuGlobal-c*v_init);
    uGlobal_prev = uGlobal - dt*v_init + (dt*dt/2)*a_init;

    std::ofstream file_handle, fHandleDisp, fHandleConn;
    file_handle.open("KE_Truss.txt");
    t=dt;
    while (t<tf)
    {
        std::cout<<t<<"\n";
        if (t<1.0)
        {
            // Applying force BCs
            for (int i = 0; i < NSSV; i++)
            {
                fGlobal(DOF_per_Node*ISSV(i,0)+ISSV(i,1),0) = t*VSSV(i,0);
            }
        }
        int bcDof=0;
        KuGlobal.setZero();
        for (int itEl = 0; itEl < mesh.num_El; itEl++) {
            calculateForceVector_Truss(itEl, mesh, mat, K_norm, KuGlobal, uGlobal, feL, GP_Data);
        }
        uGlobal_new = (1/(m/(dt*dt)+c/(2*dt)))*(fGlobal + 2/(dt*dt)*m*uGlobal - KuGlobal -1/(dt*dt)*m*uGlobal_prev + 1/(2*dt)*c*uGlobal_prev);
        // Applying displ. BCs
        for (int i = 0; i < NSPV; i++)
            {
                uGlobal_new(DOF_per_Node*ISPV(i,0)+ISPV(i,1),0) = VSPV(i,0);
            }
        // Calculate velocity
        v = 1/(2*dt)*(uGlobal_new-uGlobal_prev);
        kEner = 0.5*v.inner_product(v);
        file_handle<<t<<","<<kEner<<"\n";
        t = t + dt;
        Error = (uGlobal_new-uGlobal).norm()/uGlobal_new.norm();
        if (Error<=TOL)
        {
            std::cout<<"Soln. converged!";
            break;
        }
        
        uGlobal_prev=uGlobal;
        uGlobal = uGlobal_new;

    }
    file_handle.close();
    fHandleDisp.open("Disp_Truss.txt");
    for(int i=0; i<mesh.num_Nd; i++)
    {
        fHandleDisp<<mesh.nodes[i][0]<<","<<mesh.nodes[i][1]<<","<<mesh.nodes[i][2]<<","<<uGlobal(DOF_per_Node*i,0)<<","<<uGlobal(DOF_per_Node*i+1,0)<<"\n";
    }
    fHandleDisp.close();
    fHandleConn.open("Conn_Truss.txt");
    for(int i=0; i<mesh.num_El; i++)
    {
        fHandleConn<<mesh.elements[i][0]<<","<<mesh.elements[i][1]<<"\n";
    }
    fHandleConn.close();
    return 0;
}