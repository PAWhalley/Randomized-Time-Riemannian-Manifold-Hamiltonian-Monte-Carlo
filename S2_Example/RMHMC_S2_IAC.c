//
//  main.c
//  adaptive
//
//  Created by Peter on 8/11/2022.
//

# include <stdio.h>
# include <math.h>
# include <complex.h>
# include <stdlib.h>
# include <stdio.h>
# include <time.h>
# include "normal.h"



#define n                   3               // dimension
#define dt                  0.001           // stepsize
#define number_of_samples   100000          // final (real) time
#define burn                0
#define M                   2000            // number_of_samples // 50
#define printskip           0              // skipping for printing
#define numruns             1               // total number of trajectories

/////////////////////////////////////////////////////////////////////////////
////////////////////Basic Functions//////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

double absolute(double x){

    static double absx;

    if (x < 0){
        absx = -1.0*x;
    }
    else {
        absx = x;
    }
    return absx;
}

double min(double num1, double num2){
    double minval;
    if (num1 < num2){
        minval = num1;
    } else {
        minval = num2;
    }
    return minval;
}

double dot_product(double x[n],double v[n]){
    
    static double dot;
    static int l;

    dot = 0.;

    for (l = 0; l < n; l++){
        //Potential code
        dot += x[l]*v[l];

    }

    return dot;
}

double* matrix_vec_multiplication(double A[n][n],double x[n]){
    
    static double res[n];
    static int m,o,p;

    //setting the values to zero
    for (m = 0; m < n; m++){
        res[m] = 0.;

    }
    for (o = 0; o < n; o++){
        for (p = 0; p < n; p++){
        //Ax
        res[o] += A[o][p]*x[p];
    }

    }

    return res;
}

/////////////////////////////////////////////////////////////////////////////
////////////////////Force, Potential and Hamiltonian/////////////////////////
/////////////////////////////////////////////////////////////////////////////




void Force(double x[n],double f[n])
{
    static double* v;
    static double mean[n];
    static double A[n][n];
    static int i,j;
    //Potential Derivative
    
    mean[0] = 100.;
    for (i = 0; i < n-1; i++){
        mean[i+1] = 0.;
    }
    
    

    for (i = 0; i <n; i++){
        for (j = 0; j<n; j++){
            A[i][j] = 0.;
        }
    }
    A[0][0] = -1000.;
    A[2][2] = 1000.;

    //matrix vector multiplication
    
    v = matrix_vec_multiplication(A,x);
    

    for (i = 0; i<n; i++){
        f[i] = 2.0*v[i];
        f[i] += mean[i]; 
    } 
}


//////////////////////////////////////////////////////////////////////////////////
////////////////////////////////Constraints///////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////


double g(double x[n])
{
    static double constraint;
    static int i;

    double norm_2 = 0.;

    for (i = 0; i < n; i++){
        //Potential code
        norm_2 += x[i]*x[i];

    }
    
    constraint  = norm_2 -1.;

    return constraint;
}

void G_derv(double x[n],double constraint_derv[n])
{
    static int k;
    for (k = 0; k < n; k++){
        //Constraint Derivative code
        constraint_derv[k] = 2.0*x[k];

    }
}


double* tangent_space_gaussian(double x[n])
{
    static double v_xi[n];
    static double* v_pointer;
    static double gram_inv;
    static double proj_matrix[n][n], outer[n][n];
    static double v[n];
    static double G_derv_q[n];
    static int i,j;

    G_derv(x, G_derv_q);
    
    //outer product
    for (i = 0; i < n; i++){
        for (j = 0; j < n; j++){
        outer[i][j] = G_derv_q[i]*G_derv_q[j];
    }

    }

    //inverse of gram matrix
    gram_inv = 1.0/dot_product(G_derv_q,G_derv_q);

    for (i = 0; i<n ; i++){
        for (j = 0; j<n; j++){
            proj_matrix[i][j] = -gram_inv*outer[i][j];
        }


    }
    for (i = 0; i<n ; i++){
        proj_matrix[i][i] += 1.;

    }
    
    for (i = 0; i<n; i++){
        v[i] = r8_normal_01();

    }
    v_pointer = matrix_vec_multiplication(proj_matrix,v);//matrix mult

    for (i = 0; i<n; i++){
        v_xi[i] = v_pointer[i];

    }


    return v_xi;
}

double Potential(double x[n])
{
    static double U_res;
    static double A[n][n];
    static double mean[n];
    static double* Ax;
    static int i,j;

    U_res = 0.;

    for (i = 0; i <n; i++){
        for (j = 0; j<n; j++){
            A[i][j] = 0.;
        }
    }
    A[0][0] = -1000.;
    A[2][2] = 1000.;


    Ax = matrix_vec_multiplication(A,x);


    

    mean[0] = 100.;
    for (i = 1; i < n; i++){
        mean[i] = 0.;
    }

    U_res -= dot_product(mean,x);
    U_res -= dot_product(x,Ax);

    return U_res;
}



double Hamiltonian(double x[n],double v[n])
{
    
    static double H_res;
    static int j;

    H_res = Potential(x);
    
    for (j = 0; j < n; j++){
        //Hamiltonian code
        H_res += 0.5*v[j]*v[j];
    }
    


    return H_res;
}


//////////////////////////////////////////////////////////////////////////////////
////////////////////////////////Integrator////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

void one_RATTLE_step(double q[n], double p[n], double f[n])
{
    static double gram_inv,g_Qloop,dlambda,coeffs_v;
    static double Qloop[n], phalf[n], solve[n];
    static double G_derv_q[n];
    static double G_derv_Qloop[n];
    static int j,k;
    
    

    for(j= 0; j<n; j++)
    {
    Qloop[j] = q[j] + p[j]*dt + 0.5*dt*dt*f[j];
    }


    G_derv(q, G_derv_q);

    
    

    for(j = 0; j<50; j++)
    {
        g_Qloop = g(Qloop); 
        

        //update step
        
        G_derv(Qloop, G_derv_Qloop);

        dlambda = g_Qloop/dot_product(G_derv_Qloop,G_derv_q); //dot_product
        
        //update Q
        for(k = 0; k<n; k++){
            Qloop[k] = Qloop[k] - G_derv_q[k]*dlambda; //dot_product
        }
        if (absolute(g_Qloop)< 1e-8){
            break;
        }
    }
    
    //half step
    for (j = 0; j <n; j++){
        phalf[j] = (Qloop[j] - q[j])/dt;
        q[j] = Qloop[j];
    }
    
    //update force
    
    Force(q,f);
    
    
    //update G_q
    G_derv(q, G_derv_q);
    
    
    //inverse of gram matrix
    gram_inv = 1.0/dot_product(G_derv_q,G_derv_q);
    
    
    //linear solver for Lagrange velocity multipliers
    
    for (j = 0; j <n; j++){
        solve[j] = 2.0*phalf[j]/dt + f[j];
    }
    

    coeffs_v = gram_inv*dot_product(G_derv_q,solve);

    //full step

    for(j = 0; j<n; j++)
    {
    p[j] = phalf[j] + 0.5*dt*f[j] - 0.5*dt*coeffs_v*G_derv_q[j];
    }

}
//////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////IAC////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////


double Compute_IAC(double f_sample_list[number_of_samples])
{

    static int i,j;
    static double IAC_estimate;
    static double mean;
    static double c_f_list[M+1];
    static double c_f_current;

    mean = 0.0;

    for (i = 0; i<number_of_samples; i++){
        mean += f_sample_list[i];
    }
    mean /= number_of_samples;

    for (i = 0; i<M+1; i++){

        c_f_current = 0.0;

        for (j = 0; j< number_of_samples - i; j++){

            c_f_current +=  (f_sample_list[j] - mean)*(f_sample_list[j+i] - mean);
            
        }
        c_f_current /= (number_of_samples - i);

        c_f_list[i] = c_f_current;
    }


    IAC_estimate = 0.0;

    for (i = 1; i < M+1; i++){
        
        IAC_estimate += c_f_list[i];


    }
    
    IAC_estimate /= c_f_list[0];
    IAC_estimate *= 2.0;
    IAC_estimate += 1.0;

    return IAC_estimate;
}


//////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////MAIN///////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

void one_step_HMC(int num_leap, double qold[n], double pold[n], double q[n], double p[n], double f[n])
{
    static double u, alpha, Hnew, Hold;
    static double* random_vector;
    static int i,k;
    //
    // HMC integrator
    //
    random_vector = tangent_space_gaussian(q);
    
    for(i = 0; i<n; i++)
    {
    pold[i] = random_vector[i];
    qold[i] = q[i];
    p[i] = pold[i];
    }
    
    for (i = 0; i<num_leap; i++)
    {
    one_RATTLE_step(q,p,f);
    }
    //M-H step
    Hnew = Hamiltonian(q,p);
    Hold = Hamiltonian(qold,pold);

    u = r8_uniform_01();
    alpha = min(1.0,exp(Hold - Hnew));

    if (u > alpha)
    {
        for(k = 0; k<n; k++)
        {
        p[k] = pold[k];
        q[k] = qold[k];
        }
        Force(q,f);
    }
}

int main(int argc, const char * argv[]) {
    
    static double qold[n], pold[n], q[n], p[n], f[n];
    static int nt,counter;
    static double t;
    static double IAC_list[number_of_samples];
    static double IAC;
    static int num_leap,num;

    srand48((unsigned)time(0));

    for(nt = 0; nt<numruns; nt++)
    {
        // initialize the trajectory
        t=0;
        p[0] = 0.;
        p[1] = 1.;
        p[2] = 1.;
        q[0] = 0.;
        q[1] = 0.;
        q[2] = 1.;

        
        Force(q,f);  // force
        

    /////////////////////////////////////////////////////////////////////////////  
    ///////////////////////////////////IAC///////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////
    
    for (num_leap = 1; num_leap < 251; num_leap++){
        
            // for(nt = 0; nt<numruns; nt++)
            // {

                // initialize the trajectory
                
                p[0] = 1.;
                p[1] = 0.;
                p[2] = 0.;

                q[0] = 0.;
                q[1] = 0.;
                q[2] = 1.;

                
                Force(q,f);  // force
                
                
                

            

              
                for (num = 0; num<number_of_samples + burn; num++)
                { 
                    one_step_HMC(2*num_leap,qold,pold,q,p,f);
                    
                    //IAC computation
                    if (num > burn - 1){
                    IAC_list[num-burn] = Potential(q);                     
                    }
                }
                    
                

                
                //Compute IAC
                IAC = Compute_IAC(IAC_list);

                printf("%lf %lf ",2*num_leap*dt,IAC);
                printf("\n");
    }
            }
        

    printf("\n");

    return 0;

}


