//
//  main.c
//  
//
//

# include <stdio.h>
# include <math.h>
# include <complex.h>
# include <stdlib.h>
# include <stdio.h>
# include <time.h>
# include "normal.h"



#define d                       3              // V_{d,p} is the manifold of
#define p                       3               // dimension d*p
#define n                       9  //d*p             // dimension
#define dt                      0.01                  // stepsize
#define number_of_samples       100000             // final (real) time
#define burn                    0
#define M                       2000            // number_of_samples // 50
#define printskip               0               // skipping for printing
#define numruns                 1               // total number of trajectories
#define num_of_constraints      6 //p*(p+1)/2       // number of constraints

/////////////////////////////////////////////////////////////////////////////
////////////////////Basic Functions//////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////
////////////////////Cholesky Decomposition Inverse///////////////////////////
/////////////////////////////////////////////////////////////////////////////

void AT_A_multiplication(double A[num_of_constraints][num_of_constraints], double C[num_of_constraints][num_of_constraints]){

    static int i,j,k;

    //set C to be zeroes.
    for (i = 0; i < num_of_constraints; i++){
        for (j = 0; j < num_of_constraints; j++){

            C[i][j] = 0.;

        }

    }
    for (i = 0; i < num_of_constraints; i++){
        for (j = 0; j < num_of_constraints; j++){
            for(k = 0; k < num_of_constraints; k++){
                
                C[i][j] += A[k][i]*A[k][j]; 


            }

        }

    }



}

double *cholesky(double A[num_of_constraints][num_of_constraints], int n_chol) {
    double *L = (double*)calloc(n_chol * n_chol, sizeof(double));
    static int i,j,k;

    for (i = 0; i < n_chol; i++)
        for (j = 0; j < (i+1); j++) {
            double s = 0;
            for (k = 0; k < j; k++)
                s += L[i * n_chol + k] * L[j * n_chol + k];
            L[i * n_chol + j] = (i == j) ?
                           sqrt(A[i][i] - s) :
                           (1.0 / L[j * n_chol + j] * (A[i][j] - s));
        }

    return L;
}

void chold_inv(double a[num_of_constraints][num_of_constraints])
{
    
    static int i,j,k;
    static double sum;
    static double L_inv[num_of_constraints][num_of_constraints];

    double *chol = cholesky(a,num_of_constraints); 

    for (i = 0; i < num_of_constraints; i++) {
        for (j = 0; j < num_of_constraints; j++)
            L_inv[i][j] = chol[i * num_of_constraints + j];
    }
    free(chol);
    //back subsitution
    for (i=0;i<num_of_constraints ;i++) {
        L_inv[i][i]=1.0/L_inv[i][i];
        for (j=i+1;j < num_of_constraints;j++) {
            sum=0.0;
            for (k=i;k<j;k++){ 
                sum -= L_inv[j][k]*L_inv[k][i];
                L_inv[j][i]=sum/L_inv[j][j];
            }
        }
     }

    //finding inverse of A through multiplication
    AT_A_multiplication(L_inv, a);


}

/////////////////////////////////////////////////////////////////////////

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

    for (int l = 0; l < n; l++){
        //Potential code
        dot += x[l]*v[l];

    }

    return dot;
}

void matrix_vec_multiplication(int i, int j, double A[i][j],double x[j], double res[i]){
    static int m,k,l;    

    //setting the values to zero
    for (m = 0; m < i; m++){
        res[m] = 0.;

    }
    for (k = 0; k < i; k++){
        for (l = 0; l < j; l++){
        //Ax
        res[k] += A[k][l]*x[l];
    }

    }

}

void vec_to_matrix(double q[n], double X[d][p])
{
    static int i,j;

    for (i = 0; i < d; i++)
    {
        for (j = 0; j < p; j++){
            X[i][j] = q[j*d+i];

        }

    }

}

void vec_to_matrix_square(double q[num_of_constraints*num_of_constraints], double X[num_of_constraints][num_of_constraints])
{
    static int i,j;
    for (i = 0; i < num_of_constraints; i++)
    {
        for (j = 0; j < num_of_constraints; j++){
            X[i][j] = q[j*num_of_constraints+i];

        }

    }

}

void matrix_to_vec_square(double X[num_of_constraints][num_of_constraints],double q[num_of_constraints*num_of_constraints])
{
    static int i_index;
    static int j_index;
    static int i;

    for (i = 0; i < num_of_constraints*num_of_constraints; i++)
    {
       i_index = i%num_of_constraints;
       j_index = (i - i_index)/num_of_constraints;
       q[i] = X[i_index][j_index];
    }

}

void matrix_to_vec(double X[d][p],double q[n])
{
    static int i_index;
    static int j_index;
    static int i;

    for (i = 0; i < n; i++)
    {
       i_index = i%d;
       j_index = (i - i_index)/d;
       q[i] = X[i_index][j_index];
    }

}


void matrix_matrix_multiplication(int x, int y, int z, double A[x][y], double B[y][z], double C[x][z]){

    static int i,j,k;
    //set C to be zeroes.
    for (i = 0; i < x; i++){
        for (j = 0; j < z; j++){

            C[i][j] = 0.;

        }

    }
    for (i = 0; i < x; i++){
        for (j = 0; j < z; j++){
            for(k = 0; k < y ; k++){
                
                C[i][j] += A[i][k]*B[k][j]; 


            }

        }

    }



}

/////////////////////////////////////////////////////////////////////////////
////////////////////Force, Potential and Hamiltonian/////////////////////////
/////////////////////////////////////////////////////////////////////////////




void Force(double x[n],double f[n])
{
    static double K;
    static int i;
    K = 1.;

    static double Fmat[d][p];

    //define F

    Fmat[0][0] = 0.;
    Fmat[1][1] = 0.;
    Fmat[2][2] = 0.;
    Fmat[0][1] = 2.;
    Fmat[0][2] = -45.;
    Fmat[1][0] = -2.;
    Fmat[1][2] = -4.;
    Fmat[2][0] = 45.;
    Fmat[2][1] = 4.;


    matrix_to_vec(Fmat,f);

    for (i = 0; i < n; i++){

        f[i] *= K;

    }

}


//////////////////////////////////////////////////////////////////////////////////
////////////////////////////////Constraints///////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////


double g_ij(double x[n], int i, int j)
{
    static double constraint;
    static double MAT_X[d][p];

    vec_to_matrix(x,MAT_X);

    if (i == j){
        constraint = 0;
        //pick ith column
        for (int k; k < d; k++){
            constraint += MAT_X[k][i]*MAT_X[k][i];
        }
        constraint -= 1.;

    }

    
    else {
        //pick ith column
        constraint = 0;
        for (int k; k < d; k++){
            constraint += MAT_X[k][i]*MAT_X[k][j];
        }
    }
    
    return constraint;
}

/////////////////////////////////////

void G_derv(double x[n],double constraint_derv[num_of_constraints][n])
{   
    //Constraint Derivative code
    //requires i<=j otherwise not invertible!
    static double MAT_X[d][p];
    static int i,j,i_d,i_p,j_d;

    vec_to_matrix(x,MAT_X);
    
    for (int i = 0; i< num_of_constraints; i++){
        for (int j = 0; j < n; j++){

            constraint_derv[i][j] = 0.;

        }

    }


    for (i = 0; i < p; i++){
        //block diagonals
        for (i_d = 0; i_d < d; i_d++){
            for (i_p = 0; i_p < p-i; i_p ++){
                constraint_derv[p*i - i*(i-1)/2 + i_p][d*i + i_d] = MAT_X[i_d][i + i_p];
            }
        } 
        //vector diagonals
        for (j = 0; j < p- i; j++){
            for (j_d = 0; j_d < d ; j_d ++){

                constraint_derv[p*i - i*(i-1)/2 + j][(j + i)*d + j_d] += MAT_X[j_d][i];

            }   
        }
    }
}


double* tangent_space_gaussian(double x[n])
{
    static double v_xi[n];
    static double gram[n][n];
    static double gram_inv[num_of_constraints][num_of_constraints];
    static double gram_inv_vec[num_of_constraints*num_of_constraints];
    static double proj_matrix[n][n];
    static double v[n];
    static double G_derv_q[num_of_constraints][n];
    static double G_derv_q_T[n][num_of_constraints];
    static double gram_inv_G_q[num_of_constraints][n];
    static int i,j;

    G_derv(x, G_derv_q);


    int row, columns;

    

    //compute the transpose
    for (i = 0; i < n; i++){

        for (j = 0; j<num_of_constraints; j++){

            G_derv_q_T[i][j] = G_derv_q[j][i];

        }

    }

   
    //should be multiplication
    matrix_matrix_multiplication(num_of_constraints,n,num_of_constraints,G_derv_q,G_derv_q_T,gram_inv);

    

    //inverse of gram matrix.
    chold_inv(gram_inv);

    matrix_matrix_multiplication(num_of_constraints,num_of_constraints,n,gram_inv,G_derv_q,gram_inv_G_q);

    matrix_matrix_multiplication(n,num_of_constraints,n,G_derv_q_T,gram_inv_G_q,proj_matrix);


    for (i = 0; i<n; i++){
        for(j = 0; j<n; j++){
            proj_matrix[i][j] = -proj_matrix[i][j];
       }
    }
    for (i = 0; i<n ; i++){
        //adding identity.
        proj_matrix[i][i] += 1.;

    }
    
    for (i = 0; i<n; i++){
        v[i] = r8_normal_01();

    }
    matrix_vec_multiplication(n,n,proj_matrix,v,v_xi);//matrix mult


    return v_xi;
}

double Potential(double x[n])
{
    static double U_res;
    static double K;
    K = 1.;

    static double Fmat[d][p];
    static double Fvec[n];
    static int i;

    //define F

    Fmat[0][0] = 0.;
    Fmat[1][1] = 0.;
    Fmat[2][2] = 0.;
    Fmat[0][1] = 2.;
    Fmat[0][2] = -45.;
    Fmat[1][0] = -2.;
    Fmat[1][2] = -4.;
    Fmat[2][0] = 45.;
    Fmat[2][1] = 4.;

    matrix_to_vec(Fmat,Fvec);

    for (i = 0; i < n; i++){

        Fvec[i] *= -K;

    }

    U_res = dot_product(Fvec,x);

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

void one_RATTLE_step(double q[n], double v[n], double f[n])
{
    static double gram_inv[num_of_constraints][num_of_constraints],g_Qloop,dlambda,coeffs_v[num_of_constraints];
    static double gram_inv_vec[num_of_constraints*num_of_constraints];
    static double Qloop[n], phalf[n], solve[n];
    static double G_derv_q[num_of_constraints][n];
    static double G_derv_q_T[n][num_of_constraints];
    static double G_derv_Qloop[num_of_constraints][n];
    static double G_derv_index[n], G_derv_Qloop_index[n];
    static int    index, break_counter;
    static double residual_list[num_of_constraints];
    static double b[num_of_constraints];
    static int i,j,k, G_i,Qi;
    
    G_derv(q, G_derv_q);

    
    for(j= 0; j<n; j++)
    {
    Qloop[j] = q[j] + v[j]*dt + 0.5*dt*dt*f[j];
    
    }

    G_derv(Qloop, G_derv_Qloop);
                    
    for(k = 0; k<50; k++)
    {   
        break_counter = 0;

        for (i = 0; i < p ; i++){
            for (j = i; j < p; j++){

                g_Qloop = g_ij(Qloop,i,j); 
                

                index = i*p - i*(i-1)/2 + j - i;

                residual_list[index] = g_Qloop;

                if (absolute(g_Qloop)< 1e-8){
                    continue;
                }
                
                //update step
                G_derv(Qloop, G_derv_Qloop);

                
                for (G_i = 0; G_i < n; G_i++){
                    
                    G_derv_index[G_i] = G_derv_q[index][G_i];
                   
                }

                for (G_i = 0; G_i < n; G_i++){
                    
                    G_derv_Qloop_index[G_i] = G_derv_Qloop[index][G_i];
                   
                }
                

                dlambda = g_Qloop/dot_product(G_derv_Qloop_index,G_derv_index); //dot_product

                //update Q
                for(Qi = 0; Qi<n; Qi++){
                    Qloop[Qi] = Qloop[Qi] - G_derv_index[Qi]*dlambda; //dot_product
                }
            }
        }
        

        for (i = 0; i < num_of_constraints; i++){
            if (absolute(residual_list[i])<1e-8){
                break_counter += 1;
            }
        }
        if (break_counter == num_of_constraints){
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

    

    //compute the transpose
    for (i = 0; i < n; i++){

        for (j = 0; j<num_of_constraints; j++){

            G_derv_q_T[i][j] = G_derv_q[j][i];


        }

    }
    
    //should be multiplication
    matrix_matrix_multiplication(num_of_constraints,n,num_of_constraints,G_derv_q,G_derv_q_T,gram_inv);

    //inverse of gram matrix
    chold_inv(gram_inv);
    
    
    
    //linear solver for Lagrange velocity multipliers
    
    for (j = 0; j <n; j++){
        solve[j] = 2.0*phalf[j]/dt + f[j];
    }
    //find b
    matrix_vec_multiplication(num_of_constraints,n,G_derv_q,solve,b);

    //find coeff_v
    matrix_vec_multiplication(num_of_constraints,num_of_constraints,gram_inv,b,coeffs_v);

    //////////////
    //full step
    matrix_vec_multiplication(n,num_of_constraints,G_derv_q_T,coeffs_v,solve);


    for(j = 0; j<n; j++)
    {
    v[j] = phalf[j] + 0.5*dt*f[j] - 0.5*dt*solve[j];
    }

}
//////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////IAC////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////


double Compute_IAC(double f_sample_list[number_of_samples])
{


    static double IAC_estimate;
    static double mean;
    static double c_f_list[M+1];
    static double c_f_current;
    static int i,j;

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

void one_step_HMC(int num_leap, double qold[n], double pold[n], double q[n], double v[n], double f[n])
{
    static double u, alpha, Hnew, Hold;
    static double* random_vector;
    static int i,j,k;
    //
    // HMC integrator
    //
    random_vector = tangent_space_gaussian(q);

    for(i = 0; i<n; i++)
    {
    pold[i] = random_vector[i];
    qold[i] = q[i];
    v[i] = pold[i];
    }
    
    for (i = 0; i<num_leap; i++)
    {
    one_RATTLE_step(q,v,f);
    }
    //M-H step
    Hnew = Hamiltonian(q,v);
    Hold = Hamiltonian(qold,pold);

    u = r8_uniform_01();
    alpha = min(1.0,exp(Hold - Hnew));

    if (u > alpha)
    {
        for(k = 0; k<n; k++)
        {
        v[k] = pold[k];
        q[k] = qold[k];
        }
        Force(q,f);
    }
}

int main(int argc, const char * argv[]) {
    
    static double qold[n], pold[n], q[n], v[n], f[n];
    static double initial_mat[d][p];
    static int nt,counter,num_leap,i,j,num;
    static double IAC_list[number_of_samples];
    static double IAC;
    static double* init_pointer;

    //to change the seed. Uncomment if you want fixed seed.
    srand48((unsigned)time(0));
    


    /////////////////////////////////////////////////////////////////////////////  
    ///////////////////////////////////IAC///////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////
    
    for (num_leap = 1; num_leap < 101; num_leap++){
        
            // for(nt = 0; nt<numruns; nt++)
            // {

               // initialize the trajectory
        for (i = 0; i < d; i++){
            for (j = 0; j < p; j++){
                if (i == j){
                    initial_mat[i][j] = 1.;
                }
                else {
                    initial_mat[i][j] = 0.;

                }
            }

        }
        matrix_to_vec(initial_mat,q);

        init_pointer = tangent_space_gaussian(q);

        for (i = 0 ; i<n ;i++){

            v[i] = init_pointer[i];

        }
    
        Force(q,f);  // force
              
        for (num = 0; num<number_of_samples + burn; num++)
        { 
            one_step_HMC(num_leap,qold,pold,q,v,f);
            
            //IAC computation
            if (num > burn - 1){
            IAC_list[num-burn] = Potential(q);                     
            }
        }
            
        //Compute IAC
        IAC = Compute_IAC(IAC_list);

        printf("%lf %lf ",num_leap*dt,IAC);
        printf("\n");
    }
    
        

    printf("\n");

    

    
    
       
    return 0;

}


