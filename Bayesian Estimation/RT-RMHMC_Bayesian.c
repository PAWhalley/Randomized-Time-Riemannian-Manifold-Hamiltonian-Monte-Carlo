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



#define d                       30               // V_{d,p} is the manifold of
#define p                       5               // dimension d*p
#define dp                      150             //d*p
#define m                       40               //number of data points
#define n                       185              //d*p + p + d  // dimension
#define dt_max                  0.001            // stepsize
#define number_of_samples       10000          //100000// final (real) time
#define burn                    0
#define M                       200            //2000// number_of_samples // 50
#define printskip               0               // skipping for printing
#define numruns                 1               // total number of trajectories
#define num_of_constraints      15              //p*(p+1)/2       // number of constraints
#define sigma_1                 0.1
#define sigma_2                 0.1
#define M_PI                    3.14159265358979323846  
#define max_elim_iters          500

/////////////////////////////////////////////////////////////////////////////
////////////////////Basic Functions//////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////
////////////////////Cholesky Decomposition Inverse///////////////////////////
/////////////////////////////////////////////////////////////////////////////

void AT_A_multiplication(int dimension, double A[dimension][dimension], double C[dimension][dimension]){

    static int i,j,k;

    //set C to be zeroes.
    for (i = 0; i < dimension; i++){
        for (j = 0; j < dimension; j++){

            C[i][j] = 0.;

        }

    }
    for (i = 0; i < dimension; i++){
        for (j = 0; j < dimension; j++){
            for(k = 0; k < dimension; k++){
                
                C[i][j] += A[k][i]*A[k][j]; 


            }

        }

    }



}

double *cholesky(int n_chol,double A[n_chol][n_chol]) {

    //reference: W.H.Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery ,Numerical Recipes in C, The Art of Scientific Computing Second Edition


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

void chold_inv(int dimension, double a[dimension][dimension])
{

    //reference: W.H.Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery ,Numerical Recipes in C, The Art of Scientific Computing Second Edition
    
    static int i,j,k;
    static double sum;
    double L_inv[dimension][dimension];

    double *chol = cholesky(dimension,a); 

    for (i = 0; i < dimension; i++) {
        for (j = 0; j < dimension; j++)
            L_inv[i][j] = chol[i * dimension + j];
    }
    free(chol);
    //back subsitution
    for (i=0;i<dimension ;i++) {
        L_inv[i][i]=1.0/L_inv[i][i];
        for (j=i+1;j < dimension;j++) {
            sum=0.0;
            for (k=i;k<j;k++){ 
                sum -= L_inv[j][k]*L_inv[k][i];
                L_inv[j][i]=sum/L_inv[j][j];
            }
        }
     }

    //finding inverse of A through multiplication
    AT_A_multiplication(dimension,L_inv, a);


}

/////////////////////////////////////////////////////////////////////////////
///////////////////////////////Determinant///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

double Determinant(int dimension, double a[dimension][dimension])
{
    
    double *chol = cholesky(dimension,a); 
    double det = 1.;
    static int i;

    for (i = 0 ; i < dimension ; i++){

        det *= chol[i*dimension +i];
        

    }
    det = det*det;
    return(det);
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

double dot_product(int dimension, double x[dimension],double v[dimension]){
    
    static double dot;
    static int l;

    dot = 0.;

    for (l = 0; l < dimension; l++){
        //Potential code
        dot += x[l]*v[l];

    }

    return dot;
}

void matrix_vec_multiplication(int i, int j, double A[i][j],double x[j], double res[i]){
    
    static int k,l;
    //setting the values to zero
    for (k = 0; k < i; k++){
        res[k] = 0.;

    }
    for (k = 0; k < i; k++){
        for (l = 0; l < j; l++){
        //Ax
        res[k] += A[k][l]*x[l];
    }

    }

}

void vec_to_matrix(double q[dp], double X[d][p])
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

void matrix_to_vec(double X[d][p],double q[dp])
{
    static int i_index;
    static int j_index;
    static int i;

    for (i = 0; i < dp; i++)
    {
       i_index = i%d;
       j_index = (i - i_index)/d;
       q[i] = X[i_index][j_index];
    }

}


void matrix_matrix_multiplication(int x, int y, int z, double A[x][y], double B[y][z], double C[x][z]){


    //set C to be zeroes.
    static int i,j,k;
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

double sample_exp(double EX){
    static double sample;
    static double u;

    u = r8_uniform_01();

    sample = -log(1.0 - u)*EX;

    return sample;
}
/////////////////////////////////////////////////////////////////////////////
////////////////////Force, Potential and Hamiltonian/////////////////////////
/////////////////////////////////////////////////////////////////////////////




void Force(double x[n],double data[d][m],double mean[d],double f[n])
{
    static double MAT_X[d][p], MAT_X_T[p][d];
    static double x_dp[d*p];
    static int i,j,k,l,r;

    for (i = 0; i < d*p; i++){
        x_dp[i] = x[i];
        
    }
    

    vec_to_matrix(x_dp,MAT_X);


    for (i = 0; i < p ; i++){
        for (j = 0; j < d; j++){
            MAT_X_T[i][j] = MAT_X[j][i];
        }
    }

    static double d_1[p],d_2[d];
    static double d_1_mat[p][p], d_2_mat[d][d];

    for (i = 0; i < p ; i++){
        for (j = 0; j < p; j++){
            d_1_mat[i][j] = 0;

        }

    }
    for (i = 0; i < d ; i++){
        for (j = 0; j < d; j++){
            d_2_mat[i][j] = 0;

        }

    }

    for (i = 0; i<p; i++){

        d_1[i] = x[d*p+i];
        d_1_mat[i][i] = x[d*p+i];

    }

    for (i = 0; i<d; i++){

        d_2[i] = x[d*p+p+i];
        d_2_mat[i][i] = x[d*p+p+i];

    }

    //Σ
    static double X_D1[d][p];
    static double X_D1_XT[d][d];
    static double sigma_inv_T[d][d];
    static double M_kl[d][d];

    matrix_matrix_multiplication(d, p, p, MAT_X, d_1_mat, X_D1);
    matrix_matrix_multiplication(d, p, d, X_D1,MAT_X_T,X_D1_XT);
    
    for (i = 0; i < d ; i++){
        for (j = 0; j < d; j++){
            sigma_inv_T[i][j] = X_D1_XT[i][j] + d_2_mat[i][j];
        }
    }
    
    
    //Σ_inv_T
    chold_inv(d,sigma_inv_T);
    
    

    //Constructing M_kl
    for (i = 0; i < d ; i++){
        for (j = 0; j < d; j++){
            sigma_inv_T[i][j] = sigma_inv_T[j][i];
            M_kl[i][j] = 0.5*m*sigma_inv_T[i][j];
        }
    }

    
    //Data calculation
    static double Data_r[d];
    static double sigma_inv_T_k[d], sigma_inv_T_l[d];

    for (k = 0; k<d; k++){
        for (l = 0; l < d ;l++){
            for (r = 0; r < m; r++){

                //Data_r
                for ( i = 0; i<d;i++){
                Data_r[i] = data[i][r] - mean[i];
                sigma_inv_T_k[i] = sigma_inv_T[k][i];
                sigma_inv_T_l[i] = sigma_inv_T[l][i];
                }
                
                M_kl[k][l] -= 0.5*dot_product(d,Data_r,sigma_inv_T_k)*dot_product(d,Data_r,sigma_inv_T_l);
            } 

        }

    }

   
    

    //dUdX
    static double dUdX[d][p];

    for (i = 0; i<d;i++){
        for(j = 0; j <p;j++){
            dUdX[i][j] = 0.;
        }
    }

    static double dsigma_kl_dX_ij = 0.;

    for (i = 0; i<d; i++){
        for (j = 0; j<p; j++){
            for (k = 0; k <d; k++){
                for (l = 0; l < d;l ++){

                    if ((k == i) && (l == i)){
                        
                        dsigma_kl_dX_ij = 2.0*d_1[j]*MAT_X[i][j];

                    }else if(k == i){
                        
                        dsigma_kl_dX_ij = d_1[j]*MAT_X[l][j];

                    }else if(l == i){

                        dsigma_kl_dX_ij = d_1[j]*MAT_X[k][j];

                    }else {
                        continue;
                    }

                    dUdX[i][j] += M_kl[k][l]*dsigma_kl_dX_ij;

                }

            }

        }


    }


    static double dUd1[p];
    static double dsigma_kl_dD1_jj;

    for (i = 0; i< p ; i++){
        dUd1[i] = 0.;
    }
    
    for (j = 0; j <p; j++){
        for (k = 0; k < d; k++){
            for (l = 0; l < d ;l++){

                dsigma_kl_dD1_jj = MAT_X[k][j]*MAT_X[l][j];
                dUd1[j] += M_kl[k][l]*dsigma_kl_dD1_jj;

            }
        }

        dUd1[j] += d_1[j]/(pow(sigma_1,2));        

    }

    static double dUd2[d];

    for (i = 0; i< d ; i++){
        dUd2[i] = 0.;
    }

    for (j = 0; j<d ;j++){

        dUd2[j] += M_kl[j][j];

        dUd2[j] += d_2[j]/(pow(sigma_2,2));


    }

    //update force
    static double dUdX_vec[d*p];
    matrix_to_vec(dUdX,dUdX_vec);

    for (i = 0; i< d*p;i++){
        f[i] = -dUdX_vec[i];
    }
    for (j = 0; j < p;j++){
        f[d*p + j] = -dUd1[j];
    }
    for (k = 0; k<d ; k++){
        f[d*p+p + k] = -dUd2[k];
    }

}


//////////////////////////////////////////////////////////////////////////////////
////////////////////////////////Constraints///////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////


double g_ij(double x[n], int i, int j)
{
    static double constraint;
    static double MAT_X[d][p];
    static double x_dp[d*p];
    static int k;

    for (k = 0; k < d*p; k++){
        x_dp[k] = x[k];
    }



    vec_to_matrix(x_dp,MAT_X);

    if (i == j){
        constraint = 0;
        //pick ith column
        for (k = 0; k < d; k++){
            constraint += MAT_X[k][i]*MAT_X[k][i];
        }
        constraint -= 1.;

    }

    
    else {
        //pick ith column
        constraint = 0;
        for (k; k < d; k++){
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
    static double x_dp[d*p];
    static int i,i_d,i_p,j,j_d;

    for (i = 0; i < d*p; i++){
        x_dp[i] = x[i];
    }


    vec_to_matrix(x_dp,MAT_X);
    
    for (i = 0; i< num_of_constraints; i++){
        for (j = 0; j < n; j++){

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
    
    G_derv(x, G_derv_q);


    static int i,j;
    static int row, columns;

    

    //compute the transpose
    for (i = 0; i < n; i++){

        for (j = 0; j<num_of_constraints; j++){

            G_derv_q_T[i][j] = G_derv_q[j][i];
            
        }

    }

   
    //should be multiplication
    matrix_matrix_multiplication(num_of_constraints,n,num_of_constraints,G_derv_q,G_derv_q_T,gram_inv);

    

    //inverse of gram matrix.
    chold_inv(num_of_constraints,gram_inv);

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

double Potential(double x[n], double data[d][m], double mean[d])
{
    static double U_res;
    static double MAT_X[d][p], MAT_X_T[p][d];
    static double x_dp[d*p];
    static int i,j;

    for (i = 0; i < d*p; i++){
        x_dp[i] = x[i];
    }
    

    vec_to_matrix(x_dp,MAT_X);


    for (i = 0; i < p ; i++){
        for (j = 0; j < d; j++){
            MAT_X_T[i][j] = MAT_X[j][i];
        }
    }

    static double d_1[p],d_2[d];
    static double d_1_mat[p][p], d_2_mat[d][d];

    for (i = 0; i < p ; i++){
        for (j = 0; j < p; j++){
            d_1_mat[i][j] = 0;

        }

    }
    for (i = 0; i < d ; i++){
        for (j = 0; j < d; j++){
            d_2_mat[i][j] = 0;

        }

    }

    for (i = 0; i<p; i++){

        d_1[i] = x[d*p+i];
        d_1_mat[i][i] = x[d*p+i];

    }

    for (i = 0; i<d; i++){

        d_2[i] = x[d*p+p+i];
        d_2_mat[i][i] = x[d*p+p+i];

    }

    

    //Σ
    static double X_D1[d][p];
    static double X_D1_XT[d][d];
    static double sigma_inv[d][d];
    static double M_kl[d][d];

    matrix_matrix_multiplication(d, p, p, MAT_X, d_1_mat, X_D1);
    matrix_matrix_multiplication(d, p, d, X_D1,MAT_X_T,X_D1_XT);
    
    for (i = 0; i < d ; i++){
        for (j = 0; j < d; j++){
            sigma_inv[i][j] = X_D1_XT[i][j] + d_2_mat[i][j];
        }
    }


    
    //Σ_inv_T
    chold_inv(d,sigma_inv);


    

    U_res = -0.5*m*log(Determinant(d,sigma_inv)) + 0.5*m*p*log(2*M_PI);
    
    static double Data_i[d]; 
    static double Bx[d];

    

    for (i = 0; i < m; i++){
        for (j = 0; j < d;j++){
            Data_i[j] = data[j][i] - mean[j];
        }
        matrix_vec_multiplication(d,d,sigma_inv,Data_i,Bx);

        U_res += 0.5*dot_product(d,Data_i,Bx);

    }

    U_res += 0.5*log(2*M_PI*(sigma_1)*(sigma_1));
    U_res += 0.5*dot_product(p,d_1,d_1)/((sigma_1)*(sigma_1));

    U_res += 0.5*log(2*M_PI*(sigma_2)*(sigma_2));
    U_res += 0.5*dot_product(d,d_2,d_2)/((sigma_2)*(sigma_2));



    return U_res;
}



double Hamiltonian(double x[n],double v[n],double data[d][m], double mean[d])
{
    
    static double H_res;
    static int j;

    H_res = Potential(x, data, mean);
    
    for (j = 0; j < n; j++){

        //Hamiltonian code
        H_res += 0.5*v[j]*v[j];

    }
    


    return H_res;
}


//////////////////////////////////////////////////////////////////////////////////
////////////////////////////////Integrator////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

void one_RATTLE_step(double dt, double q[n], double v[n], double f[n], double data[d][m], double mean[d])
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
    static int i,j,k,Qi,G_i;
    
    G_derv(q, G_derv_q);

    
    for( j= 0; j<n; j++)
    {
    Qloop[j] = q[j] + v[j]*dt + 0.5*dt*dt*f[j];  
    }
    

    G_derv(Qloop, G_derv_Qloop);


                    
    for(k = 0; k<max_elim_iters; k++)
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
                

                dlambda = g_Qloop/dot_product(n,G_derv_Qloop_index,G_derv_index); //dot_product
                
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
    
    Force(q,data,mean,f);
    
    
    
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
    chold_inv(num_of_constraints, gram_inv);
    
    

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
    static double IAC_mean;
    static double c_f_list[M+1];
    static double c_f_current;
    static int i,j;

    IAC_mean = 0.0;

    for (i = 0; i<number_of_samples; i++){
        IAC_mean += f_sample_list[i];
    }
    IAC_mean /= number_of_samples;

    for (i = 0; i<M+1; i++){

        c_f_current = 0.0;

        for (j = 0; j< number_of_samples - i; j++){

            c_f_current +=  (f_sample_list[j] - IAC_mean)*(f_sample_list[j+i] - IAC_mean);
            
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

void one_step_RHMC(double dt, int num_leap, double qold[n], double pold[n], double q[n], double v[n], double f[n], double data[d][m], double mean[d])
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
    v[i] = pold[i];
    }
    
    for (i = 0; i<num_leap; i++)
    {
    one_RATTLE_step(dt,q,v,f,data,mean);
    }
    //M-H step
    Hnew = Hamiltonian(q,v,data,mean);
    Hold = Hamiltonian(qold,pold,data,mean);

    u = r8_uniform_01();
    alpha = min(1.0,exp(Hold - Hnew));

    if (u > alpha)
    {
        for(k = 0; k<n; k++)
        {
        v[k] = pold[k];
        q[k] = qold[k];
        }
        Force(q,data,mean,f);
    }
}

int main(int argc, const char * argv[]) {
    
    static double qold[n], pold[n], q[n], v[n], f[n];
    static double q_init[n];
    static double initial_mat[d][p];
    static int nt,counter;
    static double IAC_list[number_of_samples];
    static double IAC;
    static double* init_pointer;
    static double dt,t,L,T;
    static int i,j,num_leap,num;

    //to change the seed. Uncomment if you want fixed seed.
    srand48((unsigned)time(0));


    
    //Import Data Matrix
    FILE *fdata =  fopen("DIRECTORY_PATH/Data.txt", "ro");

    double data[d][m];


    for(i = 0; i < d; i++) {
        for(j = 0; j < m; j++) {
            if (fscanf(fdata, "%lf ", &data[i][j]) == 1) {
            }
            else {
                printf("Failed to read data.\n");
            }
        }
    } 
    
    fclose(fdata);

    FILE *fmean =  fopen("DIRECTORY_PATH/mean.txt", "ro");


    //import mean vector

    double mean[d];


    for(i = 0; i < d; i++) {
        
        if (fscanf(fmean, "%lf ", &mean[i]) == 1) {
        }
        else {
            printf("Failed to read data.\n");
        }
    } 
    
    fclose(fmean);

    //import initial vec

    FILE *fq =  fopen("DIRECTORY_PATH/initial_vec.txt", "ro");

    for(i = 0; i < n; i++) {
        
        if (fscanf(fq, "%lf ", &q_init[i]) == 1) {
        }
        else {
            printf("Failed to read data.\n");
        }
    } 
    
    fclose(fq);



    

    /////////////////////////////////////////////////////////////////////////////  
    ///////////////////////////////////IAC///////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////
    
    for (num_leap = 1; num_leap < 50; num_leap++){
        
        
        for(i = 0; i < n; i++) {
            q[i] = q_init[i];
        } 
    
        
        init_pointer = tangent_space_gaussian(q);

        for (i = 0 ; i<n ;i++){

            v[i] = init_pointer[i];
            
        }

    
        Force(q,data,mean,f);  // force
        
        T = 10*num_leap*dt_max;

        for (num = 0; num<number_of_samples + burn; num++)
        { 
            t = sample_exp(T);
            L = ceil(t/dt_max);
            dt = t/L;
            one_step_RHMC(dt,L,qold,pold,q,v,f,data,mean);
            
            //IAC computation
            if (num > burn - 1){
            IAC_list[num-burn] = Potential(q,data,mean);                     
            }
        }
        //
        
        
        
        

        //Compute IAC
        IAC = Compute_IAC(IAC_list);

        printf("%lf %lf \n",10*dt_max*num_leap,IAC);
    }
    
        
    
       
    return 0;

}

