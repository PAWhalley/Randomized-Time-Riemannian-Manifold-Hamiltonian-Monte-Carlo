/***
Code for covariance estimation using the NERCOME estimator
v1.0 24 Nov, 2016
Author: Benjamin Joachimi (UCL)

Licence: 

see README for details
***/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf_gamma.h>

#define NLINE 8           // number of parameters to be read from parameter file 
#define MAX_CHAR_FILE 100000  // max. number of characters per line in files
#define COMB_MIN 3         // min. number of samples to average for nercome averaging
#define COMB_THRESHOLD 3   // factor to stay below max. number of available combinations for nercome averaging
#define AVOID_FRAC 0.1     // fraction at lower and upper end to avoid for testing of optimal split location for nercome estimator
#define GENERATE_TYPE 1    // 1: id; 2: scaled id; 3: id with off-diagonals

void read_parameter_file();
void generate_input_data();
void read_matrix();
void write_matrix();
void write_eigen();
void write_split();
void subtract_mean();
void rescale_input();
void rescale_matrix();
void cov_maxlike();
void cov_shrink_nercome();
void nercome_average();
FILE *bj_fileopen();
int bj_get_file_columns();
int bj_get_file_rows();
double **bj_alloc_2D();
void bj_free_2D();
int inverse_symm();
void gaussian_correlated_withrnd();

int *sequence;
gsl_rng *rnd;

typedef struct {
  int nd,nr,np,split_steps,numb_estimates;
  char path[200],outident[200],inputfile[200],normfile[200];
} cov_est_pars;
  


int main(int argc, char *argv[])
{
  int i,flag=0,split_opt;
  cov_est_pars *parameters=malloc(sizeof(cov_est_pars));
  char filename[200],parfile[200];


  // process command line & read parameter file
  if (argc!=2) {
    printf("Syntax: %s <parameter file name (full path)>\n\n",argv[0]);
    exit(-1);
  }
  sprintf(parfile,"%s",argv[1]);
  
  read_parameter_file(parfile,parameters);


  // memory allocation
  rnd=gsl_rng_alloc(gsl_rng_default);
  gsl_rng_env_setup();  // reads seed from GSL_RNG_SEED environment variable; default 0
    
  gsl_matrix *data_matrix=gsl_matrix_alloc(parameters->nd,parameters->nr);
  gsl_matrix *cov=gsl_matrix_alloc(parameters->nd,parameters->nd);
  gsl_matrix *inv=gsl_matrix_alloc(parameters->nd,parameters->nd);
  gsl_matrix *cov_sample=gsl_matrix_alloc(parameters->nd,parameters->nd);
  gsl_matrix *inv_sample=gsl_matrix_alloc(parameters->nd,parameters->nd);
  gsl_vector *norm_data_vector=gsl_vector_alloc(parameters->nd);


  // input
  if (strcmp(parameters->normfile,"-")!=0) {
    gsl_matrix_view tmp;
    sprintf(filename,"%s/%s",parameters->path,parameters->normfile);
    tmp=gsl_matrix_view_vector(norm_data_vector,parameters->nd,1);
    read_matrix(filename,&tmp.matrix);  // writes input into norm_data_vector
  }
  else {
    gsl_vector_set_all(norm_data_vector,1.);
  }

  if (strcmp(parameters->inputfile,"-")==0) {
    generate_input_data(data_matrix);
  }
  else {
    sprintf(filename,"%s/%s",parameters->path,parameters->inputfile);
    read_matrix(filename,data_matrix);
  }
  subtract_mean(data_matrix);  // ensures data vector now has mean zero
  rescale_input(data_matrix,norm_data_vector);  // normalises data vector by analytic normalisation

  sequence=calloc(parameters->nr,sizeof(int));  // set global variable with integer sequence
  for (i=0;i<parameters->nr;i++) {
    sequence[i]=i;
  }


  // sample covariance estimate 
  cov_maxlike(data_matrix,cov_sample);


  // calculate inverse covariance estimates
  flag=inverse_symm(cov_sample,inv_sample);
  if (flag) {
    printf("Warning: sample covariance estimate not invertible for nr=%i; nd=%i\n",parameters->nr,parameters->nd);
  }
  else {
    rescale_matrix(inv_sample,norm_data_vector,1); // un-does normalisation rescaling of data matrix
    sprintf(filename,"%s/cov_shrink_sample_inverse_%s.dat",parameters->path,parameters->outident);
    write_matrix(filename,inv_sample);
  }

  rescale_matrix(cov_sample,norm_data_vector,0); // un-does normalisation rescaling of data matrix
  sprintf(filename,"%s/cov_shrink_sample_covariance_%s.dat",parameters->path,parameters->outident);
  write_matrix(filename,cov_sample);


  // write eigenvalues of covariance
  sprintf(filename,"%s/cov_shrink_sample_eigenvalues_%s.dat",parameters->path,parameters->outident);
  write_eigen(filename,cov_sample,parameters->nr);


  // NERCOME estimate
  cov_shrink_nercome(data_matrix,parameters,cov,&split_opt);
  sprintf(filename,"%s/cov_shrink_nercome_split_%s.dat",parameters->path,parameters->outident);
  write_split(filename,split_opt,parameters->nr);


  // calculate inverse covariance estimates
  flag=inverse_symm(cov,inv);
  if (flag) {
    printf("Warning: NERCOME covariance estimate not invertible for nr=%i; nd=%i\n",parameters->nr,parameters->nd);
  }
  else {
    rescale_matrix(inv,norm_data_vector,1); // un-does normalisation rescaling of data matrix
    sprintf(filename,"%s/cov_shrink_nercome_inverse_%s.dat",parameters->path,parameters->outident);
    write_matrix(filename,inv);
  }

  rescale_matrix(cov,norm_data_vector,0); // un-does normalisation rescaling of data matrix
  sprintf(filename,"%s/cov_shrink_nercome_covariance_%s.dat",parameters->path,parameters->outident);
  write_matrix(filename,cov);


  // write eigenvalues of covariance
  sprintf(filename,"%s/cov_shrink_nercome_eigenvalues_%s.dat",parameters->path,parameters->outident);
  write_eigen(filename,cov,parameters->nr);


  // clean up  
  gsl_matrix_free(data_matrix);
  gsl_matrix_free(cov);
  gsl_matrix_free(inv);
  gsl_matrix_free(cov_sample);
  gsl_matrix_free(inv_sample);
  gsl_vector_free(norm_data_vector);
  gsl_rng_free(rnd);
  free(sequence);
  free(parameters);
  return 0;
}




// read parameter file
void read_parameter_file(char filename[],cov_est_pars *par)
{
  int i=0;
  char *word,buffer[MAX_CHAR_FILE];
  FILE *dat;

  dat=bj_fileopen("r",filename);

  while((!feof(dat))&&(i<NLINE)) {
    fgets(buffer,MAX_CHAR_FILE,dat);
    word=strtok(buffer," ");     // separates buffer into words
    word[strcspn(word,"\n")]=0;  // remove trailing new line character
    if ((strcmp(word,"#")==0)||(strcmp(word,"\0")==0)) continue;  // ignore comment lines and empty lines
    else {
      switch(i) {
      case 0:
	sprintf(par->path,"%s",word);
	break;
      case 1:
	sprintf(par->inputfile,"%s",word);
	break;
      case 2:
	sprintf(par->outident,"%s",word);
	break;
      case 3:
	par->nd=atoi(word);
	break;
      case 4:
	par->nr=atoi(word);
	break;
      case 5:
	par->split_steps=atoi(word);
	break;
      case 6:
	par->numb_estimates=atoi(word);
	break;
      case 7:
	sprintf(par->normfile,"%s",word);
	break;
      default:
	printf("Error: failed to read from parameter file.\n");
	exit(-1);
      }  // end of switch
    }
    i++;
  };

  fclose(dat);
  return;
}


// generate artificial data vector for testing purposes
void generate_input_data(gsl_matrix *y)
{
  int i,j;

  for (i=0;i<y->size1;i++) {
    for (j=0;j<y->size2;j++) {
      if (GENERATE_TYPE==1) gsl_matrix_set(y,i,j,gsl_ran_gaussian(rnd,1.));
      if (GENERATE_TYPE==2) {
	if (i>y->size1/2-1) gsl_matrix_set(y,i,j,gsl_ran_gaussian(rnd,5.));  // rescaled variance
	else gsl_matrix_set(y,i,j,gsl_ran_gaussian(rnd,1.));
	//gsl_matrix_set(y,i,j,gsl_ran_gaussian(rnd,1.+0.2*i));
      }
    }
  }
  if (GENERATE_TYPE==3) {
    double corr=0.6;
    double mean[y->size1];
    double **res=bj_alloc_2D(y->size1,y->size2);
    gsl_matrix *cov=gsl_matrix_alloc(y->size1,y->size1);

    for (i=0;i<y->size1;i++) {
      mean[i]=0.0;
      for (j=0;j<y->size1;j++) {
	gsl_matrix_set(cov,i,j,pow(corr,fabs(1.*(i-j))));
      }
    }

    gaussian_correlated_withrnd(res,y->size1,y->size2,mean,cov,rnd);  // provide random sampler to avoid repetition on separate cals

    for (i=0;i<y->size1;i++) {
      for (j=0;j<y->size2;j++) {
	gsl_matrix_set(y,i,j,res[i][j]);
      }
    }

    gsl_matrix_free(cov);
    bj_free_2D(res,y->size1,y->size2);
  }
  return;
}


// eigen-decomposition of symmetric matrix
void eigen_decomp(gsl_matrix *m,gsl_vector *eval,gsl_matrix *evec)
{
  const int dim=(int)m->size1;
  gsl_matrix *u=gsl_matrix_alloc(dim,dim);
  gsl_vector *work=gsl_vector_alloc(dim);

  gsl_matrix_memcpy(u,m);   // ensure m remains unmodified
  gsl_linalg_SV_decomp(u,evec,eval,work);  // u should be same as evec

  gsl_matrix_free(u);
  gsl_vector_free(work);
  return;
}


// read matrix from ascii file written in matrix form
void read_matrix(char filename[], gsl_matrix *y)
{
  int i,j;
  int ndat,offset,read;
  double val;
  char *buffer=malloc(MAX_CHAR_FILE*sizeof(char));
  FILE *dat;

  ndat=bj_get_file_rows(filename);
  if (ndat!=(int)y->size1) {
    printf("Error in routine 'read_matrix': incorrect number of rows in file %s\n",filename);
    exit(-1);
  }
  ndat=bj_get_file_columns(filename);
  if (ndat<(int)y->size2) {
    printf("Error in routine 'read_matrix': too few columns in file %s - requested %i but found %i\n",filename,(int)y->size2,ndat);
    exit(-1);
  }

  dat=bj_fileopen("r",filename);
  for (i=0;i<y->size1;i++) {
    offset=0;
    fgets(buffer,MAX_CHAR_FILE,dat);
    for (j=0;j<y->size2;j++) {
      sscanf(buffer+offset,"%lf%n",&val,&read); 
      offset+=read;
      gsl_matrix_set(y,i,j,val);
    }
  }
  fclose(dat);
  free(buffer);
  return;
}


// write matrix to file
void write_matrix(char filename[], gsl_matrix *m)
{
  int i,j;
  FILE *dat;
  
  dat=bj_fileopen("w",filename);
  for (i=0;i<m->size1;i++) {
    for (j=0;j<m->size2;j++) {
      fprintf(dat,"%15.10e\t",gsl_matrix_get(m,i,j));
    }
    fprintf(dat,"\n");
  }
  fclose(dat);
  return;
}


// write eigenvalues of cov to file
void write_eigen(char filename[],gsl_matrix *m,int nr)
{
  int i;
  const int nd=(int)m->size1;
  gsl_matrix *eigen_vec=gsl_matrix_alloc(nd,nd);
  gsl_vector *eigen_val=gsl_vector_alloc(nd);
  FILE *dat;

  eigen_decomp(m,eigen_val,eigen_vec);

  dat=bj_fileopen("w",filename);
  fprintf(dat,"%i\t",nr);
  for (i=0;i<nd;i++) {
    fprintf(dat,"%15.10g\t",sqrt(gsl_vector_get(eigen_val,i)));
  }
  fprintf(dat,"\n");
  
  fclose(dat);
  gsl_matrix_free(eigen_vec);
  gsl_vector_free(eigen_val);
  return;
}


// write ptimal split value to file
void write_split(char filename[],int split,int nr)
{
  FILE *dat;
  dat=bj_fileopen("w",filename);
  fprintf(dat,"%i\t%i\n",nr,split);
  fclose(dat);
  return;
}


// subtract mean from every column in data matrix
void subtract_mean(gsl_matrix *y)
{
  int i;
  gsl_vector_view column;
  gsl_vector *c=gsl_vector_alloc(y->size1);
  gsl_vector *mean=gsl_vector_alloc(y->size1);
  gsl_vector_set_zero(mean);

  for (i=0;i<y->size2;i++) {
    column=gsl_matrix_column(y,i);          // extracts ith matrix column
    gsl_vector_add(mean,&column.vector);    // adds it to previous columns
  }
  gsl_vector_scale(mean,1./(1.*y->size2));  // now contains vector of the means

  for (i=0;i<y->size2;i++) {
    column=gsl_matrix_column(y,i);          // extracts ith matrix column
    gsl_vector_sub(&column.vector,mean);    // replaces every column of y with mean-subtracted version
  }

  gsl_vector_free(c);
  gsl_vector_free(mean);
  return;
}


// rescales input data matrix with analytic model to improve stability of covarianc inversion
void rescale_input(gsl_matrix *y,gsl_vector *norm)
{
  int i;
  gsl_vector_view column;

  for (i=0;i<y->size2;i++) {
    column=gsl_matrix_column(y,i);        // extracts column i of y
    gsl_vector_div(&column.vector,norm);  // replaces column i of y with (column i)/norm
  }
  return;
}


// un-does rescaling of 'rescale_input' for covariance matrices
// ops=0 for covariance; ops=1 for inverse and rescaling of target
void rescale_matrix(gsl_matrix *inv,gsl_vector *norm,int ops)
{
  int i,j;
  double val;

  for (i=0;i<inv->size1;i++) {
    for (j=0;j<inv->size2;j++) {
      val=gsl_matrix_get(inv,i,j);
      if (!ops) gsl_matrix_set(inv,i,j,val*(gsl_vector_get(norm,i)*gsl_vector_get(norm,j)));
      else gsl_matrix_set(inv,i,j,val/(gsl_vector_get(norm,i)*gsl_vector_get(norm,j)));
    }
  }
  return;
}


// apply standard max. likelihood covariance estimator (assume mean is unknown)
void cov_maxlike(gsl_matrix *y,gsl_matrix *cov)
{
  if (y->size2<=1) {
    printf("Error in routine 'cov_maxlike': cannot compute covariance with less than 2 realisations - nr=%zu\n",y->size2);
    exit(-1);
  }
  gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,y,y,0.0,cov);  // product y*y^t
  gsl_matrix_scale(cov,1./(y->size2-1.));
  return;
}


// computes the square of the Frobenius norm of a general matrix
double frobenius_norm_squared(gsl_matrix *m)
{
  int i;
  double sum;
  gsl_matrix *res=gsl_matrix_alloc(m->size2,m->size2);
  
  gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,m,m,0.0,res);

  sum=0.0;
  for (i=0;i<m->size2;i++) {
    sum+=gsl_matrix_get(res,i,i);  // take trace
  }
  
  gsl_matrix_free(res);
  return(sum);
}


// basic nercome estimator based on splitting realisations into 2 subsets
void nercome_split_estimator(gsl_matrix *y1,gsl_matrix *y2,gsl_matrix *res,gsl_matrix *sample_cov2)
{
  int i;
  const int nd=(int)y1->size1;
  gsl_matrix *sample_cov1=gsl_matrix_alloc(nd,nd);
  gsl_matrix *eigen_vec=gsl_matrix_alloc(nd,nd);
  gsl_matrix *tmp=gsl_matrix_alloc(nd,nd);
  gsl_matrix *tmp2=gsl_matrix_alloc(nd,nd);
  gsl_vector *eigen_val=gsl_vector_alloc(nd);
  
  cov_maxlike(y1,sample_cov1);
  cov_maxlike(y2,sample_cov2);  // this is returned as well for split location optimisation
  eigen_decomp(sample_cov1,eigen_val,eigen_vec);  // eigenvalues unused

  gsl_blas_dsymm(CblasLeft,CblasUpper,1.0,sample_cov2,eigen_vec,0.0,tmp);  // tmp=Sigma2*P1
  gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,eigen_vec,tmp,0.0,tmp2); // tmp2=P1^t*Sigma2*P1
 
  gsl_matrix_set_zero(tmp);
  for (i=0;i<nd;i++) {
    gsl_matrix_set(tmp,i,i,gsl_matrix_get(tmp2,i,i));  // keep only diagonal of tmp2 in tmp
  }

  gsl_blas_dsymm(CblasRight,CblasUpper,1.0,tmp,eigen_vec,0.0,tmp2);  // tmp2=P1*diag
  gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,tmp2,eigen_vec,0.0,res); // res=P1*diag*P1^t

  gsl_matrix_free(sample_cov1);
  gsl_matrix_free(eigen_vec);
  gsl_matrix_free(tmp);
  gsl_matrix_free(tmp2);  
  gsl_vector_free(eigen_val);
  return;
}


// average over permutations of realisations in subsets for nercome split estimator
void nercome_average(gsl_matrix *y,int split,int numb_estimates,gsl_matrix *average,gsl_matrix *sigma2_average)
{
  int i,j,index1,index2,n_est=0;
  int *comb=calloc(split,sizeof(int));
  gsl_vector_view column;
  gsl_matrix_view submatrix;
  gsl_matrix *y1=gsl_matrix_alloc(y->size1,split);
  gsl_matrix *y2=gsl_matrix_alloc(y->size1,y->size2-split);
  gsl_matrix *cov=gsl_matrix_alloc(y->size1,y->size1);
  gsl_matrix *sigma2=gsl_matrix_alloc(y->size1,y->size1);
  gsl_matrix_set_zero(average);
  gsl_matrix_set_zero(sigma2_average);

  const double maxcomb=gsl_sf_choose((unsigned int)y->size2,split);
  if (maxcomb<1.*numb_estimates*COMB_THRESHOLD) n_est=(int)(maxcomb/(1.*COMB_THRESHOLD));
  else n_est=numb_estimates;  // avoid too high risk of repetition in random sampling of subsets
  //printf("NERCOME average: number of estimates %i for %g possible comb.\n",n_est,maxcomb);

  if (n_est<COMB_MIN) {  // just split original data matrix
    submatrix=gsl_matrix_submatrix(y,0,0,y->size1,split);
    gsl_matrix_memcpy(y1,&submatrix.matrix);  
    submatrix=gsl_matrix_submatrix(y,0,split,y->size1,y->size2-split);
    gsl_matrix_memcpy(y2,&submatrix.matrix);

    // get covariance estimate in cov; sample covariance for subset 2 in sigma2
    nercome_split_estimator(y1,y2,cov,sigma2); 
    return;
  }

  for (i=0;i<n_est;i++) {
    gsl_ran_choose(rnd,comb,split,sequence,(int)y->size2,sizeof(int));  // comb contains random ordered combination of split column indices -- note repetition of combinations is possible if maxcomb is small (sequence is global variable)

    index1=0;
    index2=0;
    for (j=0;j<y->size2;j++) {
      column=gsl_matrix_column(y,j);               // extracts column j of y
      if ((index1<split)&&(j==comb[index1])) {   // add column to y1
	gsl_matrix_set_col(y1,index1,&column.vector); // adds columns to y1
	index1++;
      }
      else {                  // add column to y2
	gsl_matrix_set_col(y2,index2,&column.vector); // adds columns to y2
	index2++;
      }
    }

    // get covariance estimate in cov; sample covariance for subset 2 in sigma2
    nercome_split_estimator(y1,y2,cov,sigma2);

    gsl_matrix_add(average,cov);
    gsl_matrix_add(sigma2_average,sigma2);
  }
  gsl_matrix_scale(average,1./(1.*n_est));   // now contains averaged covariance
  gsl_matrix_scale(sigma2_average,1./(1.*n_est));   // now contains averaged sigma2 covariance

  free(comb);
  gsl_matrix_free(y1);
  gsl_matrix_free(y2);
  gsl_matrix_free(cov);
  gsl_matrix_free(sigma2);
  return;
}


// nercome estimator following Lam (2016)
void cov_shrink_nercome(gsl_matrix *y,cov_est_pars *par,gsl_matrix *cov,int *minsplit)
{
  int i,split;
  double criterion,dsplit,mincrit=1.e10;
  gsl_matrix *res=gsl_matrix_alloc(y->size1,y->size1);
  gsl_matrix *sigma2=gsl_matrix_alloc(y->size1,y->size1);

  *minsplit=0;
  double split_min=y->size2*AVOID_FRAC;          // set range for possible splits
  double split_max=y->size2*(1.-AVOID_FRAC);
  if (par->split_steps>1) dsplit=(split_max-split_min)/(par->split_steps-1.);
  else {
    dsplit=0.0;
    //split_min=y->size2/2.;  // place split in centre
    split_min=y->size2*2./3.;  // place split at 2/3 of range
  }

  for (i=0;i<par->split_steps;i++) {
    split=(int)(split_min+i*dsplit);
    if ((split<=1)||(split>=y->size2-1)) {
      printf("Warning in routine 'cov_shrink_nercome': subsets after splitting too small %i <= 1 or %i >= %i; choose smaller number of splits for given number of realisations.\n",split,split,(int)y->size2-1);
      continue;
    }
    nercome_average(y,split,par->numb_estimates,res,sigma2);
    gsl_matrix_sub(sigma2,res);   // sigma2 now contains difference between res and sigma2
    criterion=frobenius_norm_squared(sigma2);  // criterion for optimal split

    //printf("NERCOME minimisation: %i  %g\n",split,criterion);
    if ((criterion<mincrit)&&(criterion!=0.0)) {   // criterion=0 if no splits possible
      mincrit=criterion;
      *minsplit=split;
      gsl_matrix_memcpy(cov,res);
    }
  }

  if (*minsplit==0) gsl_matrix_memcpy(cov,res);  // assign result of last run if no splits at all
  printf("NERCOME minimisation result: %g at split %i of %i\n",mincrit,*minsplit,(int)y->size2);  
  gsl_matrix_free(res);
  gsl_matrix_free(sigma2);
  return;
}




// opens file
FILE *bj_fileopen(char *rw,char *name)
{
  FILE *dat;
  if ((dat=fopen(name,rw))==NULL) {
    printf("Could not open file %s\n",name);
    exit(-1);
  }
  return dat;
}


// returns no. of columns in file
int bj_get_file_columns(char *name)
{
  FILE *dat;
  int ncol,read,cursor=0;
  char dummy[MAX_CHAR_FILE];
  char *line = calloc(MAX_CHAR_FILE,sizeof(char));

  dat=bj_fileopen("r",name);
  ncol=0;
  fgets(line,MAX_CHAR_FILE,dat);
  while (sscanf(line+cursor, "%s%n", dummy, &read)==1) {
    cursor+=read;
    ncol++;
  }
  free(line);
  fclose(dat);
  return(ncol);
}


// returns no. of rows in file
int bj_get_file_rows(char *name)
{
  FILE *dat;
  int nrow=0;
  dat=bj_fileopen("r",name);
  while ((fscanf(dat,"%*[^\n]"), fscanf(dat,"%*c"))!=EOF) nrow++; 
  fclose(dat);
  return(nrow);
}


// allocates 2D array
double **bj_alloc_2D(int dim1,int dim2)
{
  int i;
  double **arr=calloc(dim1,sizeof(double *));
  for (i=0;i<dim1;i++) {
    arr[i]=calloc(dim2,sizeof(double));
  }
  return arr;
}


// frees 2D array
void bj_free_2D(double **arr,int dim1,int dim2)
{
  int i;
  for (i=0;i<dim1;i++) {
    free(arr[i]);
  }
  free(arr);
  return;
}


// inverse of symmetric matrix via Cholesky decomposition
int inverse_symm(gsl_matrix *in,gsl_matrix *out)
{
  #include<gsl/gsl_errno.h>
  gsl_set_error_handler_off();

  int i,flag=0;
  const int dim=(int)in->size1;
  gsl_matrix *check=gsl_matrix_alloc(dim,dim);
  gsl_vector *b=gsl_vector_alloc(dim);
  gsl_matrix_memcpy(check,in);  // do not overwrite 'in'

  flag=gsl_linalg_cholesky_decomp(check);

  if (!flag) {
    for (i=0;i<dim;i++) {
      gsl_vector_set_basis(b,i);
      gsl_linalg_cholesky_svx(check,b);
      gsl_matrix_set_col(out,i,b);
    }
  }
  gsl_matrix_free(check);
  gsl_vector_free(b);
  return(flag);
}


// correlated Gaussian variables, with random sampler provided externally
void gaussian_correlated_withrnd(double **res,int dim,int nsample,double *mean,gsl_matrix *cov,gsl_rng *rnd)
{
  int i,j,k;
  double sum;
  gsl_vector *varu=gsl_vector_alloc(dim);
  gsl_matrix *corr=gsl_matrix_alloc(dim,dim);

  // compute correlation matrix and decompose
  for (i=0;i<dim;i++) {
    for (j=0;j<dim;j++) {
      gsl_matrix_set(corr,i,j,gsl_matrix_get(cov,i,j)/sqrt(gsl_matrix_get(cov,i,i)*gsl_matrix_get(cov,j,j)));
    }
  }
  gsl_linalg_cholesky_decomp(corr);

  // sample
  for (i=0;i<nsample;i++) {
    for (j=0;j<dim;j++) {
      gsl_vector_set(varu,j,gsl_ran_gaussian(rnd,sqrt(gsl_matrix_get(cov,j,j))));  // uncorrelated variates
    }

    for (j=0;j<dim;j++) {
      sum=0.0;
      for (k=0;k<=j;k++) {
	sum+=gsl_matrix_get(corr,j,k)*gsl_vector_get(varu,k);  // correlate
      }
      res[j][i]=mean[j]+sum;
    }
  }

  // clean up
  gsl_matrix_free(corr);
  gsl_vector_free(varu);
  return;
}
