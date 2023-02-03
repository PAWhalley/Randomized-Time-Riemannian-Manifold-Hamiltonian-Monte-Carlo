//
//  normal.h
//  adaptive
//
//  Created by Ben on 27/10/2022.
//

#ifndef normal_h
#define normal_h

float complex c4_normal_01 (void );

double complex c8_normal_01 ( void);

int i4_normal_ab ( double a, double b );

float r4_normal_01 ( void);
float r4_normal_ab ( float a, float b );
float r4_uniform_01 (void );
void r4mat_print ( int m, int n, float a[], char *title );
void r4mat_print_some ( int m, int n, float a[], int ilo, int jlo, int ihi,
  int jhi, char *title );
void r4vec_print ( int n, float a[], char *title );
float *r4vec_uniform_01_new ( int n );

double r8_normal_01 (void );
double r8_normal_ab ( double a, double b );
double r8_uniform_01 (void );
void r8mat_normal_01 ( int m, int n, double x[] );
double *r8mat_normal_01_new ( int m, int n );
void r8mat_normal_ab ( int m, int n, double a, double b, double x[] );
double *r8mat_normal_ab_new ( int m, int n, double a, double b );
void r8mat_print ( int m, int n, double a[], char *title );
void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, char *title );
void r8vec_normal_01 ( int n, double x[] );
double *r8vec_normal_01_new ( int n );
void r8vec_normal_ab ( int n, double a, double b, double x[] );
double *r8vec_normal_ab_new ( int n, double a, double b );
void r8vec_print ( int n, double a[], char *title );
double *r8vec_uniform_01_new ( int n );

#endif /* normal_h */
