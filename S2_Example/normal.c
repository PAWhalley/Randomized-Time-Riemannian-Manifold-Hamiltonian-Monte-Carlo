//
//  normal.c
//  adaptive
//
//  Created by Ben on 27/10/2022.
//

# include <complex.h>
# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <time.h>

# include "normal.h"



float r4_normal_01 (void )

/******************************************************************************/
/*
  Purpose:

    r4_normal_01() returns a unit pseudonormal R4.

  Discussion:

    The standard normal probability distribution function (PDF) has
    mean 0 and standard deviation 1.

    The Box-Muller method is used, which is efficient, but
    generates two values at a time.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 September 2022

  Author:

    John Burkardt

  Output:

    float R4_NORMAL_01, a normally distributed random value.
*/
{
  float r1;
  float r2;
  const double r4_pi = 3.141592653589793;
  float x;

  r1 = ( float ) drand48 ( );
  r2 = ( float ) drand48 ( );
  x = sqrt ( - 2.0 * log ( r1 ) ) * cos ( 2.0 * r4_pi * r2 );

  return x;
}
/******************************************************************************/

float r4_normal_ab ( float a, float b )

/******************************************************************************/
/*
  Purpose:

    r4_normal_ab() returns a scaled pseudonormal R4.

  Discussion:

    The normal probability distribution function (PDF) is sampled,
    with mean A and standard deviation B.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 September 2022

  Author:

    John Burkardt

  Input:

    float A, the mean of the PDF.

    float B, the standard deviation of the PDF.

  Output:

    float R4_NORMAL_AB, a sample of the normal PDF.
*/
{
  float value;

  value = a + b * r4_normal_01 ( );

  return value;
}
/******************************************************************************/

void r4mat_print ( int m, int n, float a[], char *title )

/******************************************************************************/
/*
  Purpose:

    r4mat_print() prints an R4MAT.

  Discussion:

    An R4MAT is a doubly dimensioned array of R4 values, stored as a vector
    in column-major order.

    Entry A(I,J) is stored as A[I+J*M]

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 September 2022

  Author:

    John Burkardt

  Input:

    int M, the number of rows in A.

    int N, the number of columns in A.

    float A[M*N], the M by N matrix.

    char *TITLE, a title.
*/
{
  r4mat_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
/******************************************************************************/

void r4mat_print_some ( int m, int n, float a[], int ilo, int jlo, int ihi,
  int jhi, char *title )

/******************************************************************************/
/*
  Purpose:

    r4mat_print_some() prints some of an R4MAT.

  Discussion:

    An R4MAT is a doubly dimensioned array of R4 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 September 2022

  Author:

    John Burkardt

  Input:

    int M, the number of rows of the matrix.
    M must be positive.

    int N, the number of columns of the matrix.
    N must be positive.

    float A[M*N], the matrix.

    int ILO, JLO, IHI, JHI, designate the first row and
    column, and the last row and column to be printed.

    char *TITLE, a title.
*/
{
# define INCX 5

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "%s\n", title );

  if ( m <= 0 || n <= 0 )
  {
    fprintf ( stdout, "\n" );
    fprintf ( stdout, "  (None)\n" );
    return;
  }
/*
  Print the columns of the matrix, in strips of 5.
*/
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    if ( n < j2hi )
    {
      j2hi = n;
    }
    if ( jhi < j2hi )
    {
      j2hi = jhi;
    }

    fprintf ( stdout, "\n" );
/*
  For each column J in the current range...

  Write the header.
*/
    fprintf ( stdout, "  Col:  ");
    for ( j = j2lo; j <= j2hi; j++ )
    {
      fprintf ( stdout, "  %7d     ", j - 1 );
    }
    fprintf ( stdout, "\n" );
    fprintf ( stdout, "  Row\n" );
    fprintf ( stdout, "\n" );
/*
  Determine the range of the rows in this strip.
*/
    if ( 1 < ilo )
    {
      i2lo = ilo;
    }
    else
    {
      i2lo = 1;
    }
    if ( m < ihi )
    {
      i2hi = m;
    }
    else
    {
      i2hi = ihi;
    }

    for ( i = i2lo; i <= i2hi; i++ )
    {
/*
  Print out (up to) 5 entries in row I, that lie in the current strip.
*/
      fprintf ( stdout, "%5d:", i - 1 );
      for ( j = j2lo; j <= j2hi; j++ )
      {
        fprintf ( stdout, "  %14g", a[i-1+(j-1)*m] );
      }
      fprintf ( stdout, "\n" );
    }
  }

  return;
# undef INCX
}
/******************************************************************************/

void r4vec_print ( int n, float a[], char *title )

/******************************************************************************/
/*
  Purpose:

    r4vec_print() prints an R4VEC.

  Discussion:

    An R4VEC is a vector of R4's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 September 2022

  Author:

    John Burkardt

  Input:

    int N, the number of components of the vector.

    float A[N], the vector to be printed.

    char *TITLE, a title.
*/
{
  int i;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "%s\n", title );
  fprintf ( stdout, "\n" );
  for ( i = 0; i < n; i++ )
  {
    fprintf ( stdout, "  %8d: %14g\n", i, a[i] );
  }

  return;
}
/******************************************************************************/

float *r4vec_uniform_01_new ( int n )

/******************************************************************************/
/*
  Purpose:

    r4vec_uniform_01_new() returns a unit pseudorandom R4VEC.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 September 2022

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Springer Verlag, pages 201-202, 1983.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, pages 362-376, 1986.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, pages 136-143, 1969.

  Input:

    int N, the number of entries in the vector.

  Output:

    float R4VEC_UNIFORM_01_NEW[N], the vector of pseudorandom values.
*/
{
  int i;
  float *r;

  r = ( float * ) malloc ( n * sizeof ( float ) );

  for ( i = 0; i < n; i++ )
  {
    r[i] = drand48 ( );
  }

  return r;
}
/******************************************************************************/

double r8_normal_01 ( void)

/******************************************************************************/
/*
  Purpose:

    r8_normal_01() returns a unit pseudonormal R8.

  Discussion:

    The standard normal probability distribution function (PDF) has
    mean 0 and standard deviation 1.

    Because this routine uses the Box Muller method, it requires pairs
    of uniform random values to generate a pair of normal random values.
    This means that on every other call, the code can use the second
    value that it calculated.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 September 2022

  Author:

    John Burkardt

  Output:

    double R8_NORMAL_01, a normally distributed random value.
*/
{
  double r1;
  double r2;
  const double r8_pi = 3.141592653589793;
  double x;

  r1 = drand48 ( );
  r2 = drand48 ( );
  x = sqrt ( - 2.0 * log ( r1 ) ) * cos ( 2.0 * r8_pi * r2 );

  return x;
}
/******************************************************************************/

double r8_normal_ab ( double a, double b )

/******************************************************************************/
/*
  Purpose:

    r8_normal_ab() returns a scaled pseudonormal R8.

  Discussion:

    The normal probability distribution function (PDF) is sampled,
    with mean A and standard deviation B.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 September 2022

  Author:

    John Burkardt

  Input:

    double A, the mean of the PDF.

    double B, the standard deviation of the PDF.

  Output:

    double R8_NORMAL_AB, a sample of the normal PDF.
*/
{
  double value;

  value = a + b * r8_normal_01 ( );

  return value;
}
/******************************************************************************/

double r8_uniform_01 (void )

/******************************************************************************/
/*
  Purpose:

    r8_uniform_01() returns a unit pseudorandom R8.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 September 2022

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Springer Verlag, pages 201-202, 1983.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation
    edited by Jerry Banks,
    Wiley Interscience, page 95, 1998.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, pages 362-376, 1986.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, pages 136-143, 1969.

  Output:

    double R8_UNIFORM_01, a new pseudorandom variate, strictly between
    0 and 1.
*/
{
  double r;

  r = drand48 ( );

  return r;
}
/******************************************************************************/

void r8mat_normal_01 ( int m, int n, double r[] )

/******************************************************************************/
/*
  Purpose:

    r8mat_normal_01() returns a unit pseudonormal R8MAT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 September 2022

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Springer Verlag, pages 201-202, 1983.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, pages 362-376, 1986.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, pages 136-143, 1969.

  Input:

    int M, N, the number of rows and columns in the array.

  Output:

    double R[M*N], the array of pseudonormal values.
*/
{
  r8vec_normal_01 ( m * n, r );

  return;
}
/******************************************************************************/

double *r8mat_normal_01_new ( int m, int n )

/******************************************************************************/
/*
  Purpose:

    r8mat_normal_01_new() returns a unit pseudonormal R8MAT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 September 2022

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Springer Verlag, pages 201-202, 1983.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, pages 362-376, 1986.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, pages 136-143, 1969.

  Input:

    int M, N, the number of rows and columns in the array.

  Output:

    double R8MAT_NORMAL_01_NEW[M*N], the array of pseudonormal values.
*/
{
  double *r;

  r = r8vec_normal_01_new ( m * n );

  return r;
}
/******************************************************************************/

void r8mat_normal_ab ( int m, int n, double a, double b, double r[] )

/******************************************************************************/
/*
  Purpose:

    r8mat_normal_ab() returns a scaled pseudonormal R8MAT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 September 2022

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Springer Verlag, pages 201-202, 1983.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, pages 362-376, 1986.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, pages 136-143, 1969.

  Input:

    int M, N, the number of rows and columns in the array.

    double A, B, the mean and standard deviation.

  Output:

    double R[M*N], the array of pseudonormal values.
*/
{
  r8vec_normal_ab ( m * n, a, b, r );

  return;
}
/******************************************************************************/

double *r8mat_normal_ab_new ( int m, int n, double a, double b )

/******************************************************************************/
/*
  Purpose:

    r8mat_normal_ab_new() returns a scaled pseudonormal R8MAT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 September 2022

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Springer Verlag, pages 201-202, 1983.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, pages 362-376, 1986.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, pages 136-143, 1969.

  Input:

    int M, N, the number of rows and columns in the array.

    double A, B, the mean and standard deviation.

  Output:

    double R8MAT_NORMAL_AB_NEW[M*N], the array of pseudonormal values.
*/
{
  double *r;

  r = r8vec_normal_ab_new ( m * n, a, b );

  return r;
}
/******************************************************************************/

void r8mat_print ( int m, int n, double a[], char *title )

/******************************************************************************/
/*
  Purpose:

    r8mat_print() prints an R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    Entry A(I,J) is stored as A[I+J*M]

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 September 2022

  Author:

    John Burkardt

  Input:

    int M, the number of rows in A.

    int N, the number of columns in A.

    double A[M*N], the M by N matrix.

    char *TITLE, a title.
*/
{
  r8mat_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
/******************************************************************************/

void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, char *title )

/******************************************************************************/
/*
  Purpose:

    r8mat_print_some() prints some of an R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 September 2022

  Author:

    John Burkardt

  Input:

    int M, the number of rows of the matrix.
    M must be positive.

    int N, the number of columns of the matrix.
    N must be positive.

    double A[M*N], the matrix.

    int ILO, JLO, IHI, JHI, designate the first row and
    column, and the last row and column to be printed.

    char *TITLE, a title.
*/
{
# define INCX 5

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "%s\n", title );

  if ( m <= 0 || n <= 0 )
  {
    fprintf ( stdout, "\n" );
    fprintf ( stdout, "  (None)\n" );
    return;
  }
/*
  Print the columns of the matrix, in strips of 5.
*/
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    if ( n < j2hi )
    {
      j2hi = n;
    }
    if ( jhi < j2hi )
    {
      j2hi = jhi;
    }

    fprintf ( stdout, "\n" );
/*
  For each column J in the current range...

  Write the header.
*/
    fprintf ( stdout, "  Col:  ");
    for ( j = j2lo; j <= j2hi; j++ )
    {
      fprintf ( stdout, "  %7d     ", j - 1 );
    }
    fprintf ( stdout, "\n" );
    fprintf ( stdout, "  Row\n" );
    fprintf ( stdout, "\n" );
/*
  Determine the range of the rows in this strip.
*/
    if ( 1 < ilo )
    {
      i2lo = ilo;
    }
    else
    {
      i2lo = 1;
    }
    if ( m < ihi )
    {
      i2hi = m;
    }
    else
    {
      i2hi = ihi;
    }

    for ( i = i2lo; i <= i2hi; i++ )
    {
/*
  Print out (up to) 5 entries in row I, that lie in the current strip.
*/
      fprintf ( stdout, "%5d:", i - 1 );
      for ( j = j2lo; j <= j2hi; j++ )
      {
        fprintf ( stdout, "  %14g", a[i-1+(j-1)*m] );
      }
      fprintf ( stdout, "\n" );
    }
  }

  return;
# undef INCX
}
/******************************************************************************/

void r8vec_normal_01 ( int n, double x[] )

/******************************************************************************/
/*
  Purpose:

    r8vec_normal_01() returns a unit pseudonormal R8VEC.

  Discussion:

    The standard normal probability distribution function (PDF) has
    mean 0 and standard deviation 1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 September 2022

  Author:

    John Burkardt

  Input:

    int N, the number of values desired.

  Output:

    double X[N], a sample of the standard normal PDF.

  Local:

    double R[N+1], is used to store some uniform random values.
    Its dimension is N+1, but really it is only needed to be the
    smallest even number greater than or equal to N.

    int X_LO, X_HI, records the range of entries of
    X that we need to compute.
*/
{
  int i;
  int m;
  double *r;
  const double r8_pi = 3.141592653589793;
  int x_hi;
  int x_lo;
/*
  Record the range of X we need to fill in.
*/
  x_lo = 1;
  x_hi = n;
/*
  If we need just one new value, do that here to avoid null arrays.
*/
  if ( x_hi - x_lo + 1 == 1 )
  {
    r = r8vec_uniform_01_new ( 2 );

    x[x_hi-1] = sqrt ( - 2.0 * log ( r[0] ) ) * cos ( 2.0 * r8_pi * r[1] );

    free ( r );
  }
/*
  If we require an even number of values, that's easy.
*/
  else if ( ( x_hi - x_lo + 1 ) % 2 == 0 )
  {
    m = ( x_hi - x_lo + 1 ) / 2;

    r = r8vec_uniform_01_new ( 2*m );

    for ( i = 0; i <= 2*m-2; i = i + 2 )
    {
      x[x_lo+i-1] = sqrt ( - 2.0 * log ( r[i] ) ) * cos ( 2.0 * r8_pi * r[i+1] );
      x[x_lo+i  ] = sqrt ( - 2.0 * log ( r[i] ) ) * sin ( 2.0 * r8_pi * r[i+1] );
    }

    free ( r );
  }
/*
  If we require an odd number of values, we generate an even number,
  and handle the last pair specially, storing one in X(N), and
  saving the other for later.
*/
  else
  {
    x_hi = x_hi - 1;

    m = ( x_hi - x_lo + 1 ) / 2 + 1;

    r = r8vec_uniform_01_new ( 2 * m );

    for ( i = 0; i <= 2 * m - 4; i = i + 2 )
    {
      x[x_lo+i-1] = sqrt ( - 2.0 * log ( r[i] ) ) * cos ( 2.0 * r8_pi * r[i+1] );
      x[x_lo+i  ] = sqrt ( - 2.0 * log ( r[i] ) ) * sin ( 2.0 * r8_pi * r[i+1] );
    }

    i = 2*m - 2;

    x[x_lo+i-1] = sqrt ( - 2.0 * log ( r[i] ) ) * cos ( 2.0 * r8_pi * r[i+1] );

    free ( r );
  }

  return;
}
/******************************************************************************/

double *r8vec_normal_01_new ( int n )

/******************************************************************************/
/*
  Purpose:

    r8vec_normal_01_new() returns a unit pseudonormal R8VEC.

  Discussion:

    The standard normal probability distribution function (PDF) has
    mean 0 and standard deviation 1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 September 2022

  Author:

    John Burkardt

  Input:

    int N, the number of values desired.

  Output:

    double R8VEC_NORMAL_01_NEW[N], a sample of the standard normal PDF.

  Local:

    double R[N+1], is used to store some uniform random values.
    Its dimension is N+1, but really it is only needed to be the
    smallest even number greater than or equal to N.

    int X_LO, X_HI, records the range of entries of
    X that we need to compute.
*/
{
  int i;
  int m;
  double *r;
  const double r8_pi = 3.141592653589793;
  double *x;
  int x_hi;
  int x_lo;

  x = ( double * ) malloc ( n * sizeof ( double ) );
/*
  Record the range of X we need to fill in.
*/
  x_lo = 1;
  x_hi = n;
/*
  If we need just one new value, do that here to avoid null arrays.
*/
  if ( x_hi - x_lo + 1 == 1 )
  {
    r = r8vec_uniform_01_new ( 2 );

    x[x_hi-1] = sqrt ( - 2.0 * log ( r[0] ) ) * cos ( 2.0 * r8_pi * r[1] );

    free ( r );
  }
/*
  If we require an even number of values, that's easy.
*/
  else if ( ( x_hi - x_lo + 1 ) % 2 == 0 )
  {
    m = ( x_hi - x_lo + 1 ) / 2;

    r = r8vec_uniform_01_new ( 2*m );

    for ( i = 0; i <= 2*m-2; i = i + 2 )
    {
      x[x_lo+i-1] = sqrt ( - 2.0 * log ( r[i] ) ) * cos ( 2.0 * r8_pi * r[i+1] );
      x[x_lo+i  ] = sqrt ( - 2.0 * log ( r[i] ) ) * sin ( 2.0 * r8_pi * r[i+1] );
    }
    free ( r );
  }
/*
  If we require an odd number of values, we generate an even number,
  and handle the last pair specially, storing one in X(N), and
  saving the other for later.
*/
  else
  {
    x_hi = x_hi - 1;

    m = ( x_hi - x_lo + 1 ) / 2 + 1;

    r = r8vec_uniform_01_new ( 2 * m );

    for ( i = 0; i <= 2*m-4; i = i + 2 )
    {
      x[x_lo+i-1] = sqrt ( - 2.0 * log ( r[i] ) ) * cos ( 2.0 * r8_pi * r[i+1] );
      x[x_lo+i  ] = sqrt ( - 2.0 * log ( r[i] ) ) * sin ( 2.0 * r8_pi * r[i+1] );
    }

    i = 2 * m - 2;

    x[x_lo+i-1] = sqrt ( - 2.0 * log ( r[i] ) ) * cos ( 2.0 * r8_pi * r[i+1] );

    free ( r );
  }

  return x;
}
/******************************************************************************/

void r8vec_normal_ab ( int n, double b, double c, double x[] )

/******************************************************************************/
/*
  Purpose:

    r8vec_normal_ab() returns a scaled pseudonormal R8VEC.

  Discussion:

    The scaled normal probability distribution function (PDF) has
    mean A and standard deviation B.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 September 2022

  Author:

    John Burkardt

  Input:

    int N, the number of values desired.

    double B, C, the mean and standard deviation.

  Output:

    double X[N], a sample of the standard normal PDF.

  Local:

    double R[N+1], is used to store some uniform random values.
    Its dimension is N+1, but really it is only needed to be the
    smallest even number greater than or equal to N.

    int X_LO, X_HI, records the range of entries of
    X that we need to compute.
*/
{
  int i;
  int m;
  double *r;
  const double r8_pi = 3.141592653589793;
  int x_hi;
  int x_lo;
/*
  Record the range of X we need to fill in.
*/
  x_lo = 1;
  x_hi = n;
/*
  If we need just one new value, do that here to avoid null arrays.
*/
  if ( x_hi - x_lo + 1 == 1 )
  {
    r = r8vec_uniform_01_new ( 2 );

    x[x_hi-1] = sqrt ( - 2.0 * log ( r[0] ) ) * cos ( 2.0 * r8_pi * r[1] );

    free ( r );
  }
/*
  If we require an even number of values, that's easy.
*/
  else if ( ( x_hi - x_lo + 1 ) % 2 == 0 )
  {
    m = ( x_hi - x_lo + 1 ) / 2;

    r = r8vec_uniform_01_new ( 2 * m );

    for ( i = 0; i <= 2*m-2; i = i + 2 )
    {
      x[x_lo+i-1] = sqrt ( - 2.0 * log ( r[i] ) ) * cos ( 2.0 * r8_pi * r[i+1] );
      x[x_lo+i  ] = sqrt ( - 2.0 * log ( r[i] ) ) * sin ( 2.0 * r8_pi * r[i+1] );
    }

    free ( r );
  }
/*
  If we require an odd number of values, we generate an even number,
  and handle the last pair specially, storing one in X(N), and
  saving the other for later.
*/
  else
  {
    x_hi = x_hi - 1;

    m = ( x_hi - x_lo + 1 ) / 2 + 1;

    r = r8vec_uniform_01_new ( 2 * m );

    for ( i = 0; i <= 2*m-4; i = i + 2 )
    {
      x[x_lo+i-1] = sqrt ( - 2.0 * log ( r[i] ) ) * cos ( 2.0 * r8_pi * r[i+1] );
      x[x_lo+i  ] = sqrt ( - 2.0 * log ( r[i] ) ) * sin ( 2.0 * r8_pi * r[i+1] );
    }

    i = 2 * m - 2;

    x[x_lo+i-1] = sqrt ( - 2.0 * log ( r[i] ) ) * cos ( 2.0 * r8_pi * r[i+1] );

    free ( r );
  }

  for ( i = 0; i < n; i++ )
  {
    x[i] = b + c * x[i];
  }

  return;
}
/******************************************************************************/

double *r8vec_normal_ab_new ( int n, double b, double c )

/******************************************************************************/
/*
  Purpose:

    r8vec_normal_ab_new() returns a scaled pseudonormal R8VEC.

  Discussion:

    The scaled normal probability distribution function (PDF) has
    mean A and standard deviation B.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 September 2022

  Author:

    John Burkardt

  Input:

    int N, the number of values desired.

    double B, C, the mean and standard deviation.

  Output:

    double R8VEC_NORMAL_AB_NEW[N], a sample of the standard normal PDF.

  Local:

    double R[N+1], is used to store some uniform random values.
    Its dimension is N+1, but really it is only needed to be the
    smallest even number greater than or equal to N.

    int X_LO, X_HI, records the range of entries of
    X that we need to compute.
*/
{
  int i;
  int m;
  double *r;
  const double r8_pi = 3.141592653589793;
  double *x;
  int x_hi;
  int x_lo;

  x = ( double * ) malloc ( n * sizeof ( double ) );
/*
  Record the range of X we need to fill in.
*/
  x_lo = 1;
  x_hi = n;
/*
  If we need just one new value, do that here to avoid null arrays.
*/
  if ( x_hi - x_lo + 1 == 1 )
  {
    r = r8vec_uniform_01_new ( 2 );

    x[x_hi-1] = sqrt ( - 2.0 * log ( r[0] ) ) * cos ( 2.0 * r8_pi * r[1] );

    free ( r );
  }
/*
  If we require an even number of values, that's easy.
*/
  else if ( ( x_hi - x_lo + 1 ) % 2 == 0 )
  {
    m = ( x_hi - x_lo + 1 ) / 2;

    r = r8vec_uniform_01_new ( 2 * m );

    for ( i = 0; i <= 2*m-2; i = i + 2 )
    {
      x[x_lo+i-1] = sqrt ( - 2.0 * log ( r[i] ) ) * cos ( 2.0 * r8_pi * r[i+1] );
      x[x_lo+i  ] = sqrt ( - 2.0 * log ( r[i] ) ) * sin ( 2.0 * r8_pi * r[i+1] );
    }

    free ( r );
  }
/*
  If we require an odd number of values, we generate an even number,
  and handle the last pair specially, storing one in X(N).
*/
  else
  {
    x_hi = x_hi - 1;

    m = ( x_hi - x_lo + 1 ) / 2 + 1;

    r = r8vec_uniform_01_new ( 2 * m );

    for ( i = 0; i <= 2*m-4; i = i + 2 )
    {
      x[x_lo+i-1] = sqrt ( - 2.0 * log ( r[i] ) ) * cos ( 2.0 * r8_pi * r[i+1] );
      x[x_lo+i  ] = sqrt ( - 2.0 * log ( r[i] ) ) * sin ( 2.0 * r8_pi * r[i+1] );
    }

    i = 2*m - 2;

    x[x_lo+i-1] = sqrt ( - 2.0 * log ( r[i] ) ) * cos ( 2.0 * r8_pi * r[i+1] );

    free ( r );
  }

  for ( i = 0; i < n; i++ )
  {
    x[i] = b + c * x[i];
  }

  return x;
}
/******************************************************************************/

void r8vec_print ( int n, double a[], char *title )

/******************************************************************************/
/*
  Purpose:

    r8vec_print() prints an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 September 2022

  Author:

    John Burkardt

  Input:

    int N, the number of components of the vector.

    double A[N], the vector to be printed.

    char *TITLE, a title.
*/
{
  int i;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "%s\n", title );
  fprintf ( stdout, "\n" );
  for ( i = 0; i < n; i++ )
  {
    fprintf ( stdout, "  %8d: %14g\n", i, a[i] );
  }

  return;
}
/******************************************************************************/

double *r8vec_uniform_01_new ( int n )

/******************************************************************************/
/*
  Purpose:

    r8vec_uniform_01_new() returns a unit pseudorandom R8VEC.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 September 2022

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Springer Verlag, pages 201-202, 1983.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, pages 362-376, 1986.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, pages 136-143, 1969.

  Input:

    int N, the number of entries in the vector.

  Output:

    double R8VEC_UNIFORM_01_NEW[N], the vector of pseudorandom values.
*/
{
  int i;
  double *r;

  r = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    r[i] = drand48 ( );
  }

  return r;
}
