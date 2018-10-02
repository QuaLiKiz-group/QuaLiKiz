/*
 * Copyright (c) 2005-2013 Steven G. Johnson
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

/* This stand-alone program, which should be compiled and linked against
   FFTW (www.fftw.org) version 3 or later, is used to generate the clencurt.h
   file for pcubature.c.  You only need to run it if you want to do
   p-adaptive cubature with more than 8193 points per dimension.  See
   the README file for more information. */


#include <stdlib.h>
#include <stdio.h>
#include <fftw3.h>

extern long double cosl(long double x);

/* Program to generate tables of precomputed points and weights
   for Clenshaw-Curtis quadrature on an interval [-1, 1] with

        3, 5, 9, ..., 2^(m+1)+1, ..., 2^(M+1)+1

   points up to some given M.  Because the quadrature rules are
   mirror-symmetric, we only need to store 2^m+1 weights for each rule.

   Furthermore, the rules are nested, so we only need to store the
   points for the M rule and the points for the other rules are a subset.
   We store the points and weights in a permuted order P corresponding to a
   usage where we first evaluate m=0, then m=1, etc. until it is converged.

   In particular, for the m rule (2^m+1 weights w[j], j=0,1,...,2^m),
   the corresponding points are

       x[j] = +/- cos(pi * j / 2^(m+1))

   (Note that for x[2^m] = 0; this point must be specially handled
    so that it is not counted twice.)

   So, we precompute an array clencurt_x of length 2^M storing

       clencurt_x[j] = cos(pi * P_M(j) / 2^(M+1))

   for j = 0,1,...,2^M-1.  Then, for a given rule m, we use

      x[P_m(j)] = clencurt_x[j]

   for j = 0,1,...,2^m-1 and x=0 for j = 2^m.  P_m is the permutation
 
      P_0(j) = j
      P_m(j) = 2 * P_{m-1}(j)          if j < 2^(m-1)
               2 * (j - 2^(m-1)) + 1   otherwise 

   The 2^m+1 weights w are stored for m=0,1,..,M in the same order in an array
   clencurt_w of length M+2^(M+1), in order of m.  So, the weights for
   a given m start at clencurt_w + (m + 2^m - 1), in the same order as
   clencurt_x except that it starts with the weight for x=0.
*/

static int P(int m, int j)
{
     if (m == 0) return j;
     else if (j < (1<<(m-1))) return 2 * P(m-1,j);
     else return 2 * (j - (1<<(m-1))) + 1;
}

/***************************************************************************/
/* The general principle is this: in Fejer and Clenshaw-Curtis quadrature,
   we take our function f(x) evaluated at f(cos(theta)), discretize
   it in theta to a vector f of points, compute the DCT, then multiply
   by some coefficients c (the integrals of cos(theta) * sin(theta)).
   In matrix form, given the DCT matrix D, this is:

             c' * D * f = (D' * c)' * f = w' * f

   where w is the vector of weights for each function value.  It
   is obviously much nicer to precompute w if we are going to be
   integrating many functions.   Since w = D' * c, and the transpose
   D' of a DCT is another DCT (e.g. the transpose of DCT-I is a DCT-I,
   modulo some rescaling in FFTW's normalization), to compute w is just
   another DCT.

   There is an additional wrinkle, because the odd-frequency DCT components
   of f integrate to zero, so every other entry in c is zero.  We can
   take advantage of this in computing w, because we can essentially do
   a radix-2 step in the DCT where one of the two subtransforms is zero.
   Therefore, for 2n+1 inputs, we only need to do a DCT of size n+1, and
   the weights w are a nice symmetric function.

   The weights are for integration of functions on (-1,1).
*/

void clencurt_weights(int n, long double *w)
{
     fftwl_plan p;
     int j;
     long double scale = 1.0 / n;
     
     p = fftwl_plan_r2r_1d(n+1, w, w, FFTW_REDFT00, FFTW_ESTIMATE);
     for (j = 0; j <= n; ++j) w[j] = scale / (1 - 4*j*j);
     fftwl_execute(p);
     w[0] *= 0.5;
     fftwl_destroy_plan(p);
}
/***************************************************************************/

#define KPI 3.1415926535897932384626433832795028841971L

int main(int argc, char **argv)
{
  
  
  /* 
  The program constructs a Fortran module file.
  Because of limits of how many continuation lines and lines per row are allowed
  in old standard Fortran, we split the array up into sections, and then construct the final array
  at the very end. If you wish to skip this, then simply change the if statements below
  wherever it says "We can fit the data comfortably into one array" such that it's guaranteed
  to evaluate as true. This statement appears twice in the program. 
  */

  FILE *fp;
  fp = freopen("clencurt.f90", "w", stdout);
  int M = argc > 1 ? atoi(argv[1]) : 11;
  long double *w;
  int i, j, m, n, n_count, d_count, m_count;
  long double k;
  i = 1;
  /* Note, max_row_length and max_cont lines are chosen carefuly to be powers of 2 */
  /* It is recommended that these are *not* changed */
  int max_row_b = 2;
  int max_cont_b = 5;
  int max_row_length = 1<<max_row_b; /* 4 */
  int max_cont_lines = 1<<max_cont_b; /* 32 */
  int max_array_size = max_row_length * max_cont_lines; /*128 */
  int max_array_size_b = max_row_b + max_cont_b; /* 7 */
  char str[1024]; 
  strcpy(str, "_");
  
  int max_w_m = 5; /* Magic number; M can go up to 5 before we started needing extra arrays for w */
  
  if (argc > 2 || M < 0) 
  {
	  fprintf(stderr, "clencurt_gen usage: clencurt_gen [M]\n");
	  return EXIT_FAILURE;
  }

  w = (long double *) fftwl_malloc(sizeof(long double) * ((1<<(M+2)) - 2));
  if (!w) 
    {
	  fprintf(stderr, "clencurt_gen: out of memory\n");
	  return EXIT_FAILURE;
    }	  

  printf("!AUTOMATICALLY GENERATED -- DO NOT EDIT \n\n");
  printf("MODULE CLENCURT \n");
  printf("USE KIND \n");
  printf("IMPLICIT NONE \n");
  printf("INTEGER, PARAMETER :: clencurt_M = %d;\n\n", M);
  
  /* We can fit the data comfortably into one array */
  if((1<<M) <= max_array_size)
  {
    
    printf("REAL(KIND=DBL), DIMENSION(%d), PARAMETER :: clencurt_x = (/ & ! length 2^M \n", 1<<M);
    k = KPI / ((long double) (1<<(M+1)));
    for (j = 0; j < (1<<M); ++j)
    {
      printf("%s%#0.18LG_DBL%s", i==1 ? "& " : "", cosl(k*P(M,j)), 
        j == (1<<M)-1 ? " &\n" : (i==max_row_length ? ", &\n" : ", "));
        
      if(i==max_row_length) i = 0;
      i++;
    }
    printf("& /)\n\n");
  }
  
  /* Cannot fit comfortably in one array */
  else
  {
    /* First, we construct the values necessary into small arrays of size max_array_size */
    k = KPI / ((long double) (1<<(M+1)));
    for(n = 0; n < (1<<M)/max_array_size; n++)
    {
      i = 1;
      printf("REAL(KIND=DBL), DIMENSION(%d), PARAMETER :: clencurt_x_%d = (/ & ! length %d \n", 
        max_array_size, n, max_array_size);
      
      for (j = 0; j < max_array_size; ++j)
      {
	      printf("%s%#0.18LG_DBL%s", i==1 ? "& " : "", cosl(k*P(M,j + n*max_array_size)), 
          j == max_array_size-1 ? " &\n" : (i==max_row_length ? ", &\n" : ", "));
          
        if(i==max_row_length) i = 0;
        i++;
        
      }
      printf("& /)\n\n");
    }
    
    /* Finished constructing all arrays, now need to contactenate them */
    n_count = (1<<M)/max_array_size;
    d_count = max_array_size * max_array_size;
    /* For high M, we need to contactenate and contactenate until we reach the desired length */
    /* Each layer of arrays will be of size max_array_size^a, where a is an integer */
    /* Because we've picked powers of two, this is easy */
    while(n_count > max_array_size)
    {
      for (n = 0; n < n_count/max_array_size; n++)
      {
        printf("REAL(KIND=DBL), DIMENSION(%d), PARAMETER :: clencurt_x_%s%d = (/ & ! length %d \n", 
          d_count, str, n, d_count);
        i = 1;
        for (j = 0; j < max_array_size; ++j)
        {
          printf("%sclencurt_x%s%d%s", i==1 ? "& " : "", str, j + n*max_array_size, 
            j == max_array_size-1 ? " &\n" : (i==max_row_length ? ", &\n" : ", " ));
            
          if(i==max_row_length) i = 0;
          i++;
        }
        printf("& /)\n\n");
      }
      strcat(str, "_");
      n_count = n_count / max_array_size;
      d_count = d_count * max_array_size;
    }
    /* Finally we can construct the final array */
    printf("REAL(KIND=DBL), DIMENSION(%d), PARAMETER :: clencurt_x = (/ & ! length 2^M \n", 1<<M);
    i = 1;
    for (n=0; n < n_count; n++)
    {
      printf("%sclencurt_x%s%d%s", i==1 ? "& " : "", str, 
        n, n== n_count - 1 ? " &\n" : (i==max_row_length ? ", &\n" : ", "));
        
      if(i==max_row_length) i = 0;
      i++;
    }
    printf("& /)\n\n");
    
  }
     
  /* We can fit the data comfortably into one array */
  if((M + (1<<(M+1)))<=max_array_size)
  {
    printf("REAL(KIND=DBL), DIMENSION(%d), PARAMETER :: clencurt_w = (/ & ! length M+2^(M+1) \n",
      M + (1<<(M+1)));
    
    for (m = 0; m <= M; ++m) 
    {
      clencurt_weights(1 << m, w);
      printf("& %#0.18LG_DBL, & ! m = %d\n", w[1 << m], m);
      i = 1;
      for (j = 0; j < (1 << m); ++j)
      {
        printf("%s%#0.18LG_DBL%s", i==1 ? "& ":"", w[P(m,j)], 
          j == (1<<m)-1 && m == M ? " &\n" : ((i==max_row_length || j ==((1<<m)-1))  ?  ", &\n" : ", "));
      
        if(i==max_row_length) i = 0;
        i++;
      }
     
    }
     printf("& /) \n\n");
     
  }
  /* Cannot fit comfortably into one array */
  else
  {
    /* Can easily fit up to m = max_w_m, so let's do that first */
    printf("REAL(KIND=DBL), DIMENSION(%d), PARAMETER :: clencurt_w_0 = (/ & ! length %d+2^(%d+1) \n",
      max_w_m + (1<<(max_w_m+1)), max_w_m, max_w_m);
    
    for (m = 0; m <= max_w_m; ++m) 
    {
      clencurt_weights(1 << m, w);
      printf("& %#0.18LG_DBL, & ! m = %d\n", w[1 << m], m);
      i = 1;
      for (j = 0; j < (1 << m); ++j)
      {
        printf("%s%#0.18LG_DBL%s", i==1 ? "& ":"", w[P(m,j)], 
          j == (1<<m)-1 && m == max_w_m ? " &\n" : ((i==max_row_length || j ==((1<<m)-1))  ?  ", &\n" : ", "));
      
        if(i==max_row_length) i = 0;
        i++;
      }
     
    }
    printf("& /) \n\n");
     
    /* Now, we note that for a given m, we need to construct a 2^m + 1 sized array */
    /* Up until 2^m = max_array_size, it's manageable for each m to be its own array */
    /* After that, we need to split them up */
    int m_up = max_array_size_b;
    if(M < m_up) m_up = M;
    for (m = max_w_m+1; m <= m_up; m++)
    {
      printf("REAL(KIND=DBL), DIMENSION(%d), PARAMETER :: clencurt_w_%d = (/ & ! length 1+2^(%d) \n",
       1 + (1<<(m)), m, m);
        
      clencurt_weights(1<<m, w);
      i = 1; 
      printf("& %#0.18LG_DBL, & ! m = %d\n", w[1 << m], m);
      for (j=0; j < (1<<m); j++)
      {
        printf("%s%#0.18LG_DBL%s", i==1 ? "& ":"", w[P(m,j)], 
          j == (1<<m)-1 ? " &\n" : ((i==max_row_length || j ==((1<<m)-1))  ?  ", &\n" : ", "));
      
        if(i==max_row_length) i = 0;
        i++;
      }
      printf("& /) \n\n");
      
      
    }
    m_count = m;
    /* Now m is too big to fit into one array */
    /* clencurt_w_m will be of size 1 + 2^(m) */
    /* therefore we split it just like we did before */
    for (m = m_count; m <= M; m++)
    {
      clencurt_weights(1<<m, w);
      /* Calculate how many times we split */
      /* Note by default max_array_size = 2^7 */
      int split = 1<<(m - max_array_size_b); /* Number of times we have to split */
      for (n = 0; n < split; n++)
      {
        if (n==0)
        {
          printf("REAL(KIND=DBL), DIMENSION(%d), PARAMETER :: clencurt_w_%d_%d = (/ & ! length 1+2^(%d) \n", 
            1 + max_array_size, m, n, max_array_size_b);
           printf("& %#0.18LG_DBL, & ! m = %d\n", w[1 << m], m);
        }
        else
        {
          printf("REAL(KIND=DBL), DIMENSION(%d), PARAMETER :: clencurt_w_%d_%d = (/ & ! length 2^(%d) \n", 
           max_array_size, m, n, max_array_size_b);
        }
        i = 1;
        for (j=0; j < max_array_size; j++)
        {
          printf("%s%#0.18LG_DBL%s", i==1 ? "& ":"", w[P(m,j+n*max_array_size)], 
            j == max_array_size-1 ? " &\n" : ((i==max_row_length || j ==((1<<m)-1))  ?  ", &\n" : ", "));
      
          if(i==max_row_length) i = 0;
          i++;
        }
        printf("& /) \n\n");
        
      }
      /* Splitting is done for a given m, now contactenate */
      /* For high m this needs to be contactenated in layers */
      d_count = max_array_size * max_array_size;
      int b_count = max_array_size_b *2;
      char str2[1024];
      strcpy(str2, "_");
      while(split > max_array_size)
      {
        
        for(n=0; n < split/max_array_size; n++)
        {
          if(n==0)
          {
            printf("REAL(KIND=DBL), DIMENSION(%d), PARAMETER :: clencurt_w_%d_%s%d = (/ & ! length 1+2^(%d) \n", 
              1 + d_count, m, str2, n, b_count);
          }
          else
          {
            printf("REAL(KIND=DBL), DIMENSION(%d), PARAMETER :: clencurt_w_%d_%s%d = (/ & ! length 2^(%d) \n", 
              d_count, m, str2, n, b_count);
          }
          i = 1;
          for (j=0; j < max_array_size; j++)
          {
            printf("%sclencurt_w_%d%s%d%s", i==1 ? "& " : "", m, str2, j+n*max_array_size,
              j==max_array_size-1 ? " &\n" : (i==max_row_length ? ", &\n" : ", "));
            if(i==max_row_length) i = 0;
            i++;
          }
          printf("& /) \n\n");
          
           
        }
        split = split/max_array_size;
        d_count = d_count * max_array_size;
        b_count = b_count * 2;
        strcat(str2, "_");
         
      }
       
      printf("REAL(KIND=DBL), DIMENSION(%d), PARAMETER :: clencurt_w_%d = (/ & ! length 1+2^(%d) \n", 
        1 + (1<<m), m, m);
      i = 1;
      for (n=0; n<split; n++)
      {
        printf("%sclencurt_w_%d%s%d%s", i==1 ? "& " : "", m, str2,
          n, n == split-1 ? " &\n" : (i==max_row_length ? ", &\n" : ", "));
        
        if(i==max_row_length) i = 0;
        i++;
      }
      printf("& /) \n\n");
    }
     
    /* Finally combined everything together */
    printf("REAL(KIND=DBL), DIMENSION(%d), PARAMETER :: clencurt_w = (/ & ! length M+2^(M+1) \n",
      M + (1<<(M+1)));
    printf("& clencurt_w_0, ");
    i = 2;
    for (j=max_w_m+1; j<=M; j++)
    {
      printf("%sclencurt_w_%d%s", i==1 ? "& " : "", j,
          j == M ? " &\n" : (i==max_row_length ? ", &\n" : ", "));
      if(i==max_row_length) i = 0;
      i++;
    }
    
    printf("& /) \n\n");
    
  }
  /* Finished */
  printf("\n\n END MODULE");

  fclose(fp);
  return EXIT_SUCCESS;
}
