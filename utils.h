/** 
 * utils.h : Basic non-specialty utilities for the program. 
 * Author: Maxie D. Schmidt (maxieds@gmail.com, mschmidt34@gatech.edu)
 * Created: 2017.10.18 
 **/ 

#ifndef __UTILS_H__
#define __UTILS_H__

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <unistd.h>
#include <limits.h>

#define MAX(x, y)          (x >= y ? x : y)
#define MIN(x, y)          (x <= y ? x : y)
#define FLOOR(x)           (floor(x))
#define CEILING(x)         (ceil(x))
#define SQRT(x)            (sqrt(x))
#define LOG(x)             (log(x))
#define LOG10(x)           (log10(x))
#define POW(x, ppow)       (pow(x, ppow))
#define DBL(j)             ((double) j)

typedef int lint_t; 

class Utils { 

     public: 
          static bool IsInteger(double x); 
          static bool IsInteger(int p, int q);
     
          static inline int compare(const void *e1, const void *e2, bool reverse = false); 
          static inline int compare_max(const void *e1, const void *e2); 
          static inline int compare_min(const void *e1, const void *e2); 

          static lint_t find_array_max(lint_t *arr, int asize); 
          static lint_t find_array_min(lint_t *arr, int asize); 
          
          static void print_array(lint_t *arr, int asize); 
          static inline void stamp(FILE *ofile = stdout); 

          static int array2mfile(const char *fpath, const lint_t * &arr, int asize); 
          static int array2mfile(const char *fpath, const double * &arr, int asize); 

}; 

#endif

