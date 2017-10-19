/** 
 * utils.cpp : Basic non-specialty utilities for the program. 
 * Author: Maxie D. Schmidt (maxieds@gmail.com, mschmidt34@gatech.edu)
 * Created: 2017.10.18 
 **/ 

#include "utils.h" 

bool Utils::IsInteger(double x) { 
     
     double ipart; 
     return modf(x, &ipart) == 0.0; 
     
} 

bool Utils::IsInteger(int p, int q) { 

     return Utils::IsInteger(DBL(p) / DBL(q)); 

}

int Utils::compare(const void *ep1, const void *ep2, bool reverse) { 
     
     lint_t e1 = *((lint_t *) ep1); 
     lint_t e2 = *((lint_t *) ep2); 
     int cmp_value = 0;
     if(e1 < e2)
          cmp_value = -2; 
     else if(e1 == e2) 
          cmp_value = 0; 
     else
          cmp_value = 1; 
          
     if(!reverse || cmp_value == 0) 
          return cmp_value; 
     else
          return ~cmp_value; 
          
} 

int Utils::compare_max(const void *ep1, const void *ep2) { 
     
     return Utils::compare(ep1, ep2, true); 
     
} 

int Utils::compare_min(const void *ep1, const void *ep2) { 
     
     return Utils::compare(ep1, ep2, false); 
     
} 

lint_t Utils::find_array_max(lint_t *arr, int asize) { 

     if(asize == 0) 
          return INT_MAX; 
     
     qsort(arr, asize, sizeof(lint_t), Utils::compare_max); 
     return arr[0]; 
     
}

lint_t Utils::find_array_min(lint_t *arr, int asize) { 

     if(asize == 0) 
          return INT_MIN; 
     
     qsort(arr, asize, sizeof(lint_t), Utils::compare_min); 
     return arr[0]; 
     
}

void Utils::print_array(lint_t *arr, int asize) { 

     fprintf(stdout, "{"); 
     for(int i = 0; i < asize - 1; i++) { 
          fprintf(stdout, "%d, ", arr[i]); 
     } 
     if(asize - 1 >= 0) 
          fprintf(stdout, "%d", arr[asize - 1]); 
     fprintf(stdout, "}"); 

} 

void Utils::stamp(FILE *ofile) { 

     fprintf(ofile, "Stamp: \"%s:%d\" (in func %s)\n", __FILE__, __LINE__, __FUNCTION__);
     fflush(ofile); 

} 


