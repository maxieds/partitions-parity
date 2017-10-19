/** 
 * prime-utils.h : Prime number utilities for the program. 
 * Author: Maxie D. Schmidt (maxieds@gmail.com, mschmidt34@gatech.edu)
 * Created: 2017.10.09 
 **/ 

#ifndef __PRIME_UTILS_H__
#define __PRIME_UTILS_H__

#include <math.h> 
#include <string.h>
#include <stdlib.h>
#include <time.h>

#define DEFAULTP          ((double) 1.96)
#define DEFAULT_TMAX      (12)
#define P1                (0)
#define P2                (1)
#define QX                (0)
#define IEX               (1)
#define ERATIO            (2)
#define PRIME             (2)

class PrimesInInterval { 

     public: 
          static double ppow, epsilon; 
                    
          class Indicator { 
          
               public: 
                    static inline int I1(lint_t j, lint_t pgap); 
                    static inline int I2(lint_t j, lint_t pgap); 
                    static inline int I3(lint_t j, lint_t pgap); 
                    static inline int I4(lint_t j, lint_t pgap); 
                    static inline int indicator_sum(lint_t j, lint_t pgap); 
               
          }; 
          
          const int is_prime(int n, int minp = 2); 

          static double p2epsilon(double p); 
          static double epsilon2p(double eps); 
          static double t2x(int t); 
          static int x2t(double x); 
          int Qx_size_from_x(double x); 
          int Qx_size_from_t(int t); 

     private: 
          bool initQ; 
          int *esieve_array; 
          lint_t *primes_array; 
          int esarr_size, parr_size, nsize; 
          
          bool init_arrays(int nsize, bool copy_elts = false); 
          bool copy(const PrimesInInterval &src); 
          bool clear(); 
          
          lint_t * compute_prime_interval(int t, int &isize); 
          lint_t * compute_prime_gaps(int maxt, int &asize, lint_t prime_prospect); 
          lint_t compute_weight(int maxt, lint_t prime_prospect); 
          lint_t set_next_prime(int t); 
          
     public: 
          PrimesInInterval(); 
          PrimesInInterval(const PrimesInInterval &rhs); 
          ~PrimesInInterval(); 
          
          bool initialize(int tmax = DEFAULT_TMAX); 
          PrimesInInterval & operator=(const PrimesInInterval &rhs); 
          const lint_t operator[](int t); 
          
          int size() const; 
          lint_t * get_primes(int tmax, int &asize); 
          lint_t * get_primes(int tmin, int tmax, int &asize); 
          
          lint_t compute_Iex(int t); 
          
          void print_stats() const; 

}; 

extern PrimesInInterval PrimesInIntervalInst; 

class PIntRecord { 

     public: 
          int t; 
          double approx_x, logx; 
          lint_t last_primes[2]; 
          lint_t Iex, Qx_size; 
          long double eratio; 
          
          static void *base_ref; 
          void *ref_data[2]; 
          
          PIntRecord(); 
          PIntRecord(int t, lint_t Iex, lint_t Qx, lint_t *lastp = NULL); 
          
          bool set_last_primes(const lint_t &p1, const lint_t &p2); 
          bool set_last_primes(const lint_t * &primes, int asize); 
          const lint_t get_last_prime(int pidx = P1) const; 
          const long double compute_ratio(); 
          
          void print() const; 
          void print_latex() const; 

          static PIntRecord build_record(int t, lint_t *last_primes = NULL); 
          static PIntRecord * build_records(int tmin, int tmax, int &asize); 
          static void * extract_records(const PIntRecord * &precs, int rid, int &asize); 

}; 

#endif

