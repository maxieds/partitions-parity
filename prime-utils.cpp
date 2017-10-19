/** 
 * prime-utils.cpp : Prime number utility workhorse for the program. 
 * Author: Maxie D. Schmidt (maxieds@gmail.com, mschmidt34@gatech.edu)
 * Created: 2017.10.18 
 **/ 

#include "utils.h"
#include "prime-utils.h"

double PrimesInInterval::ppow = 2.0;
double PrimesInInterval::epsilon = PrimesInInterval::p2epsilon(PrimesInInterval::ppow); 
PrimesInInterval PrimesInIntervalInst(PrimesInInterval::ppow); 

int PrimesInInterval::Indicator::I1(lint_t j, lint_t pgap) { 

     //return pgap % (72 * j) == 12 * j * (3 * j + 1) % (72 * j);
     return Utils::IsInteger(pgap - 12 * j * (3 * j + 1), 72 * j); 

} 

int PrimesInInterval::Indicator::I2(lint_t j, lint_t pgap) { 

     //return pgap % (36 * j + 12) == 0; 
     return Utils::IsInteger(pgap - 12 * j * (3 * j + 1), 24 * (3 * j + 1)); 

} 

int PrimesInInterval::Indicator::I3(lint_t j, lint_t pgap) { 

     //return pgap % (72 * j - 2) == 6 * j * (3 * j + 1) % (72 * j - 2); 
     return Utils::IsInteger(pgap - 12 * j * (3 * j + 1), 2 * (36 * j - 2)); 

} 

int PrimesInInterval::Indicator::I4(lint_t j, lint_t pgap) { 

     //return pgap % (864 * j + 312) == 6 * j * (3 * j + 1) % (864 * j + 312); 
     return Utils::IsInteger(pgap - 12 * j * (3 * j + 1), 24 * (36 * j + 13)); 

} 

int PrimesInInterval::Indicator::indicator_sum(lint_t j, lint_t pgap) { 

     return PrimesInInterval::Indicator::I1(j, pgap) + \
            PrimesInInterval::Indicator::I2(j, pgap) + \
            PrimesInInterval::Indicator::I3(j, pgap) + \
            PrimesInInterval::Indicator::I4(j, pgap); 
            
} 

const bool PrimesInInterval::SieveE::is_prime(int n, int *esieve_array, int minp) { 

     if(n < 1)
          return false; 
          
     for(int p = MAX(2, minp); p * p <= n; p++) { 
     
          if(esieve_array[p - 1]) { 
               for(int i = 2 * p; i <= n; i += p)
                   esieve_array[i - 1] = 0;
          }
          
     }
     return esieve_array[n - 1]; 

}

double PrimesInInterval::p2epsilon(double p) { 

     return 1.0 / p - 0.5; 
     
} 

double PrimesInInterval::epsilon2p(double eps) { 

     return 2.0 / (1.0 + 2.0 * eps); 
     
} 

double PrimesInInterval::t2x(int t) { 

     return CEILING(POW(t + 1, ppow)); 
     
} 

int PrimesInInterval::x2t(double x) { 

     return (int) FLOOR(POW(x, 0.5 + epsilon) - 1.0); 
     
} 

int PrimesInInterval::Qx_size_from_x(double x) { 

     return (int) FLOOR(POW(x, 0.5 + epsilon)); 
     
} 

int PrimesInInterval::Qx_size_from_t(int t) { 

     return PrimesInInterval::Qx_size_from_x(PrimesInInterval::t2x(t)); 
     
} 

bool PrimesInInterval::init_arrays(int asize, bool copy_elts) { 

     int *old_esieve_array = esieve_array; 
     esieve_array = (int *) malloc(asize * sizeof(int)); 
     memset(esieve_array, 1, asize * sizeof(int)); 
     if(copy_elts && old_esieve_array != NULL) { 
          for(int n = 0; n < esarr_size; n++) 
               esieve_array[n] = old_esieve_array[n];
     } 
     if(old_esieve_array != NULL)
          free(old_esieve_array); 
     esarr_size = asize; 
     
     lint_t *old_primes_array = primes_array; 
     primes_array = (lint_t *) malloc(asize * sizeof(lint_t)); 
     memset(primes_array, 0, asize * sizeof(lint_t)); 
     if(copy_elts && old_primes_array != NULL) { 
          for(int n = 0; n < parr_size; n++) 
               primes_array[n] = old_primes_array[n];
     } 
     if(old_primes_array != NULL) 
          free(old_primes_array); 
     parr_size = 0; 
     nsize = asize; 
     
     return true; 

}

bool PrimesInInterval::copy(const PrimesInInterval &src) { 

     if(nsize > 0) { 
          free(esieve_array); 
          free(primes_array); 
     }
     int tmax = src.nsize; 
     esieve_array = (int *) malloc(tmax * sizeof(int)); 
     memcpy(esieve_array, src.esieve_array, tmax * sizeof(int)); 
     esarr_size = src.esarr_size; 
     primes_array = (lint_t *) malloc(tmax * sizeof(lint_t)); 
     memcpy(primes_array, src.primes_array, tmax * sizeof(lint_t)); 
     parr_size = src.parr_size; 
     nsize = tmax; 

     return true; 
     
}

bool PrimesInInterval::clear() { 

     if(esieve_array == NULL || primes_array == NULL)
          return false; 
     
     free(esieve_array);
     esieve_array = NULL;  
     esarr_size = 0; 
     
     free(primes_array);
     primes_array = NULL;  
     parr_size = 0; 
     nsize = 0; 

     return true; 

}

lint_t * PrimesInInterval::compute_prime_interval(int t, int &isize) { 

     int pmin = (int) CEILING(POW(t, ppow)), pmax = (int) FLOOR(POW(t + 1, ppow)); 
     isize = pmax - pmin + 1; 
     if(nsize < isize + pmin) { // resize the storage arrays: 
          PrimesInInterval::init_arrays(2 * (isize + pmin), true); 
     } 
     int qidx = 0; 
     lint_t *primes = (lint_t *) malloc(isize * sizeof(lint_t)); 
     for(int q = 0; q < MIN(isize, nsize - pmin); q++) { 
          bool primeQ = PrimesInInterval::SieveE::is_prime(pmin + q, esieve_array); 
          if(primeQ)
               primes[qidx++] = pmin + q; 
     } 
     
     if(qidx == 0) { 
          isize = 0; 
          free(primes); 
          return NULL; 
     } 
     
     // resize the returned interval: 
     lint_t *rprimes = (lint_t *) malloc(qidx * sizeof(lint_t)); 
     memcpy(rprimes, primes, qidx * sizeof(lint_t));
     free(primes); 
     isize = qidx; 
     
     return rprimes; 

} 

lint_t * PrimesInInterval::compute_prime_gaps(int tmax, int &asize) { 

     if(tmax > parr_size) 
          fprintf(stderr, "BIG PROBLEM!\n"); 
     int pgidx = 0; 
     asize = tmax / 2.0 * (tmax - 1); 
     lint_t *prime_gaps = (lint_t *) malloc(asize * sizeof(lint_t)); 
     fprintf(stdout, "&& "); Utils::print_array(primes_array, parr_size); 
     fprintf(stdout, "\n"); 
     for(int i = 0; i < tmax; i++) { 
          for(int s = i + 1; s < tmax; s++) 
               prime_gaps[pgidx++] = primes_array[s] - primes_array[i]; 
     } 
     
     fprintf(stdout, "Prime Gaps(%d): ", tmax); Utils::print_array(prime_gaps, asize); 
     fprintf(stdout, "\n"); 
     
     return prime_gaps; 
     
} 

lint_t PrimesInInterval::compute_weight(int maxt, lint_t prime_prospect) { 

     primes_array[maxt] = prime_prospect; 
     int pgarr_size = 0; 
     lint_t *prime_gaps = PrimesInInterval::compute_prime_gaps(maxt + 1, pgarr_size); 
     int isum = 0; 
     //lint_t dmax = Utils::find_array_max(prime_gaps, pgarr_size); 
     for(int g = 0; g < pgarr_size; g++) { 
          for(int j = 1; j <= FLOOR(SQRT(prime_gaps[g])); j++) 
               isum += PrimesInInterval::Indicator::indicator_sum(j, prime_gaps[g]); 
     } 
     
     return isum; 
     
} 

lint_t PrimesInInterval::set_next_prime(int t) { 

     if(t < 1) 
          return 0; 
     else if(t == 1) { 
          primes_array[0] = 3; 
          ++parr_size; 
          return 3; 
     } 
     else if(t >= nsize) { 
          PrimesInInterval::init_arrays(2 * nsize, true); 
     } 
     
     int ptarr_size = 0; 
     lint_t *primest = PrimesInInterval::compute_prime_interval(t, ptarr_size); 
     lint_t *pweights = (lint_t *) malloc(ptarr_size * sizeof(lint_t)); 
     lint_t *pweights_temp = (lint_t *) malloc(ptarr_size * sizeof(lint_t)); // ERROR HERE!! 
     for(int p = 0; p < ptarr_size; p++) { 
          pweights[p] = PrimesInInterval::compute_weight(t - 1, primest[p]); 
     } 
     memcpy(pweights_temp, pweights, ptarr_size * sizeof(lint_t)); 
     lint_t min_weight = Utils::find_array_min(pweights, ptarr_size); 
     lint_t best_prime = 0; 
     for(int i = 0; i < ptarr_size; i++) { 
          if(pweights_temp[i] == min_weight) { 
               best_prime = primest[i]; 
               break; 
          } 
     } 
     free(primest); 
     free(pweights); 
     free(pweights_temp); 
     primes_array[parr_size++] = best_prime; 
     
     return best_prime; 
     
} 

PrimesInInterval::PrimesInInterval(double ppowp, int tmax) : 
     esieve_array(NULL), primes_array(NULL), esarr_size(0), parr_size(0), nsize(0) { 
     
     ppow = ppowp; 
     epsilon = PrimesInInterval::p2epsilon(ppow); 
     PrimesInInterval::init_arrays(tmax, false); 
     for(int t = 1; t <= tmax; t++) 
          PrimesInIntervalInst.set_next_prime(t); 
     
} 

PrimesInInterval::PrimesInInterval(const PrimesInInterval &rhs) { 

     PrimesInInterval::copy(rhs); 
     
} 

PrimesInInterval::~PrimesInInterval() { 

     PrimesInInterval::clear(); 
     
} 

PrimesInInterval & PrimesInInterval::operator=(const PrimesInInterval &rhs) { 

     if(this != &rhs) { 
          PrimesInInterval::clear(); 
          PrimesInInterval::copy(rhs); 
     } 
     return *this; 

} 

const lint_t PrimesInInterval::operator[](int t) { 

     if(t < 1) 
          return 0; 
     else if(t <= nsize) 
          return primes_array[t - 1]; 
     else { 
          PrimesInInterval::init_arrays(2 * nsize, true); 
          return primes_array[t - 1]; 
     } 

} 

int PrimesInInterval::size() const { 

     return parr_size; 
     
} 

lint_t * PrimesInInterval::get_primes(int tmax, int &asize) { 
     
     return PrimesInInterval::get_primes(1, tmax, asize); 

} 

lint_t * PrimesInInterval::get_primes(int tmin, int tmax, int &asize) { 
     
     if(tmin > tmax || tmin <= 0 || tmax <= 0)
          return NULL; 
     asize = tmax - tmin + 1; 
     lint_t *primes = (lint_t *) malloc(asize * sizeof(lint_t)); 
     memcpy(primes, &(primes_array[tmin - 1]), asize * sizeof(lint_t)); 
     return primes; 

} 

lint_t PrimesInInterval::compute_Iex(int t) { 

     int pgarr_size = 0; 
     lint_t *prime_gaps = PrimesInIntervalInst.get_primes(t, pgarr_size); 
     //lint_t dmax = Utils::find_array_max(prime_gaps, pgarr_size); 
     lint_t isum = 0; 
     for(int pidx = 0; pidx < pgarr_size; pidx++) { 
          for(int j = 1; j <= FLOOR(SQRT(prime_gaps[pidx])); j++) { 
               isum += PrimesInInterval::Indicator::indicator_sum(j, prime_gaps[pidx]); 
          } 
     } 
     return isum; 

} 

PIntRecord::PIntRecord() {} 

PIntRecord::PIntRecord(int tp, lint_t Iexp, lint_t Qx, lint_t *lastp) { 
     
     t = tp; 
     Iex = Iexp; 
     Qx_size = Qx; 
     ref_data[0] = ref_data[1] = NULL; 
     approx_x = PrimesInInterval::t2x(t); 
     logx = LOG10(approx_x); 
     PIntRecord::compute_ratio(); 
     if(lastp == NULL) 
          last_primes[0] = last_primes[1] = 0; 
     else 
          PIntRecord::set_last_primes(lastp[0], lastp[1]); 
     
} 

bool PIntRecord::set_last_primes(const lint_t &p1, const lint_t &p2) { 

     last_primes[0] = p1; 
     last_primes[1] = p2; 
     return true; 
     
} 

bool PIntRecord::set_last_primes(const lint_t * &primes, int asize) { 

     if(asize < 2) 
          return false; 
     PIntRecord::set_last_primes(primes[asize - 2], primes[asize - 1]); 
     return true; 
     
} 

const lint_t PIntRecord::get_last_prime(int pidx) const { 

     if(pidx != P1 && pidx != P2)
          return 0; 
     else 
          return last_primes[pidx]; 
          
}

const long double PIntRecord::compute_ratio() { 

     eratio = (long double) ((double) Qx_size) / (Iex == 0 ? 0.0 : Iex); 
     return eratio; 
     
} 

void PIntRecord::print() const  { 

     fprintf(stdout, "|| T: % 8d || X: % 12e || Qx: % 8d || Iex: % 10d || ERAT: % 10Lf ", 
                     t, logx, Qx_size, Iex, eratio); 
     fprintf(stdout, "|| ERATREC: % 10Lf || PIVAL: % 10f || P1: % 6Lf || P2: % 6Lf||\n", 
                     1.0 / eratio, approx_x / log(approx_x), 
                     (long double) Qx_size / POW(Iex, 0.5 + PrimesInInterval::epsilon), 
                     (long double) Qx_size / POW(Iex, PrimesInInterval::ppow)); 

} 

void PIntRecord::print_latex() const { 

     fprintf(stdout, "TODO: LaTeX printing not yet implemented.\n"); 
     
} 

PIntRecord PIntRecord::build_record(int t, lint_t *last_primes) { 
     
     PIntRecord prec(t, PrimesInIntervalInst.compute_Iex(t), 
                     PrimesInIntervalInst.Qx_size_from_t(t), last_primes); 
     return prec; 

} 

PIntRecord * PIntRecord::build_records(int tmin, int tmax, int &asize) { 

     if(tmin < 1 || tmax < 1 || tmin > tmax) 
          return NULL; 
     asize = tmax - tmin + 1; 
     PIntRecord *pints = (PIntRecord *) malloc(asize * sizeof(PIntRecord)); 
     
     for(int t = tmin; t <= tmax; t++) { 
          int lpasize = 0; 
          lint_t *last_primes = PrimesInIntervalInst.get_primes(MAX(2, t - 1), t, lpasize); 
          pints[t] = PIntRecord::build_record(t, last_primes); 
          free(last_primes); 
     } 
     return pints; 

} 



