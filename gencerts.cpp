/** 
 * gencerts.cpp : Main program to generate the certificates for our parity 
 *                problem method of selecting |Q_x| primes in intervals. 
 * Author: Maxie D. Schmidt (maxieds@gmail.com, mschmidt34@gatech.edu)
 * Created: 2017.10.09 
 **/ 
 
#include <stdlib.h> 
#include <stdio.h> 
#include <unistd.h>
#include <getopt.h>
#include <string.h>
#include <time.h>
#include <vector>

#include "utils.h" 
#include "prime-utils.h"
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;
using std::vector; 

class runopts_t { 
     
     private: 
          void print_usage(char *pname) { 
               fprintf(stderr, "  Usage: %s [--tmin TMIN] [--tmax TMAX] [--fout FPATH]\n", pname); 
               memset(pname, 0x20, strlen(pname)); 
               fprintf(stderr, "         %s [--quiet] [--latex]\n", pname); 
          } 
          
          int parse_options(int argc, char **argv) { 
               int goc = -1; 
               while(true) { 
                    static struct option longopts[] = { 
                         {"tmin", required_argument, 0, '0'}, 
                         {"tmax", required_argument, 0, '1'}, 
                         {"fout", required_argument, 0, 'f'}, 
                         {"epsilon", required_argument, 0, 'e'}, 
                         {"quiet", no_argument, &quiet, 1}, 
                         {"latex", no_argument, &latex, 1}, 
                         {0, 0, 0, 0}
                    }; 
                    int option_index = 0; 
                    goc = getopt_long(argc, argv, "0:1:f:", longopts, &option_index); 
                    if(goc == -1)
                         break; 
                    bool print_usage_exit = false; 
                    switch(goc) { 
                    
                         case '0': 
                              tmin = atoi(optarg);
                              if(tmin < 1 || tmin > tmax) { 
                                   fprintf(stderr, "--tmin TMIN: argument must be >= 1, <= TMAX.\n"); 
                                   print_usage_exit = true; 
                              } 
                              break; 
                         case '1': 
                              tmax = atoi(optarg);
                              if(tmax < 1) { 
                                   fprintf(stderr, "--tmax TMAX: argument must be >= 1.\n"); 
                                   print_usage_exit = true; 
                              } 
                              break; 
                         case 'f': { 
                              char date_str[128], suffix_str[256]; 
                              const char *dash_str = (strlen(optarg) == 0 ? "" : "-"); 
                              time_t t = time(NULL); 
                              struct tm *tp = localtime(&t); 
                              strftime(date_str, sizeof(date_str), "", tp); 
                              snprintf(suffix_str, sizeof(suffix_str), 
                                       "%sgencert-%s.out", dash_str, date_str); 
                              int fpathlen = strlen(suffix_str) + strlen(optarg) + 1; 
                              fpath = (char *) malloc(fpathlen * sizeof(char)); 
                              //(*fpath)[0] = '\0'; 
                              strncpy(fpath, optarg, fpathlen); 
                              strncat(fpath, suffix_str, fpathlen); 
                              break; 
                         } 
                         case 'e': { 
                              double epsilon = atof(optarg); 
                              PrimesInInterval::epsilon = epsilon; 
                              PrimesInInterval::ppow = PrimesInInterval::epsilon2p(epsilon); 
                              break; 
                         } 
                         default: 
                              fprintf(stderr, "Unrecognized option (%s).\n", optarg); 
                              print_usage_exit = true; 
                         
                    } 
                    if(print_usage_exit) { 
                         print_usage(argv[0]); 
                         exit(-1); 
                    } 
               
               } 
               return 0; 
          } 
     
     public: 
          int tmin, tmax; 
          char *fpath; 
          int quiet, latex;
     
          runopts_t(int argc, char **argv) : 
               tmin(1), tmax(DEFAULT_TMAX), fpath(NULL), 
               quiet(false), latex(false) { 
               
               parse_options(argc, argv); 
          
          }
          
          ~runopts_t() { 
               if(fpath != NULL)
                    free(fpath); 
          } 
}; 

int main(int argc, char **argv) { 

     runopts_t runopts(argc, argv); 
     
     int tmin = runopts.tmin, tmax = runopts.tmax; 
     PrimesInIntervalInst.initialize(tmax); 
     int record_size = tmax - tmin + 1; 
     PIntRecord *precs = new PIntRecord[record_size]; 
     vector<double> tvals(record_size), error_ratios(record_size); 
     
     for(int t = tmin; t <= tmax; t++) { 
          precs[t - 1] = PIntRecord::build_record(t); 
          error_ratios.at(t - 1) = (double) precs[t - 1].eratio; 
          tvals.at(t - 1) = (double) t; 
          if(!runopts.quiet) { 
               precs[t - 1].print(); 
               //PrimesInIntervalInst.print_stats(); 
          } 
     } 
     
     // TODO: file path saving / output 
     
     if(!runopts.quiet) 
          fprintf(stdout, "\n"); 
     
     if(runopts.latex) { 
          fprintf(stdout, "LaTeX printing option is currently not implemented. (TODO)\n"); 
     } 
     
     plt::named_plot("|Q_x| / I_{e,x}", tvals, error_ratios, "g-");
     plt::title("Error Ratio Terms (Bdd. and Oscillating)");
     plt::xlim(350, tmax);
     plt::legend(); 
     plt::save("./error-ratios.png"); 
     
     return 0; 
     
} // main 
