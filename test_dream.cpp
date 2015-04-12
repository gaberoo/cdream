#include <cstdlib>
#include <iostream>
using namespace std;

#include "dream.h"
#include <rng/GSLStream.h>

typedef struct {
  double a[10];
} fpars;

double f(const double* state, const void* pars) {
  fpars& p = *(fpars*) pars;
  double ll = 0.0;
  for (int i = 0; i < 10; ++i) ll -= p.a[i]*state[i]*state[i];
  return ll;
}

int main() {
  size_t n = 10;

  dream_pars p;
  dream_pars_default(&p);
  dream_pars_init_vars(&p,n);

  fpars A;
  for (int i = 0; i < n; ++i) A.a[i] = 1.0;
  
  p.fun = &f;
  p.funPars = &A;

  p.vflag = 1;
  p.out_fn = "test";
  p.noise = 1.0;

  char str[100];
  for (int i = 0; i < n; ++i) {
    sprintf(str,"X%d",i);
    cout << str << endl;
    p.varName[i] = str;
    p.varLo[i] = -5.0;
    p.varHi[i] = 5.0;
  }

  p.gelmanEvals = 100;
  p.maxEvals = 30000;
  p.burnIn = 5000;

  rng::GSLStream rng;
  rng.alloc();

  dream(&p,&rng);

  dream_pars_free_vars(&p);
  return 0;
}

