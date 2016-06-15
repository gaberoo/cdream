#include <cstdlib>
#include <iostream>
#include <fstream>
using namespace std;

#include "dream.h"
#include <rng/GSLStream.h>

#include <rapidjson/document.h>

typedef struct {
  double a[10];
} fpars;

double f(int chain, int gen, const double* state, const void* pars, bool recalc) {
  fpars& p = *(fpars*) pars;
  double ll = 0.0;
  for (int i = 0; i < 10; ++i) ll -= p.a[i]*state[i]*state[i];
  return ll;
}

int main() {
  ifstream in("tests/test.json");
  string json_input;
  in.seekg(0,ios::end);
  json_input.reserve(in.tellg());
  in.seekg(0,ios::beg);
  json_input.assign(istreambuf_iterator<char>(in), 
                    istreambuf_iterator<char>());

  rapidjson::Document jpars;
  dream_pars p;
  dream_pars_default(&p);

  try {
    jpars.Parse<0>(json_input.c_str());
    if (! jpars.IsObject()) throw 10;
    dream_pars_read_json(&p,jpars);
  } catch (int e) {
    cerr << "Coudln't create JSON document:" << endl;
    cerr << json_input << endl;
  } catch (const char* str) {
    cerr << "JSON exception: " << str << endl;
  }

  size_t n = p.nvar;

  fpars A;
  for (int i = 0; i < n; ++i) A.a[i] = 1.0;
  
  p.fun = &f;
  p.funPars = &A;

  rng::GSLStream rng;
  rng.alloc();

  dream(&p,&rng);

  dream_pars_free_vars(&p);

  return 0;
}

