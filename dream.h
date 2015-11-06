#ifndef __DREAM_H__
#define __DREAM_H__

#include <cmath>
#include <cstdlib>
#include <string>
#include <map>
#include <vector>
#include <algorithm>
using namespace std;

#include <gsl/gsl_statistics.h>
#include <gsl/gsl_math.h>

#include <rapidjson/document.h>
#include <rng/RngStream.h>
#include "array.h"

typedef double (*LikFun)(const double* state, const void* pars);

typedef struct t_dream_pars {
  int vflag;                 /* vebose flag */
  int maxEvals;              /* max number of function evaluations */
  int optimAlg;              /* lock variables */
  int numChains;    
  string fn;                 /* input filename */
  string out_fn;             /* output filename */
  int appendFile;            /* continue from previous state */
  int report_interval;       /* report interval for state */
  int diagnostics;           /* report diagnostics at the end of the run */
  int burnIn;                /* number of steps for which to run an adaptive proposal size */
  int recalcLik;             /* recalculate likelihood of previously evaluated states */

  // DREAM variables
  int collapseOutliers;
  int gelmanEvals;
  int loopSteps;
  int deltaMax;
  int pCR_update;
  int nCR;

  double noise;
  double bstar_zero;
  double scaleReductionCrit;
  double reenterBurnin;

  int nfree;
  size_t nvar;
  double* varLo;
  double* varHi;
  double* varInit;
  int* varLock;
  string* varName;
  char* scale;

  LikFun fun;
  void* funPars;
} dream_pars;

void dream_pars_default(dream_pars* p);
void dream_pars_init_vars(dream_pars* p, size_t n);
void dream_pars_free_vars(dream_pars* p);

size_t dream_par_by_name(const dream_pars* p, string name);

// jpars.Parse<0>(json_input.c_str());
// assert(jpars.IsObject());
void dream_pars_read_json(dream_pars* p, rapidjson::Value& jpars);

int dream_restore_state(const dream_pars* p, Array3D<double>& state, Array2D<double>& lik, vector<double>& pCR, int& inBurnIn);
void dream_initialize(const dream_pars* p, rng::RngStream* rng, Array2DView<double>& state, ArrayView<double>& lik);
int dream(const dream_pars* p, rng::RngStream* rng);

void check_outliers(int t, Array2D<double>& lik, vector<double>& meanlik,
                    vector<bool>& outliers);
void gen_CR(rng::RngStream* rng, const vector<double>& pCR, 
            Array2D<int>& CRm, vector<unsigned>& L);
void gelman_rubin(Array3D<double>& state, vector<double>& scaleReduction, 
                  const int* lockVar, int first = 0, int numIter = -1, 
                  int adjustDF = 0);

#endif

