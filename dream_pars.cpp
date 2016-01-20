#include "dream.h"

void dream_pars_default(dream_pars* p) {
  p->vflag = 0;
  p->maxEvals = 100000;
  p->optimAlg = 1;
  p->numChains = 5;
  p->fn = "";
  p->out_fn = "";
  p->appendFile = 0;
  p->report_interval = 1;
  p->diagnostics = 0;
  p->burnIn = 0;
  p->recalcLik = 0;
  p->noise = 0.05;
  p->bstar_zero = 1e-3;
  p->collapseOutliers = 1;
  p->gelmanEvals = 5000;
  p->loopSteps = 10;
  p->scaleReductionCrit = 1.01;
  p->deltaMax = 2;
  p->pCR_update = 1;
  p->nCR = 3;
  p->reenterBurnin = 0.2;
  p->fun = NULL;
  p->funPars = NULL;
}

// ---------------------------------------------------------------------------

void dream_pars_init_vars(dream_pars* p, size_t n) {
  p->nvar = n;
  p->nfree = n;
  p->varLo = (double*) calloc(n,sizeof(double));
  p->varHi = (double*) calloc(n,sizeof(double));
  p->varInit = (double*) calloc(n,sizeof(double));
  p->varLock = (int*) calloc(n,sizeof(int));
  p->varName = new string[n];
  p->scale = new char[n];
}

// ---------------------------------------------------------------------------

void dream_pars_free_vars(dream_pars* p) {
  free(p->varLo);
  free(p->varHi);
  free(p->varInit);
  free(p->varLock);
  delete[] p->varName;
  delete[] p->scale;
  p->nvar = 0;
  p->nfree = 0;
}

// ---------------------------------------------------------------------------

void dream_set_init(dream_pars* p, int n, 
                    const double* init, const string* name,
                    const int* lock, const double* lo, 
                    const double* hi, const char* scale) 
{
  p->nvar = n;
  p->nfree = n;
  memcpy(p->varInit, init,  n*sizeof(double));
  memcpy(p->varName, name,  n*sizeof(string));
  memcpy(p->varLock, lock,  n*sizeof(int));
  memcpy(p->varLo,   lo,    n*sizeof(double));
  memcpy(p->varHi,   hi,    n*sizeof(double));
  memcpy(p->scale,   scale, n*sizeof(char));
  for (int i = 0; i < n; ++i) if (lock[i]) --(p->nfree);
  p->deltaMax = (p->nfree-1)/2;
}

// ---------------------------------------------------------------------------

void dream_pars_from_json(dream_pars* p, rapidjson::Value& jpars) {
  // DEPRECATED !!!!!! 
  //  remove in next release...

  // reading priors from JSON
  rapidjson::Value::MemberIterator m1;
  rapidjson::SizeType _rj; 

  // find variable ranges
  rapidjson::Value::MemberIterator _d = jpars.FindMember("pars");
  if (_d == jpars.MemberEnd()) throw "pars";

  size_t nshifts = 0;
  rapidjson::Value::MemberIterator _s = jpars.FindMember("shifts");
  if (_s == jpars.MemberEnd()) throw "shifts";
  else nshifts = _s->value.Size();

  rapidjson::Value& d = _d->value;
  dream_pars_init_vars(p,(nshifts+1)*d.Size());

  // variables
  size_t j = 0;
  string name = "";
  double lo = 0.0;
  double hi = 0.0;
  char scale = 'n';

  for (rapidjson::SizeType i = 0; i < d.Size(); ++i) {
    m1 = d[i].FindMember("name"); 
    if (m1 != d[i].MemberEnd()) {
      if (! m1->value.IsString()) {
        name = "";
        throw "Bad varible name.";
      } else {
        name = m1->value.GetString();
      }
    }

    m1 = d[i].FindMember("limits"); 
    if (m1 != d[i].MemberEnd()) {
      if (! m1->value.IsArray()) throw "Bad variables limits.";
      if (m1->value.Size() != 2) throw "Bad variables limits.";
      _rj = 0; lo = m1->value[_rj].GetDouble();
      _rj = 1; hi = m1->value[_rj].GetDouble();
    }

    p->varName[j] = name;
    p->varLo[j] = lo;
    p->varHi[j] = hi;
    p->scale[i] = scale;

    m1 = d[i].FindMember("scale");
    if (m1 != d[i].MemberEnd()) {
      if (! m1->value.IsString()) throw "Bad scale.";
      scale = m1->value.GetString()[0];
    } else {
    }

    m1 = d[i].FindMember("lock");
    if (m1 != d[i].MemberEnd()) {
      if (m1->value.IsArray()) {
        if (m1->value.Size() != nshifts) throw "varSize";
        for (size_t k = 0; k < nshifts; ++k) {
          p->varName[j] = name;
          p->varLock[j] = 1;
          p->varHi[j] = m1->value[k].GetDouble();
          p->varLo[j] = p->varHi[i];
          p->varInit[j] = p->varHi[j];
          --(p->nfree);
          ++j;
        }
      } else if (m1->value.IsDouble()) {
        p->varName[j] = name;
        p->varLock[j] = 1;
        p->varHi[j] = m1->value.GetDouble();
        p->varLo[j] = p->varHi[j];
        p->varInit[j] = p->varHi[j];
        --(p->nfree);
        ++j;
      } else {
        throw "Locked variable isn't a double or an array.";
      }
    } else {
      m1 = d[i].FindMember("init");
      p->varInit[i] = m1->value.GetDouble();
      if (m1->value.IsArray()) {
        if (m1->value.Size() != nshifts) throw "varSize";
        for (size_t k = 0; k < nshifts; ++k) {
          p->varName[j] = name;
          p->varLock[j] = 0;
          p->varHi[j] = hi;
          p->varLo[j] = lo;
          p->varInit[j] = m1->value[k].GetDouble();
          ++j;
        }
      } else if (m1->value.IsDouble()) {
        p->varName[j] = name;
        p->varLock[j] = 0;
        p->varHi[j] = hi;
        p->varLo[j] = lo;
        p->varInit[j] = m1->value.GetDouble();
        ++j;
      } else {
        throw "Init variable isn't a double or an array.";
      }
    }
  }

  p->deltaMax = (p->nfree-1)/2;

}

// ---------------------------------------------------------------------------

void dream_pars_read_json(dream_pars* p, rapidjson::Value& jpars) {
  // reading priors from JSON
  rapidjson::Value::MemberIterator m1;
  rapidjson::SizeType _rj; 

  rapidjson::Value::MemberIterator dream = jpars.FindMember("dream");
  if (dream == jpars.MemberEnd()) throw "No dream section in JSON.";

  rapidjson::Value& dv = dream->value;

  m1 = dv.FindMember("prefix");
  if (m1 != dv.MemberEnd()) {
    if (! m1->value.IsString()) throw "prefix";
    p->out_fn = m1->value.GetString();
  }

  m1 = dv.FindMember("num_chains");
  if (m1 != dv.MemberEnd()) {
    if (! m1->value.IsInt()) throw "num_chains";
    p->numChains = m1->value.GetInt();
  }

  m1 = dv.FindMember("max_evals");
  if (m1 != dv.MemberEnd()) {
    if (! m1->value.IsInt()) throw "max_evals";
    p->maxEvals = m1->value.GetInt();
  }

  m1 = dv.FindMember("burn_in");
  if (m1 != dv.MemberEnd()) {
    if (! m1->value.IsInt()) throw "burn_in";
    p->burnIn = m1->value.GetInt();
  }

  m1 = dv.FindMember("recalc_lik");
  if (m1 != dv.MemberEnd()) {
    if (! m1->value.IsInt()) throw "recalc_lik";
    p->recalcLik = m1->value.GetInt();
  }

  m1 = dv.FindMember("gelman_evals");
  if (m1 != dv.MemberEnd()) {
    if (! m1->value.IsInt()) throw "gelman_evals";
    p->gelmanEvals = m1->value.GetInt();
  }

  m1 = dv.FindMember("vflag");
  if (m1 != dv.MemberEnd()) {
    if (! m1->value.IsInt()) throw "vflag";
    p->vflag = m1->value.GetInt();
  }

  m1 = dv.FindMember("noise");
  if (m1 != dv.MemberEnd()) {
    if (! m1->value.IsDouble()) throw "noise";
    p->noise = m1->value.GetDouble();
  }
}

// ---------------------------------------------------------------------------

size_t dream_par_by_name(const dream_pars* p, string name) {
  size_t i = 0;
  while (i < p->nvar) {
    if (p->varName[i] == name) break;
    ++i;
  }
  return i;
}

