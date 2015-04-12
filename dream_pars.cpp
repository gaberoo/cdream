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

void dream_pars_read_json(dream_pars* p, rapidjson::Value& jpars) {
  // reading priors from JSON
  rapidjson::Value::MemberIterator m1;
  rapidjson::SizeType _rj; 

  // find variable ranges
  rapidjson::Value::MemberIterator _d = jpars.FindMember("pars");
  if (_d == jpars.MemberEnd()) throw "pars";

  rapidjson::Value& d = _d->value;
  dream_pars_init_vars(p,d.Size());

  // variables
  for (rapidjson::SizeType i = 0; i < d.Size(); ++i) {
    m1 = d[i].FindMember("name"); 
    if (m1 != d[i].MemberEnd()) {
      if (! m1->value.IsString()) throw "Bad varible name.";
      else p->varName[i] = m1->value.GetString();
    }

    m1 = d[i].FindMember("limits"); 
    if (m1 != d[i].MemberEnd()) {
      if (! m1->value.IsArray()) throw "Bad variables limits.";
      if (! m1->value.Size() == 2) throw "Bad variables limits.";
      _rj = 0; p->varLo[i] = m1->value[_rj].GetDouble();
      _rj = 1; p->varHi[i] = m1->value[_rj].GetDouble();
    }

    m1 = d[i].FindMember("lock");
    if (m1 != d[i].MemberEnd()) {
      if (! m1->value.IsDouble()) {
        throw "Locked variable isn't a double.";
      } else {
        p->varLock[i] = 1;
        p->varLo[i] = p->varHi[i] = m1->value.GetDouble();
        --(p->nfree);
      }
    }

    m1 = d[i].FindMember("scale");
    if (m1 != d[i].MemberEnd()) {
      if (! m1->value.IsString()) throw "Bad scale.";
      p->scale[i] = m1->value.GetString()[0];
    } else {
      p->scale[i] = 'n';
    }
  }

  p->deltaMax = (p->nfree-1)/2;

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

size_t dream_par_by_name(const dream_pars* p, string name) {
  size_t i = 0;
  while (i < p->nvar) {
    if (p->varName[i] == name) break;
    ++i;
  }
  return i;
}

