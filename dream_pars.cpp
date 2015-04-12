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
}

void dream_pars_free_vars(dream_pars* p) {
  free(p->varLo);
  free(p->varHi);
  free(p->varInit);
  free(p->varLock);
  delete[] p->varName;
  p->nvar = 0;
  p->nfree = 0;
}

// jpars.Parse<0>(json_input.c_str());
// assert(jpars.IsObject());
void dream_pars_read_json(dream_pars* p, rapidjson::Document& jpars) {
  /* Reading priors from JSON */
  rapidjson::Value::MemberIterator m1;
  rapidjson::SizeType _rj; 

  // find variable ranges
  rapidjson::Value::MemberIterator _d = jpars.FindMember("vars");
  assert(_d != jpars.MemberEnd());

  rapidjson::Value& d = _d->value;
  dream_pars_init_vars(p,d.Size());

  // get names
  map<string,int> nameMap;
  for (rapidjson::SizeType i = 0; i < d.Size(); ++i) {
    rapidjson::Value& a = d[i];
    m1 = a.FindMember("name"); 
    if (m1 != jpars.MemberEnd()) {
      assert(m1->value.IsString());
      p->varName[i] = m1->value.GetString();
      nameMap.insert(make_pair(p->varName[i],i));
    }
  }

  for (rapidjson::SizeType i = 0; i < d.Size(); ++i) {
    rapidjson::Value& a = d[i];

    m1 = a.FindMember("limits"); 
    if (m1 != jpars.MemberEnd()) {
      assert(m1->value.IsArray());
      assert(m1->value.Size() == 2);
      _rj = 0; p->varLo[i] = m1->value[_rj].GetDouble();
      _rj = 1; p->varHi[i] = m1->value[_rj].GetDouble();
    }

    m1 = a.FindMember("lock");
    if (m1 != jpars.MemberEnd()) {
      assert(m1->value.IsBool());
      p->varLock[i] = m1->value.GetBool();
      if (p->varLock[i]) --(p->nfree);
    }

    m1 = a.FindMember("fixed_to");
    if (m1 != jpars.MemberEnd()) {
      p->varLo[i] = m1->value.GetDouble();
      p->varHi[i] = m1->value.GetDouble();
    }
  }

  p->deltaMax = (p->nfree-1)/2;

  m1 = jpars.FindMember("prefix");
  if (m1 != jpars.MemberEnd()) {
    assert(m1->value.IsString());
    p->out_fn = m1->value.GetString();
  }

  m1 = jpars.FindMember("num_chains");
  if (m1 != jpars.MemberEnd()) {
    assert(m1->value.IsInt());
    p->numChains = m1->value.GetInt();
  }

  m1 = jpars.FindMember("max_evals");
  if (m1 != jpars.MemberEnd()) {
    assert(m1->value.IsInt());
    p->maxEvals = m1->value.GetInt();
  }

  m1 = jpars.FindMember("burn_in");
  if (m1 != jpars.MemberEnd()) {
    assert(m1->value.IsInt());
    p->burnIn = m1->value.GetInt();
  }

  m1 = jpars.FindMember("gelman_evals");
  if (m1 != jpars.MemberEnd()) {
    assert(m1->value.IsInt());
    p->gelmanEvals = m1->value.GetInt();
  }

}

