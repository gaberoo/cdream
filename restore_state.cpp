#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
using namespace std;

#include "dream.h"
#include "array.h"

int dream_restore_state(
    const dream_pars* p, 
    Array3D<double>& state, 
    Array2D<double>& lik,
    vector<double>& pCR,
    int& inBurnIn) 
{
  int prevLines = 0;
  int genNumber = 0;
  if (p->appendFile && (p->out_fn != "" || p->out_fn != "-")) {
    if (p->vflag) cerr << "Restoring previous state... ";
    prevLines = p->maxEvals;
    for (int i = 0; i < p->numChains; ++i) {
      if (p->vflag) cerr << i << " ";
      int line(0);
      ostringstream chain_fn("");
      chain_fn << p->out_fn << "." << i << ".txt";
      ifstream ifile(chain_fn.str().c_str());
      if (! ifile) {
        prevLines = 0;
        cerr << "files don't exists: " << chain_fn.str() << endl;
        break;
      } else {
        string input;
        while (! getline(ifile,input).eof()) {
          if (line >= p->maxEvals) break;
          if (input.length() == 0) continue;
          if (input[0] == '#') continue;
          istringstream istr(input);
          for (int j = 0; j < p->nvar; ++j) istr >> state(line,i,j); 
          istr >> lik(line,i) >> inBurnIn >> genNumber;
          for (int j(0); j < p->nCR; ++j) istr >> pCR[j];
          ++line;
        }
        ifile.close();
        if (prevLines > line) prevLines = line-1;
      }
      for (int j(0); j < p->nCR; ++j) cerr << pCR[j] << " ";
      if (p->vflag) cerr << endl;
    }
  }
  return prevLines;
}


