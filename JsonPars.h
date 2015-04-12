#ifndef __JSONPARS_H__
#define __JSONPARS_H__

#include <rapidjson/document.h>
#include <string>

class JSONPars {
  JSONPars(std::string in) { parse(in); }
  virtual ~JSONPars() {}

  void parse(std::string in);

  protected:
    rapidjson::Document json;
};

inline void JSONPars::parse(std::string in) {
  json.Parse<0>(in.c_str());
  // assert(jpars.IsObject());
}

#endif

