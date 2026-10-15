#include "ReplayResourceAdmission.h"
#include <iostream>
#include <stdexcept>
int main() {
  if(replay::validateMemoryAndSwap("25769803776", "0")!=25769803776ULL)
    throw std::runtime_error("valid cap rejected");
  for(const auto& pair : {std::pair<const char*,const char*>{"max","0"},
      {"25769803777","0"},{"1024","max"},{"1024","1"},
      {"0","0"},{"-1","0"},{"1.0","0"},{"18446744073709551616","0"}}) {
    bool rejected=false;
    try { replay::validateMemoryAndSwap(pair.first,pair.second); }
    catch(const std::runtime_error&) { rejected=true; }
    if(!rejected)throw std::runtime_error("unsafe resource limit accepted");
  }
  std::cout<<"native resource admission cases pass\n";
}
