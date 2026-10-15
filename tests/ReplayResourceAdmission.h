#pragma once
#include <cstdint>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>

namespace replay {
inline std::uint64_t validateMemoryAndSwap(const std::string& memory,
                                          const std::string& swap) {
  std::uint64_t value=0;
  if(memory.empty())throw std::runtime_error("missing memory limit");
  for(char c:memory) {
    if(c<'0'||c>'9'||value>(std::numeric_limits<std::uint64_t>::max()-(c-'0'))/10)
      throw std::runtime_error("finite memory limit required");
    value=value*10+(c-'0');
  }
  if(value==0 || value>24ULL*1024*1024*1024 || swap!="0")
    throw std::runtime_error("requires memory <= 24 GiB and swap = 0");
  return value;
}
inline void requireCurrentResourceAdmission() {
  std::ifstream membership("/proc/self/cgroup");
  std::string line,path;
  unsigned count=0;
  while(std::getline(membership,line))if(line.rfind("0::",0)==0) {
    path=line.substr(3);++count;
  }
  if(count!=1||path.empty()||path.front()!='/'||path.find("..")!=std::string::npos)
    throw std::runtime_error("unambiguous cgroup v2 path required");
  const auto read=[](const std::string& file) {
    std::ifstream stream(file);std::string value;
    if(!(stream>>value))throw std::runtime_error("cannot read cgroup resource limit");
    return value;
  };
  const std::string root="/sys/fs/cgroup"+path;
  validateMemoryAndSwap(read(root+"/memory.max"),read(root+"/memory.swap.max"));
}
}
