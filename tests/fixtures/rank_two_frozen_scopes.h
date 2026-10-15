#pragma once
#include <array>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace rank_two_frozen {
struct Scope {
    std::string recordName;
    std::vector<int> states;
    std::size_t expectedRoots=0;
};

// Test scheduling only: the frozen geometry catalogue is not caller-extensible.
// A single-state record is never an all-state acceptance record.
inline Scope selectScope(const std::string& name) {
    constexpr std::array<int,8> ids{0,1,2,18,56,57,58,59};
    constexpr std::array<std::size_t,8> roots{45,48,46,46,47,42,42,44};
    if(name=="first")return {"state-0-cost-calibration",{0},45};
    if(name=="all")return {"eight-states",{ids.begin(),ids.end()},360};
    for(std::size_t i=0;i<ids.size();++i)
        if(name=="state-"+std::to_string(ids[i]))
            return {"single-frozen-state-"+std::to_string(ids[i]),{ids[i]},roots[i]};
    throw std::invalid_argument("unsupported or noncanonical frozen-state scope: "+name);
}
} // namespace rank_two_frozen
