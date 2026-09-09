#pragma once
#include <cstddef>
#include <map>
#include <memory>
#include <string>

namespace verification {
struct ExactNefBudgetFacts {
    std::string definition="sum of vertices+halfedges+halffacets+volumes of live explicit application Nef objects; includes shared auditor operands/results; excludes CGAL allocator internals and triangle-distance storage";
    std::size_t maximumComplexity=0,currentComplexity=0,peakComplexity=0;
    std::size_t preflightChecks=0,lastRequestedComplexity=0;
    bool exceeded=false;
    std::string lastOperation;
    std::map<std::string,std::size_t> observationsByOperation;
};
/** Opt-in context for one operation chain; no process-global state. Region
 * solids retain their charge until the last region owner releases them.
 * Do not share one mutable context between concurrent operation chains.
 * Consumers publishing immutable products must COPY the computed snapshot.
 */
class ExactNefBudget {
public:
    explicit ExactNefBudget(std::size_t maximumComplexity);
    ExactNefBudget(const ExactNefBudget&)=delete;
    ExactNefBudget& operator=(const ExactNefBudget&)=delete;
    const ExactNefBudgetFacts& facts() const {return facts_;}
private:
    ExactNefBudgetFacts facts_;
    void preflight(std::size_t additional,const std::string& operation);
    void acquire(std::size_t count,const std::string& operation);
    void release(std::size_t count) noexcept;
    friend class ExactNefBudgetLease;
};
using ExactNefBudgetPtr=std::shared_ptr<ExactNefBudget>;
std::string hashExactNefBudgetFacts(const ExactNefBudgetFacts& facts);
}
