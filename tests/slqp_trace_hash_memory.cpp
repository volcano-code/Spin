// Synthetic trace storage benchmark, never a paper optimization case.
#include "spin2014/SlqpAuditTrace.h"
#include <iostream>
#include <sys/resource.h>
int main(int argc,char** argv){try {
    if(argc>1 && std::string(argv[1])=="--limited") {
        rlimit limit{192ULL*1024*1024,192ULL*1024*1024};
        if(setrlimit(RLIMIT_AS,&limit)!=0){std::cerr<<"cannot enforce memory test cap\n";return 1;}
    }
    spin2014::SlqpLbfgsResult r;
    r.status=spin2014::SlqpStatus::MaximumIterations;r.diagnostic="retained full synthetic trace";
    r.trace.resize(64);
    for(int i=0;i<64;++i){
        auto& t=r.trace[i];t.iteration=i;t.objective=1./(i+1);
        t.state.resize(100000);
        for(int j=0;j<t.state.size();++j)t.state(j)=double(j+1)/100001.;
    }
    spin2014::SlqpLbfgsSettings s;
    s.algorithm=spin2014::SlqpAlgorithm::AuditedKktRestorationV1;
    const auto h=spin2014::hashAuditedSlqp(r,s);
    if(h!="f6fe10c20e2aec7feaa90232ff30e76a15ca1eeb14aa96fa635784959a27b03c") {
        std::cerr<<"historical retained trace digest changed\n";return 1;
    }
    rusage u{};getrusage(RUSAGE_SELF,&u);
    auto rss=u.ru_maxrss;
#ifdef __APPLE__
    rss/=1024;
#endif
    std::cout<<h<<" peak_RSS_KiB="<<rss<<" retained_doubles=6400000\n";
    return 0;
}catch(const std::bad_alloc&){std::cerr<<"HASH_MEMORY_LIMIT: retained trace plus hash temporary exceeded memory\n";return 2;}}
