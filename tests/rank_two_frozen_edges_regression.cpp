#include "spin/detail/RankTwoWholeArrangement.h"
#include "spin/detail/RankTwoDiagnostics.h"
#include "verification/Sha256.h"
#if __has_include("fixtures/rank_two_frozen_scopes.h")
#include "fixtures/rank_two_frozen_scopes.h"
#endif
#include <algorithm>
#include <chrono>
#include <cstring>
#include <cmath>
#include <fstream>
#include <iostream>
#include <locale>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>

namespace {
using Q=CGAL::Gmpq;
using E=CORE::Expr;
using Edge=std::array<std::size_t,2>;
void require(bool yes,const std::string& message){if(!yes)throw std::runtime_error(message);}
void scopeTests(const char* output){
#if __has_include("fixtures/rank_two_frozen_scopes.h")
    using rank_two_frozen::selectScope;
    const std::vector<int> ids{0,1,2,18,56,57,58,59};
    const std::vector<std::size_t> roots{45,48,46,46,47,42,42,44};
    const auto first=selectScope("first"),all=selectScope("all");
    require(first.states==std::vector<int>{0} && first.expectedRoots==45
        && first.recordName=="state-0-cost-calibration","legacy first scope changed");
    require(all.states==ids && all.expectedRoots==360 && all.recordName=="eight-states",
        "all scope omits or duplicates a required frozen state");
    std::set<int> seen;
    for(std::size_t i=0;i<ids.size();++i){
        const auto one=selectScope("state-"+std::to_string(ids[i]));
        require(one.states==std::vector<int>{ids[i]} && one.expectedRoots==roots[i],
            "single-state selector silently selected another case");
        require(seen.insert(one.states.front()).second,"single-state selectors overlap");
    }
    require(std::vector<int>(seen.begin(),seen.end())==ids,"single-state schedule is incomplete");
    for(const std::string bad:{"","state-3","state-00","state--1","state-1junk",
            "STATE-1","first,all","all,all","stages","scopes"}){
        bool refused=false;try{selectScope(bad);}catch(const std::invalid_argument&){refused=true;}
        require(refused,"malformed or unsupported scope was accepted");
    }
    std::ifstream prior(output);require(!prior.good(),"refuse overwrite scope contract");
    std::ofstream out(output);
    out<<"rank-two-frozen-scope-contract-v1\nrequired_states 0 1 2 18 56 57 58 59\n"
          "required_edges 3840\nrequired_roots 360\ngeometry_evaluation NOT_EXECUTED\n";
    out.close();require(bool(out),"cannot save scope contract");
    std::cout<<"PASS explicit frozen scopes: complete catalogue, no aliases or silent fallback\n";
#else
    throw std::runtime_error("MISSING_FEATURE: explicit frozen-state scopes");
#endif
}
std::string checkedFile(const char* path,const char* sha){
    std::ifstream in(path,std::ios::binary);require(bool(in),"missing frozen input");
    const std::string bytes((std::istreambuf_iterator<char>(in)),std::istreambuf_iterator<char>());
    verification::Sha256 hash;hash.update(reinterpret_cast<const unsigned char*>(bytes.data()),bytes.size());
    require(hash.finish()==sha,"frozen input identity mismatch");return bytes;
}
E exact(const Q& q){
    std::ostringstream a,b;a<<q.numerator();b<<q.denominator();
    return E(CORE::BigRat(CORE::BigInt(a.str()),CORE::BigInt(b.str())));
}
void end(std::istringstream& row){std::string extra;require(!(row>>extra),"extra frozen record tokens");}
struct Expected {std::vector<std::array<Q,2>> roots;std::vector<bool> inside;};
struct State {spin::RankTwoExactCoefficients q;std::map<Edge,Expected> edges;};
std::map<int,State> readStates(const std::string& input){
    std::istringstream lines(input);lines.imbue(std::locale::classic());
    std::string line,kind;std::map<int,State> result;State* state=nullptr;
    while(std::getline(lines,line)){
        std::istringstream row(line);row.imbue(std::locale::classic());require(bool(row>>kind),"empty frozen state record");
        if(kind=="state"){
            int id;require(bool(row>>id),"missing state ID");end(row);
            require(result.emplace(id,State{}).second,"duplicate state");state=&result.at(id);
            require(bool(std::getline(lines,line)),"missing coefficients");
            std::istringstream coefficients(line);coefficients.imbue(std::locale::classic());
            require(bool(coefficients>>kind) && kind=="q","missing q record");
            for(auto& q:state->q){double value;require(bool(coefficients>>value) && std::isfinite(value),"invalid q");q=Q(value);}
            end(coefficients); // stored binary64 authority, no renormalization
        }else if(kind=="edge"){
            require(state!=nullptr,"edge before state");Edge edge;int certified;std::size_t roots,inside;
            require(bool(row>>edge[0]>>edge[1]>>certified>>roots>>inside),"bad edge record");
            require(edge[0]<edge[1] && certified==1 && roots<=4 && inside==roots+1,"uncertified/malformed old edge");
            Expected expected;
            for(std::size_t i=0;i<inside;++i){int value;require(bool(row>>value) && (value==0 || value==1),"bad inside state");expected.inside.push_back(value);}
            end(row);
            for(std::size_t i=0;i<roots;++i){
                require(bool(std::getline(lines,line)),"missing old root");std::istringstream root(line);
                std::array<Q,2> bounds;require(bool(root>>kind>>bounds[0]>>bounds[1]) && kind=="root","bad root record");end(root);
                require(bounds[0]>=0 && bounds[0]<bounds[1] && bounds[1]<=1,"bad root bracket");expected.roots.push_back(bounds);
            }
            require(state->edges.emplace(edge,std::move(expected)).second,"duplicate old edge");
        }else throw std::runtime_error("unexpected old certificate record");
    }
    const std::map<int,std::size_t> counts{{0,45},{1,48},{2,46},{18,46},{56,47},{57,42},{58,42},{59,44}};
    require(result.size()==counts.size(),"incomplete state catalogue");
    for(const auto& entry:counts){
        const auto& state=result.at(entry.first);require(state.edges.size()==480,"incomplete old edge set");
        std::size_t n=0;for(const auto& edge:state.edges)n+=edge.second.roots.size();require(n==entry.second,"wrong old root total");
    }
    return result;
}
spin::PreparedRankTwoNormalMesh readNormals(const std::string& input){
    std::istringstream lines(input);std::string line;spin::PreparedRankTwoNormalMesh::Rays rays;
    spin::PreparedRankTwoNormalMesh::Faces faces;
    while(std::getline(lines,line)){
        std::istringstream row(line);row.imbue(std::locale::classic());char kind;require(bool(row>>kind),"empty ray record");
        if(kind=='v'){
            std::array<Q,3> ray;for(auto& q:ray){double d;require(bool(row>>d) && std::isfinite(d),"bad normal");q=Q(d);}rays.push_back(ray);
        }else if(kind=='f'){
            std::array<std::size_t,3> f;require(bool(row>>f[0]>>f[1]>>f[2]),"bad face");faces.push_back(f);
        }else throw std::runtime_error("unknown normal record");end(row);
    }
    require(rays.size()==162 && faces.size()==320,"wrong fixed grid");
    return spin::PreparedRankTwoNormalMesh::fromExact(rays,faces);
}
// Diagnostic only: no whole-atlas or frozen-edge acceptance is inferred from
// a completed phase. Timing is kept outside each scientific canonical record.
void stages(const char* normalPath,const char* oldPath,const char* output){
    std::ifstream prior(output);require(!prior.good(),"refuse overwrite stage diagnostic");
    std::ofstream saved(output);require(bool(saved),"cannot create stage diagnostic");
    saved<<"rank-two-frozen-stage-diagnostic-v1\nstate 0 radius 4\n"
           "selected_faces 0 1 159 319\nwhole_atlas NOT_CHECKED\n"<<std::flush;
    using Clock=std::chrono::steady_clock;
    struct Trace {std::ofstream* output;Clock::time_point origin;bool failed=false;};
    Trace trace{&saved,Clock::now()};
    const auto observe=[](const char* phase,void* context) noexcept {
        auto& t=*static_cast<Trace*>(context);
        try {
            const auto ns=std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now()-t.origin).count();
            *t.output<<"TRACE "<<ns<<' '<<phase<<'\n'<<std::flush;
            std::cout<<"TRACE "<<ns<<' '<<phase<<'\n'<<std::flush;
            if(!*t.output || !std::cout)t.failed=true;
        }catch(...){t.failed=true;}
    };
    spin::detail::ScopedRankTwoDiagnostics observer(observe,&trace);
    auto begin=[&](const std::string& name){
        saved<<"BEGIN "<<name<<'\n'<<std::flush;
        std::cout<<"BEGIN "<<name<<'\n'<<std::flush;return Clock::now();
    };
    auto finish=[&](const std::string& name,Clock::time_point start){
        const auto ns=std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now()-start).count();
        saved<<"END "<<name<<" nanoseconds "<<ns<<'\n'<<std::flush;
        std::cout<<"END "<<name<<" nanoseconds "<<ns<<'\n'<<std::flush;
        require(bool(saved),"stage diagnostic persistence failure");
    };
    auto start=begin("read-frozen-states");
    const auto states=readStates(checkedFile(oldPath,"9d7c823fc966a4d91016ca549d04a245d617427d57ebda89566a0191ab2eca5f"));
    finish("read-frozen-states",start);
    start=begin("source-factory");
    const auto source=spin::PreparedRankTwoSource::fromExact(states.at(0).q,"frozen-rank-two-target-v1");
    finish("source-factory",start);saved<<"input "<<source.inputHash()<<'\n'<<std::flush;
    start=begin("normal-mesh-admission");
    const auto normals=readNormals(checkedFile(normalPath,"5ad29232268675e93c7f8be0de0709fefbb242310ed7ca2cfc7797b2ada21140"));
    finish("normal-mesh-admission",start);
    start=begin("global-arrangement");
    const auto global=spin::analyzeRankTwoGlobalArrangement(source,Q(4));
    finish("global-arrangement",start);
    saved<<"global_hash "<<global.scientificHash()<<'\n'<<global.canonicalRecord()<<std::flush;
    for(const std::size_t face:{0U,1U,159U,319U}){
        const auto name="normal-face-"+std::to_string(face);start=begin(name);
        const auto local=spin::analyzeRankTwoNormalFace(source,global,normals.face(face));
        finish(name,start);
        saved<<"face_hash "<<face<<' '<<local.scientificHash()<<'\n'<<local.canonicalRecord()<<std::flush;
    }
    saved<<"source_material_derivative_completion NOT_GRANTED\n";
    require(!trace.failed,"phase trace persistence failed");
    saved.close();require(bool(saved),"cannot finalize stage diagnostic");
}
std::size_t checkEdge(const spin::PreparedRankTwoSource& source,
    const spin::PreparedRankTwoNormalMesh& normals,
    const spin::RankTwoSharedEdgeRecord& actual,const Expected& old){
    std::vector<E> roots;
    for(const auto& event:actual.events){
        const bool radius=std::find(event.zeroPolynomials.begin(),event.zeroPolynomials.end(),"F")!=event.zeroPolynomials.end();
        require(radius==(event.radiusOrdinal>=0),"inconsistent radius event identity");
        if(!radius)continue;
        require(event.edge==actual.edge && event.radiusOrdinal==static_cast<std::ptrdiff_t>(roots.size()),"wrong radius ordinal");
        require(event.parameter>0 && event.parameter<1 && (roots.empty() || roots.back()<event.parameter),"noncanonical root ordering");
        roots.push_back(event.parameter);
    }
    require(roots.size()==old.roots.size(),"old/new radius root count mismatch");
    for(std::size_t i=0;i<roots.size();++i)
        require(exact(old.roots[i][0])<roots[i] && roots[i]<exact(old.roots[i][1]),"new exact root outside old isolator");
    Q u0=0,u1=0;
    for(int d=0;d<3;++d){u0+=source.graph().kernel[d]*normals.rays()[actual.edge[0]][d];u1+=source.graph().kernel[d]*normals.rays()[actual.edge[1]][d];}
    std::vector<E> points={E(0),E(1)};points.insert(points.end(),roots.begin(),roots.end());
    for(const auto& part:actual.intervals){points.push_back(part.lowerParameter);points.push_back(part.upperParameter);}
    if(u0!=u1){const Q ideal=u0/(u0-u1);if(ideal>0 && ideal<1)points.push_back(exact(ideal));}
    std::sort(points.begin(),points.end());points.erase(std::unique(points.begin(),points.end()),points.end());
    for(std::size_t i=1;i<points.size();++i){
        const E& left=points[i-1];const E& right=points[i];require(left>=0 && right<=1,"interval outside normal edge");
        // The ideal breakpoint is included. Signs at interval boundaries
        // therefore determine the complete open interval, without sampling.
        const int sign0=(exact(u0)+(exact(u1-u0)*left)).sign();
        const int sign1=(exact(u0)+(exact(u1-u0)*right)).sign();
        require(sign0*sign1>=0,"missed ideal breakpoint");const bool positive=sign0>0 || sign1>0;
        std::size_t coverage=0;bool inside=false;
        for(const auto& part:actual.intervals)
            if(part.lowerParameter<=left && right<=part.upperParameter){++coverage;inside=part.originalSigns[0]<0;}
        require(coverage==(positive?1U:0U),"missing/overlapping positive normal arc");
        std::size_t bin=0;while(bin<roots.size() && roots[bin]<=left)++bin;
        require(inside==old.inside.at(bin),"old/new radius inside interval mismatch");
    }
    return roots.size();
}
void run(const char* normalPath,const char* oldPath,const char* output,const std::string& scopeName){
    const auto scope=rank_two_frozen::selectScope(scopeName);
    const auto normals=readNormals(checkedFile(normalPath,"5ad29232268675e93c7f8be0de0709fefbb242310ed7ca2cfc7797b2ada21140"));
    const auto states=readStates(checkedFile(oldPath,"9d7c823fc966a4d91016ca549d04a245d617427d57ebda89566a0191ab2eca5f"));
    std::ifstream prior(output);require(!prior.good(),"refuse overwrite crosscheck");std::ofstream saved(output);require(bool(saved),"cannot create crosscheck");
    saved<<"rank-two-frozen-edge-crosscheck-v1\nscope "<<scope.recordName<<'\n'<<std::flush;
    std::size_t totalEdges=0,totalRoots=0;
    for(const auto& state:states){
        if(std::find(scope.states.begin(),scope.states.end(),state.first)==scope.states.end())continue;
        std::cout<<"BEGIN frozen state "<<state.first<<" 320 actual overlays\n"<<std::flush;
        const auto source=spin::PreparedRankTwoSource::fromExact(state.second.q,"frozen-rank-two-target-v1");
        saved<<"state_begin "<<state.first<<" input "<<source.inputHash()<<'\n'<<std::flush;
        struct Progress {int state;std::size_t faces=0,edges=0;bool failed=false;
            std::chrono::steady_clock::time_point start=std::chrono::steady_clock::now();};
        Progress progress{state.first};
        const auto observe=[](const char* phase,void* context) noexcept {
            auto& p=*static_cast<Progress*>(context);
            const bool face=std::strcmp(phase,"overlay:end")==0;
            const bool edge=std::strcmp(phase,"atlas-shared-edge:end")==0;
            if(face)++p.faces;if(edge)++p.edges;
            if((!face || p.faces%20) && (!edge || p.edges%20)
                && std::strcmp(phase,"atlas-shared-edges:begin")!=0
                && std::strcmp(phase,"atlas-shared-edges:end")!=0)return;
            try {
                const auto ns=std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::steady_clock::now()-p.start).count();
                std::cout<<"PROGRESS state "<<p.state<<" faces "<<p.faces<<" edges "<<p.edges
                    <<" elapsed_ns "<<ns<<" phase "<<phase<<'\n'<<std::flush;
                if(!std::cout)p.failed=true;
            }catch(...){p.failed=true;}
        };
        spin::detail::ScopedRankTwoDiagnostics observer(observe,&progress);
        const auto result=spin::analyzeRankTwoArrangement(source,Q(4),normals);
        require(!progress.failed && progress.faces==320 && progress.edges==480,"incomplete whole-atlas diagnostic receipt");
        const auto& atlas=spin::wholeArrangementRecords(result);
        require(atlas.sharedEdges.size()==480 && atlas.normalFaces.size()==320,"incomplete new normal atlas");
        std::size_t roots=0;std::set<Edge> seen;
        for(const auto& edge:atlas.sharedEdges){
            require(seen.insert(edge.edge).second,"duplicate new edge");
            try{roots+=checkEdge(source,normals,edge,state.second.edges.at(edge.edge));}
            catch(const std::exception& e){throw std::runtime_error("state "+std::to_string(state.first)+" edge "+std::to_string(edge.edge[0])+" "+std::to_string(edge.edge[1])+": "+e.what());}
        }
        // Verify a wrong inside label cannot pass merely because root counts match.
        auto wrong=state.second.edges.at(atlas.sharedEdges.front().edge);wrong.inside[0]=!wrong.inside[0];
        bool rejected=false;try{checkEdge(source,normals,atlas.sharedEdges.front(),wrong);}catch(const std::runtime_error&){rejected=true;}
        require(rejected,"checker accepted changed old inside state");
        saved<<"state "<<state.first<<" edges 480 roots "<<roots<<" hash "<<result.scientificHash()<<'\n'<<result.canonicalRecord()<<std::flush;
        require(bool(saved),"crosscheck persistence failure");totalEdges+=480;totalRoots+=roots;
        std::cout<<"PASS frozen state "<<state.first<<" edges 480 roots "<<roots<<'\n'<<std::flush;
    }
    require(totalEdges==480*scope.states.size() && totalRoots==scope.expectedRoots,"wrong checked scope");
    saved<<"checked_edges "<<totalEdges<<" checked_roots "<<totalRoots<<"\nsource_material_derivative_completion NOT_GRANTED\n";
    saved.close();require(bool(saved),"cannot finalize crosscheck");
}
}
int main(int argc,char** argv){try{
    if(argc==5 && std::string(argv[4])=="scopes"){scopeTests(argv[3]);return 0;}
    if(argc!=5)
        throw std::invalid_argument("normals, old certificate, new output, explicit frozen scope required");
    if(std::string(argv[4])=="stages"){stages(argv[1],argv[2],argv[3]);return 0;}
    run(argv[1],argv[2],argv[3],argv[4]);return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
