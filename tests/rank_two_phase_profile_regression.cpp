#include "rank_two_phase_profile.h"
#include "fixtures/rank_two_phase_catalogue_v1.h"
#include <iostream>
#include <sstream>
void require(bool v,const char* why){if(!v)throw std::runtime_error(why);}
int main(){try{
    using namespace rank_two_test;
    // Input observations with literal hand-derived accounting, not a fake
    // algorithm or replacement system clock used by the actual geometry.
    PhaseAccumulator a({100,50});a.advance("one",{100,50});
    a.advance("two",{110,55});a.advance("one",{140,65});a.advance("end",{145,68});
    const auto& one=a.totals().at("one");const auto& two=a.totals().at("two");
    require(one.spans==2&&one.wall==15&&one.cpu==8&&one.maxWall==10&&one.maxCpu==5,"lost repeated phase span");
    require(two.spans==1&&two.wall==30&&two.cpu==10,"incorrect separate phase accounting");
    auto snapshot=a;snapshot.advance("snapshot",{150,70});
    require(a.last().wall==145&&a.last().cpu==68,"snapshot modified live accumulator");
    std::uint64_t wall=0,cpu=0;for(const auto& p:snapshot.totals()){wall+=p.second.wall;cpu+=p.second.cpu;}
    require(wall==50&&cpu==20,"open snapshot interval omitted from totals");
    bool rejected=false;try{a.advance("bad",{144,68});}catch(const std::runtime_error&){rejected=true;}
    require(rejected,"nonmonotone clock accepted");
    require(phaseLabel("algebraic-root:validated-certificate id 2 1/1 2/1 huge")=="algebraic-root:validated-certificate","unbounded certificate used as label");
    require(phaseLabel("registry:subdivision:abc:coordinates-ready")=="unattributed:after-subdivision-coordinates-ready","completed operation mislabeled internal work");
    require(phaseLabel("regular-edge:attempt:abc:0:1/2:3/4")=="regular-edge:attempt","strip search hash generated separate histogram keys");
    require(phaseLabel("algebraic-root:reload-complete")=="unattributed:after-root-complete","root completion gap claimed as internal root time");
    require(phaseLabel("native-fiber:factor-root-set-100")=="native-fiber:factor-root-set"
        &&phaseLabel("native-fiber:factor-root-set-112")=="native-fiber:factor-root-set",
        "native root counts exhaust bounded phase label budget");
    // The actual whole trace reached the original 128-label cap when the
    // complete-product completion added a histogram name. Raw diagnostics
    // remain distinct; only AFTER-event timing categories share an identity.
    PhaseAccumulator bounded({0,0});std::uint64_t tick=0;
    bounded.advance(phaseLabel("native-fiber:raw-resultant-complete"),{tick,tick});
    for(int i=0;i<126;++i){++tick;bounded.advance("existing-"+std::to_string(i),{tick,tick});}
    ++tick;bounded.advance("observer:accounting",{tick,tick});
    ++tick;bounded.advance(phaseLabel("native-fiber:raw-resultant-full-product-verified"),{tick,tick});
    ++tick;bounded.advance(phaseLabel("native-fiber:raw-resultant-complete"),{tick,tick});
    ++tick;bounded.advance("observer:accounting",{tick,tick});
    require(bounded.totals().size()==128,"complete-product observation altered original label cap");
    require(phaseLabel("native-fiber:raw-resultant-full-product-direct-fallback")
        ==phaseLabel("native-fiber:raw-resultant-begin"),"product fallback loses direct-operation timing category");
    PhaseAccumulator componentBounded({0,0});tick=0;
    componentBounded.advance(phaseLabel("algebraic-root:native-factor-product-verified"),{tick,tick});
    for(int i=0;i<126;++i){++tick;componentBounded.advance("component-existing-"+std::to_string(i),{tick,tick});}
    ++tick;componentBounded.advance("observer:accounting",{tick,tick});
    ++tick;componentBounded.advance(phaseLabel("algebraic-root:resultant-component-product-verified"),{tick,tick});
    ++tick;componentBounded.advance("observer:accounting",{tick,tick});
    require(componentBounded.totals().size()==128,"component factor observation exhausts original label cap");
    PhaseAccumulator memoBounded({0,0});tick=0;
    memoBounded.advance(phaseLabel("regular-edge:attempt"),{tick,tick});
    for(int i=0;i<126;++i){++tick;memoBounded.advance("memo-existing-"+std::to_string(i),{tick,tick});}
    ++tick;memoBounded.advance("observer:accounting",{tick,tick});
    ++tick;memoBounded.advance(phaseLabel("regular-strip:memo-hit"),{tick,tick});
    ++tick;memoBounded.advance(phaseLabel("regular-strip:memo-miss"),{tick,tick});
    ++tick;memoBounded.advance("observer:accounting",{tick,tick});
    require(memoBounded.totals().size()==128,"strip reuse observation exhausts original label cap");
    // The actual new whole trace exhausted its original cap after adding
    // complete-component SF success/fallback. Keep raw events distinct but
    // group AFTER-event timing with the existing SF end/begin categories.
    PhaseAccumulator sfBounded({0,0});tick=0;
    sfBounded.advance(phaseLabel("native-fiber:square-free-certified"),{tick,tick});
    ++tick;sfBounded.advance(phaseLabel("native-fiber:square-free-begin"),{tick,tick});
    for(int i=0;i<125;++i){++tick;sfBounded.advance("sf-existing-"+std::to_string(i),{tick,tick});}
    ++tick;sfBounded.advance("observer:accounting",{tick,tick});
    ++tick;sfBounded.advance(phaseLabel("native-fiber:component-square-free-certified"),{tick,tick});
    ++tick;sfBounded.advance(phaseLabel("native-fiber:component-square-free-original-fallback"),{tick,tick});
    ++tick;sfBounded.advance("observer:accounting",{tick,tick});
    require(sfBounded.totals().size()==128,"complete component SF observation exhausts original label cap");
    require(sfBounded.totals().at("native-fiber:square-free-certified").spans==2
        &&sfBounded.totals().at("native-fiber:square-free-begin").spans==2,
        "new SF observations lost timing spans instead of grouping");
    std::cout<<"PASS complete-product timing categories retain original128label cap\n";
    std::cout<<"PASS literal wall/CPU accounting; snapshots, repeated spans and unknown gaps\n";
    // Break caught: stage() must not introduce an unbounded extra timing
    // identity after the actual first-face category set has reached its cap.
    // Populate ALL real saved identities, not synthetic existing-N filler.
    // Exercise the real stage path and system clocks used by the harness.
    std::ostringstream capture;PhaseProfile transition(capture);
    for(const auto* name:firstFacePhaseCatalogue)transition.stage(name);
    transition.stage("harness:validate-coverage-output");
    transition.stage("harness:import-and-validate");
    transition.finish();
    const auto text=capture.str();
    require(text.find("RANK_TWO_PHASE_ERROR")==std::string::npos,
        "actual first-face completion invalidates observation");
    std::istringstream records(text);std::string line;
    std::uint64_t beginWall=0,beginCpu=0,endWall=0,endCpu=0,sumWall=0,sumCpu=0,count=0;
    bool final=false,closed=false,hasValidation=false;
    while(std::getline(records,line)){
        std::istringstream row(line);std::string kind;row>>kind;
        if(kind=="RANK_TWO_PHASE_SNAPSHOT_BEGIN"){
            int version;std::uint64_t seq,events;std::string reason;
            row>>version>>seq>>events>>beginWall>>beginCpu>>endWall>>endCpu>>reason;
            require(bool(row)&&version==1,"invalid actual snapshot format");
            final=reason=="final";closed=false;sumWall=0;sumCpu=0;count=0;
        }else if(final&&kind=="RANK_TWO_PHASE_BUCKET"){
            std::uint64_t seq,spans,w,c,mw,mc;std::string phase;
            row>>seq>>phase>>spans>>w>>c>>mw>>mc;
            require(bool(row)&&spans>0,"dropped actual timing span");
            sumWall+=w;sumCpu+=c;++count;
            if(phase=="harness:import-and-validate")hasValidation=true;
        }else if(final&&kind=="RANK_TWO_PHASE_SNAPSHOT_END")closed=true;
    }
    require(final&&closed&&count==128&&hasValidation,
        "first-face validation changed original bounded catalogue");
    require(sumWall==endWall-beginWall&&sumCpu==endCpu-beginCpu,
        "first-face stage grouping dropped wall/CPU time");
    // Break caught: diagnostic replays AFTER four complete source owners
    // add a 129th phase and invalidate the real observer, although source
    // construction already finished. Actual R179 snapshot8557 contains this
    // same 128-name catalogue; sampleRankTwoArc's first event is the literal
    // endpoint-reload below. Exercise actual stage/event/clocks, not a mock.
    std::ostringstream replayCapture;PhaseProfile replay(replayCapture);
    for(const auto* name:firstFacePhaseCatalogue)replay.stage(name);
    replay.stage("harness:source-projection-replay");
    replay.event("arc-sample:endpoint-reload");
    replay.stage("harness:source-subdivision-replay");
    replay.event("registry:validate-plan");
    replay.finish();
    const auto replayText=replayCapture.str();
    require(replayText.find("RANK_TWO_PHASE_ERROR")==std::string::npos,
        "post-construction replay invalidates bounded observer");
    require(replayText.find("WHOLE_COVERAGE_DIAGNOSTIC arc-sample:endpoint-reload\n")!=std::string::npos
        &&replayText.find("WHOLE_COVERAGE_DIAGNOSTIC registry:validate-plan\n")!=std::string::npos,
        "diagnostic stage grouping discarded raw native events");
    std::istringstream replayRecords(replayText);
    final=false;closed=false;hasValidation=false;
    while(std::getline(replayRecords,line)){
        std::istringstream row(line);std::string kind;row>>kind;
        if(kind=="RANK_TWO_PHASE_SNAPSHOT_BEGIN"){
            int version;std::uint64_t seq,events;std::string reason;
            row>>version>>seq>>events>>beginWall>>beginCpu>>endWall>>endCpu>>reason;
            require(bool(row)&&version==1,"invalid replay snapshot format");
            final=reason=="final";closed=false;sumWall=0;sumCpu=0;count=0;
        }else if(final&&kind=="RANK_TWO_PHASE_BUCKET"){
            std::uint64_t seq,spans,w,c,mw,mc;std::string phase;
            row>>seq>>phase>>spans>>w>>c>>mw>>mc;
            require(bool(row)&&spans>0,"dropped diagnostic replay timing span");
            sumWall+=w;sumCpu+=c;++count;
            if(phase=="harness:import-and-validate")hasValidation=true;
        }else if(final&&kind=="RANK_TWO_PHASE_SNAPSHOT_END")closed=true;
    }
    require(final&&closed&&count==128&&hasValidation,
        "post-construction replay changed original bounded catalogue");
    require(sumWall==endWall-beginWall&&sumCpu==endCpu-beginCpu,
        "diagnostic replay grouping dropped wall/CPU time");
    std::cout<<"PASS actual128category post-construction replay; raw events retained; time conserved\n";
    // Genuine distinct-category overflow must still fail; do not weaken cap.
    PhaseAccumulator cap({0,0});tick=0;
    for(const auto* name:firstFacePhaseCatalogue){cap.advance(name,{tick,tick});++tick;}
    cap.advance("genuinely-new-category",{tick,tick});++tick;
    rejected=false;try{cap.advance("observer:accounting",{tick,tick});}
    catch(const std::runtime_error&){rejected=true;}
    require(rejected,"real distinct-category overflow no longer rejected");
    std::cout<<"PASS actual128category first-face stage transition; time conserved; overflow rejected\n";
    PhaseProfile real(std::cout);real.event("algebraic-root:certify");real.event("algebraic-root:reload-complete");real.finish();
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
