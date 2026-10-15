#include "detail/RankTwoFinitePartition.h"
#include "detail/RankTwoChordSweep.h"
#include "detail/RankTwoRegistryRotation.h"

namespace spin {
namespace {
void requirePartition(const RankTwoArrangement& a,const RankTwoFinitePartition& p){
    if(p.arrangementIdentity()!=a.scientificHash()
        ||p.registry().arrangementIdentity()!=a.scientificHash()
        ||p.registryIdentity()!=p.registry().scientificHash())
        throw RankTwoSourceError("INPUT_IDENTITY_MISMATCH","finite mesh scope requires the actual immutable native partition");
}
}
RankTwoChordIncidence certifyRankTwoChordIncidenceSweep(const RankTwoArrangement& a,
    const RankTwoFinitePartition& p,const RankTwoChordSweepSettings& s){
    requirePartition(a,p);
    return RankTwoChordIncidence::sweep(a,p.registry(),p.scientificHash(),p.pointIndices(),p.edgeIndices(),s);
}
RankTwoRegistryRotation certifyRankTwoRegistryRotation(const RankTwoArrangement& a,
    const RankTwoFinitePartition& p,const RankTwoChordIncidence& i,const RankTwoRegistryRotationSettings& s){
    requirePartition(a,p);
    return RankTwoRegistryRotation::certify(a,p.registry(),i,p.scientificHash(),p.pointIndices(),p.edgeIndices(),s);
}
} // namespace spin
