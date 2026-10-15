#pragma once
#include "RankTwoDiagnostics.h"
#include <CGAL/Cartesian.h>
#include <CGAL/CORE_Expr.h>
#include <array>
#include <functional>
#include <map>
#include <optional>
#include <sstream>

namespace spin::detail {
// Cartesian<CORE::Expr> predicates with a certified structural identity case
// and optional observation. No rounded coordinates or proximity threshold.
// All nonidentical-point triples still use the original CGAL functor.
struct RankTwoPredicateObservation {
    using Base=CGAL::Cartesian<CORE::Expr>;
    using Key=std::array<const CORE::ExprRep*,2>;
    struct Less {
        bool operator()(const Key& a,const Key& b)const noexcept {
            const std::less<const CORE::ExprRep*> less;
            return less(a[0],b[0]) || (a[0]==b[0] && less(a[1],b[1]));
        }
    };
    std::map<Key,std::size_t,Less> points;
    std::optional<std::size_t> constraint;
    void bind(const Base::Point_2& point,std::size_t index){
        points.emplace(Key{point.x().Rep(),point.y().Rep()},index);
    }
    std::string index(const Base::Point_2& point)const {
        const auto found=points.find(Key{point.x().Rep(),point.y().Rep()});
        return found==points.end()?"UNREGISTERED":std::to_string(found->second);
    }
    void emit(bool begin,const Base::Point_2& a,const Base::Point_2& b,
              const Base::Point_2& c,int result=0)const {
        if(!constraint || !rankTwoDiagnosticTarget.sink)return;
        std::ostringstream record;record<<"finite-cdt:predicate-"<<(begin?"begin":"end")
            <<" ORIENTATION2 "<<*constraint<<' '<<index(a)<<' '<<index(b)<<' '<<index(c);
        if(!begin)record<<' '<<result;
        rankTwoDiagnostic(record.str().c_str());
    }
};
struct RankTwoObservedCartesian : CGAL::Cartesian<CORE::Expr> {
    using Base=CGAL::Cartesian<CORE::Expr>;
    RankTwoPredicateObservation* observation=nullptr;
    explicit RankTwoObservedCartesian(RankTwoPredicateObservation* state=nullptr):observation(state){}
    struct Orientation_2 : Base::Orientation_2 {
        using Base::Orientation_2::operator();
        RankTwoPredicateObservation* observation=nullptr;
        explicit Orientation_2(RankTwoPredicateObservation* state):observation(state){}
        static bool sameExpressionPoint(const Base::Point_2& a,const Base::Point_2& b)noexcept {
            return a.x().Rep()==b.x().Rep() && a.y().Rep()==b.y().Rep();
        }
        result_type operator()(const Base::Point_2& a,const Base::Point_2& b,const Base::Point_2& c)const {
            if(observation)observation->emit(true,a,b,c);
            // The ACTUAL constraint23 replay repeats its second point. CORE
            // otherwise builds two equal high-degree determinant DAGs whose
            // zero separation bound is prohibitively costly. Identical x/y
            // expression representations prove point equality exactly. This
            // covers all duplicate placements, never nearby distinct points;
            // different representations of an equal root still fall through.
            const auto result=(sameExpressionPoint(a,b) || sameExpressionPoint(a,c) || sameExpressionPoint(b,c))
                ?CGAL::COLLINEAR:Base::Orientation_2::operator()(a,b,c);
            if(observation)observation->emit(false,a,b,c,static_cast<int>(result));
            return result;
        }
    };
    Orientation_2 orientation_2_object()const{return Orientation_2(observation);}
};
} // namespace spin::detail
