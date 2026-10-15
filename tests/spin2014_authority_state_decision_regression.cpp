#include "spin2014/QuadricAuthorityNewtonState.h"
#include <stdexcept>
#include <iostream>
int main() {
 using namespace spin2014;
 spin::SpinSolverSettings settings;settings.deterministicSeed=0;
 const auto first=makeSpin2014AuthorityInitialSequence(spin::SpinMode::SpinningTop,settings);
 const auto repeat=makeSpin2014AuthorityInitialSequence(spin::SpinMode::SpinningTop,settings);
 if(first.size()!=8 || repeat.size()!=8)throw std::runtime_error("restart starts omitted");
 for(unsigned i=0;i<8;++i)if(first[i].geometry()!=repeat[i].geometry() || first[i].dual()!=repeat[i].dual())
   throw std::runtime_error("seed initialization not deterministic");
 settings.deterministicSeed=1;
 const auto other=makeSpin2014AuthorityInitialSequence(spin::SpinMode::SpinningTop,settings);
 if(first[0].geometry()==other[0].geometry())throw std::runtime_error("seed ignored");
 settings.maximumRestarts=9;
 bool rejectedSettings=false;
 try{(void)makeSpin2014AuthorityInitialSequence(spin::SpinMode::SpinningTop,settings);}
 catch(const std::invalid_argument&){rejectedSettings=true;}
 if(!rejectedSettings)throw std::runtime_error("unregistered restart budget allowed");
 QuadricAuthorityNewtonState state;
 state.multiplierMu="1/1";
 state.rawEquationResidual="100000000000000000000000000000000000001/100000000000000000000000000000000000000000000000000";
 for(auto& row:state.nativeMomentJacobian)row.fill("0/1");
 state.nativeMomentJacobian[0][0]="1/100000000000000000000000000000000000000000000000000";
 if(classifySpin2014AuthorityNumericalState(state,false)!="continue")
   throw std::runtime_error("above-threshold sub-double residual falsely converged");
 state.rawEquationResidual="1/1000000000000";
 if(classifySpin2014AuthorityNumericalState(state,false)!="stationary_candidate")
   throw std::runtime_error("threshold endpoint wrongly rejected");
 state.multiplierMu="-1/1";state.rawEquationResidual="1/1";
 if(classifySpin2014AuthorityNumericalState(state,true)!="initial_multiplier_nonpositive" ||
    classifySpin2014AuthorityNumericalState(state,false)!="continue")
   throw std::runtime_error("initial rule leaked into intermediate iterates");
 state.rawEquationResidual="0/1";
 if(classifySpin2014AuthorityNumericalState(state,false)!="stationary_multiplier_nonpositive_restart")
   throw std::runtime_error("nonminimum stationary state passed");
 state.multiplierMu="1/1";
 for(auto& row:state.nativeMomentJacobian)row.fill("0/1");
 if(classifySpin2014AuthorityNumericalState(state,true)!="initial_zero_sensitivity")
   throw std::runtime_error("zero sensitivity accepted");
 state.rawEquationResidual="-1/1";
 bool rejected=false;
 try {(void)classifySpin2014AuthorityNumericalState(state,false);}
 catch(const std::invalid_argument&){rejected=true;}
 if(!rejected)throw std::runtime_error("negative norm accepted");
 std::cout<<"exact authority initial/intermediate/stationary decisions passed\n";
}
