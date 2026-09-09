# Independent controller small-problem derivations

These are analytic correctness fixtures, not original paper experiments, author trajectories, or global nonlinear feasibility proofs. No production checker/solver functions used in derivations.

Convention: r=gradf+J^Tlambda-alpha+eta, bounded prefix0<=x<=1; allbound multipliers on unboundedcoordinates zero. Component∞norms are equality, bounds, stationarity, dualsign, and complementarity. FullKKT is maximum of these, not projectedgradient alone.

| Fixed case | Independent expectation |
|---|---|
| f=x, x>=0, x=0 | alpha=1,eta=0 gives r=0 andcomplementarity0. |
| f=-x, x<=1, x=1 | alpha=0,eta=1 gives r=0 andcomplementarity0. |
| f=x, c=x-.75, solutionx=.75 | lambda=-1,alpha=eta=0; allcomponents0. |
| Same, initialx=0 | lambda=-1 can give stationarity0, but equality=.75 so fullKKT=.75, neverconverged. Existingworkingboundmask mayfreeze x despite requiredpositive restorationdirection. |
| Duplicateconstraints(x-.75,2x-1.5) | x=.75, anylambda1+2lambda2=-1 givesvalidstationarity. Do not compare nonunique individualmultipliers. |
| Inconsistentconstraints(x-.75,x-.25) | triangleinequality implies |x-.75|+|x-.25|>=.5, attainedforallx in[.25,.75]. A localLPcan reduceinitialL1=1 but cannot makebothconstraintszero. Nonzeroequationresidual mustsurvive evenwithzeroobjective stationarity. |
| Singular nonlinear c=x²-.25,start0 | J=0, so everylinearizedstep leavespredictedconstraint-.25. No predictedrestoration improvement. Nonlinear feasiblex=.5 exists, so stall is a particularlyclear counterexample to solverfailure=>infeasibility. |
| f=.5[(x-.8)²+(y+.2)²], box[0,1]² | minimizer(.8,0),gradient(0,.2),alpha=(0,.2),eta=(0,0). Strictconvexity confirmsuniqueglobalquadraticminimum, not a property grantedbygeneralKKTchecker. |

RestorationLP for singlelinearcase atx0 withDelta1: d in[0,1], min|d-.75|, uniquezero-residualstep d=.75. Objective increases from0 to.75 but actualconstraintL1 fallsfrom.75to0, so feasibilityrestoration mustnotbe rejected solely becausef increases.

ForKKT minimaxwitness, let R=gradf+J^Tlambda. Lowerfacealpha=max(R,0): stationarityresidual=min(R,0), henceR>=-t suffices. Upperfaceeta=max(-R,0): stationarityresidual=max(R,0), henceR<=t suffices. Interior/unbounded mustsatisfy|R|<=t. Nearlowerboundx>=0 complementarity=x*max(R,0); imposingx*R<=t controlsit wheneverR>0, automaticallytrueifR<0. Nearupperbound complementarity=(1-x)*max(-R,0), controlledby-(1-x)*R<=t. These relations assumeactualboxfeasibility, whichmustbe checked separately; no freeboundmultipliers onphi. Derivedvectors must be checked independently afterLPstatus.

First-orderKKT is not second-order/local-minimum proof: f=-(x-.5)² atx=.5 haszeroKKT but is a maximum. Avoid any localminimum/globaloptimal/physicalstable language for a generic KKT convergence result.
