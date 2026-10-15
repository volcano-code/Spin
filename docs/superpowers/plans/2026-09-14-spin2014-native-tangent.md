# Native metric-tangent rescue for Spin-It 2014

User-approved next step from docs/spin2014_direction_audit_2026-09-14.md.

## Global constraints

Preserve dirty work and all failed runs. No commit/tag/push/cleanup. No Faster
edits or edits in the main worktree. Keep old default algorithm and bytes.
New rescue is explicitly IMPLEMENTATION_DEFAULT, not author-exact. Do not alter
objectives, shell, pose, rounding, KKT or binary tolerances. Compile <=4 jobs;
large experiments sequential, <=24 GiB address space and <=3600 seconds each.

### Task 1: Native tangent QP kernel

Own only include/spin2014/SlqpTangentQp.h, src/spin2014/SlqpTangentQp.cpp,
tests/spin_it_2014_tangent_qp.cpp. Do not edit CMake or solver/history/bindings.
Use TDD (missing-header red permitted), standalone compile with SlqpHistory.cpp
and Eigen. Prototype reference experiments/teapot_spin2014/direction_replay.py.

Exact API in namespace spin2014:
struct TangentQpSettings { int maximumFaces=256; int maximumIterations=1028;
 double activeTolerance=1e-10; double residualTolerance=1e-9; };
struct TangentQpStep { std::vector<int> workingSet; std::string action;
 int changedIndex=-1; double equalityResidual=0.; double stationarity=0.; };
struct TangentQpResult { bool attempted=false; bool checked=false;
 std::string status="NOT_ATTEMPTED"; Eigen::VectorXd direction;
 std::vector<int> workingSet; Eigen::VectorXd multipliers;
 double equalityResidual=0.; double stationarity=0.;
 std::vector<TangentQpStep> trace; };
TangentQpResult solveTangentQp(const Eigen::VectorXd& x,const Eigen::VectorXd& g,
 const Eigen::MatrixXd& a,const Eigen::VectorXd& c,Eigen::Index bounded,
 const LbfgsHistory& history,const TangentQpSettings& settings={});

Input validation: finite/dimensions/box prefix/settings/history consistency.
Use the existing inverse L-BFGS action; independent direct BFGS Hessian recursion
for stationarity. Never dense n*n Hessian. Working set: initialize boundary
faces, add outward inactive faces, release wrong-sign multipliers; deterministic
ties; record every working set/action. Fixed direction entries exactly zero,
then recheck equalities and stationarity. Failure/cycle/rank/resource are UNKNOWN,
not infeasibility. Face/iteration limits before allocation and loops. No silent
repair of free variables. Bounded coordinates only, unbounded Givens coordinates
are not faces. Status checked only with finite verified residuals and face signs.
QP multipliers order: original equalities followed by returned workingSet rows.

Tests: analytic H=[[1,2],[2,9]] pair s=[1,2],y=[1,0], lower-bound solution
d=[0,-5]; upper counterpart; release wrong identity face; equality gives
[5/7,-5/7]; inconsistent tangent equality UNKNOWN; fixed-face roundoff; invalid
input/history; face/iteration cap; determinism. Do not claim whole-spin success.
No subagents/commits. Report test commands and red/green evidence.

### Task 2: Controller-owned integration and provenance

Add explicit opt-in direction policy METRIC_TANGENT_RESCUE_V1 alongside unchanged
default. Only try metric QP when old QN direction is blocked by a <=1e-14 bound
step and equality is already within existing tolerance. Unknown QP or failed
merit retains logged existing first-order fallback. Accepted rescue keeps history;
still use actual nonlinear merit and finite box step. Record/hash policy, budgets,
QP transitions/final direction/multipliers/residuals and rescue selection, separate
from old projected membership. Old controls/trace bytes remain unchanged when off.
Expose in strict source Python settings; opt-in configurations explicit. Do not
broaden unsupported cage/legacy interfaces. Tests for compatibility, unknown input,
hash sensitivity, and actual native use before experiments.

### Task 3: Fresh builds and complete experiments

Preserve old build. New native build, relevant native and Python tests, independent
review. G6O6 old-policy replay checks old full archive; new policy G6O6 and G7O7
at SLQP960 each twice. No warm starts. Keep all failures and exact source/binary/
input/output hashes. Report QP attempts, nonlinear acceptance, fallback counts,
KKT, binary result, memory and time. Success of local directions is not convergence.
If verified binary result exists, examine geometric sealing prerequisites before
physical simulation. Otherwise report actual new blocker and next bounded task.
