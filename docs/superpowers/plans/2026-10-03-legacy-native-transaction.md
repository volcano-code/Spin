# LEGACY Native Transaction Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans for inline execution or superpowers:subagent-driven-development if the user chooses delegation. Execute task-by-task; do not skip failed gates.

**Goal:** 接通已批准的可选 LEGACY 原生预测事务入口，完成真实校正、完整回滚、不可回滚成本与限定工程验证，不改变严格 DUAL 求解。

**Architecture:** 独立公开入口，共用现有 `solveSpinMomentSystemImpl` 内核。数值状态与成本状态分离，单步认证通过后再启用至多三步暂算；真实几何入口始终拥有认证和终点材料权限。

**Tech Stack:** C++17、Eigen、现有 CGAL/Cherchi 几何路径、pybind11、Python 3.13、CMake/Ninja、pytest、现有 source-job/运行证据系统。

**Spec:** `docs/superpowers/specs/2026-10-03-legacy-native-transaction-design.md`，2026-10-03 用户“批准”的版本 SHA-256：`2f4ba33b4de40c3acca717f344f6e43271070dc1017f141f41ea0570c30e8bcd`。先完整阅读该版本；以下路径相对仓库根目录。

## Global Constraints

- 工作树：`/data/xvxinyao/moment-guided-spin-v1/code/make-it-balance-paper-sim`。保留脏工作树、旧构建和所有失败；不自动提交、暂存用户修改或替换运行库。
- 严格 `DUAL_MONOMIAL_RAW_V1` 两层守卫不变；新入口只接受 `LEGACY_NATIVE_NORMALIZED_V1`。
- 不改原目标、容差、线搜索、重启次序及预算；操作方向矩阵不是数学 Jacobian。
- 固定材料真实积分后变换一次、加一次；响应仅包含已乘可变密度的目标坐标质量矩。
- 原生初始化产生检查点；缓存记录不是历史 Newton 状态。暂算不能授予几何、收敛或制造通过。
- 单步先验收，之后最多三步；两块连续失败禁用预测；失败成本、查询、逻辑预算和重启计数不回滚。
- 冻结八状态面板最多 16 次工程运行，来源映射准入后才预约；失败／崩溃不退额度。不训练、不生成教师标签、不做 40 组合加速验收、不外发。
- 编译最多 4 并发，大型作业顺序；24 GiB、swap 0、60 分钟，数值库单线程；开始至少 100 GiB、运行低于 20 GiB 安全停止。
- 代码改变后运行受影响测试；设计／报告编辑不触发科学重跑。本任务不关闭论文曲面、CAD、原资产、制造和其他独立缺口。

## Review Focus

1. predictor 重入求解器或修改请求：拒绝重入，原生请求拥有值副本；Task 1/4 测试。
2. Python `KeyboardInterrupt`／作业取消不能被异常回退吞掉：Task 4 测试终止传播。
3. 真实检查耗尽最后额度时不得免费 fallback 或以暂算结果结束：Task 3/5 测试。
4. 真实透传自身调用真实求值器，成本不得藏在预测 timing：Task 3/6 测试统一调用 ID。
5. 预算预约后进程崩溃、账本部分写入或并发预约：Task 6 测试额度不恢复、拒绝破损账本和重复身份。

## 文件责任与构建约定

| 文件 | 责任 |
| --- | --- |
| `include/spin/LegacyPredictionTransaction.h`（新） | 预测请求／响应、事务设置、身份及事件类型；不依赖 solver 头，避免循环 include |
| `src/spin/LegacyPredictionTransaction.cpp`（新） | 身份验证、事件及成本结构的纯规则；不另写 Newton |
| `include/spin/QuadricHollowingSolver.h` | 新增命名入口、版本化结果附件；旧签名保持 |
| `src/spin/QuadricHollowingSolver.cpp` | 原循环内部状态封装、事务与真实校正；保留默认调用顺序 |
| `bindings/pybind_module.cpp`、`python/make_it_stand/spin.py` | request-backed 新接口、异常和类型转换；复用真实固定材料积分 |
| `tests/legacy_prediction_transaction_regression.cpp`（新） | 原生内核、状态、预算测试；解析映射仅作内核检查 |
| `tests/python/test_legacy_prediction_transaction.py`（新） | Python 准入与真实几何集成 |
| `python/make_it_stand/verification/legacy_transaction.py`（新） | 版本化附件独立检查、预算预约；依附现有证据接口 |
| `tests/python/test_legacy_transaction_evidence.py`（新） | 账本、失配、确定性、崩溃与并发负例 |
| `experiments/spin_it_reproduction/strict/run_legacy_transaction_engineering.py`（新） | 预登记方法测试与八状态工程面板；不是论文运行器替代 |
| `CMakeLists.txt` | 注册上述新源／测试，与现有 Faster 开关一致 |

测试 executable 定名 `legacy-prediction-transaction-regression`，链接 `make_it_stand::core`。
测试注册名 `spin-faster-legacy-transaction-<case>`，case 为 `contract`、`disabled`、
`single-step`、`rollback`、`budget`、`three-step`。缺少依赖必须失败，不能 skip。

定义本次隔离构建目录 `build/legacy-transaction-20261003`；不得在旧 bundle 中编译。
配置复用当前精确构建的实际依赖路径和开关，不盲用缺精确几何的 headless preset。
先读取 `build/strict-algorithms-exact-current-20260922/CMakeCache.txt`，将有效显式
选项写入新证据目录的 configure 清单，保留新 configure 的输出和自己的 cache。
不得复制旧 cache 伪装新构建来源。依赖路径不存在时停止该构建并报告。

以下为作业内执行的验证命令；全部通过现有 `run_source_job` 封装，不能直接在
无资源约束的 shell 跑大型构建。最小 driver 如下，`identity_files` 必须在启动
前覆盖此次有效源码、配置、依赖、可执行文件和检查器，不能只填一个占位文件：

```python
from make_it_stand.verification.source_job import run_source_job
result = run_source_job(job_dir, argv, identity_files,
    maximum_wall_ns=3600_000_000_000,
    maximum_memory_bytes=24*1024**3,
    minimum_start_free_bytes=100*1024**3,
    minimum_free_bytes=20*1024**3)
assert result.returncode == 0, result
```

原生循环开发命令（新构建配置成功后）：

```sh
/usr/bin/cmake --build build/legacy-transaction-20261003 --parallel 4 --target legacy-prediction-transaction-regression quadric-hollowing-solver-regression spin-moment-subspace-regression
/usr/bin/ctest --test-dir build/legacy-transaction-20261003 -R 'spin-faster-(legacy-transaction|newton|subspace)' --output-on-failure
```

新增 target 必須真正存在，执行前检查 `ctest -N` 注册清单；旧记录不证明新代码通过。

## Task 1：定义隔离契约并锁定旧入口基线

**Files:** 新契约头／源、新 native 测试、solver 头、CMake。

**Interfaces:** 定义以下结构；`Vector10` 和 `Matrix10` 为 Eigen 固定维类型。
字符串身份使用现有规范摘要实现，不用 `std::hash`。`contextIdentity` 绑定
完整请求、输入网格、目标坐标、单位、密度、求值器配置和教师 profile；
`stateIdentity` 绑定 q、检查点、状态版本、事务和步号。

```cpp
using Vector10 = Eigen::Matrix<double, 10, 1>;
using Matrix10 = Eigen::Matrix<double, 10, 10>;
struct LegacyPredictionRequest {
    int schemaVersion = 1;
    std::string contextIdentity, stateIdentity, checkpointIdentity;
    std::string predictorIdentity, previousTentativeIdentity;
    std::uint64_t stateVersion = 0, transactionId = 0;
    int tentativeStep = 0;
    Vector10 q;
};
struct LegacyPredictionResponse {
    Vector10 variableMoments;
    Matrix10 operationalDirectionMatrix;
    bool valid = false;
    std::string contextIdentity, stateIdentity, predictorIdentity;
    std::string errorType, diagnostic;
    double reportedSeconds = 0;
};
struct LegacyTransactionSettings {
    int schemaVersion = 1;
    int maximumTentativeSteps = 1;
};
using LegacyPredictor = std::function<LegacyPredictionResponse(
    const LegacyPredictionRequest&)>;
```

- [ ] 添加契约失败测试，覆盖非 LEGACY、settings 0/4 步、身份缺失、q 非有限、
  响应 NaN、错误响应矩阵维数（Python 层）、同 q 不同状态版本、同材料输出
  hash 不同输入身份。原生测试以 `require` 抛异常模式沿用现有 native 测试。

```cpp
LegacyTransactionSettings bad;
bad.maximumTentativeSteps = 4;
bool rejected = false;
try { validateLegacyTransactionSettings(bad); }
catch (const std::invalid_argument&) { rejected = true; }
require(rejected, "four-step block admitted");
```

- [ ] 注册并运行 `--contract`，记录缺符号／缺验证的 RED；实现
  `void validateLegacyTransactionSettings(const LegacyTransactionSettings&)` 和
  `std::string validateLegacyPrediction(const LegacyPredictionRequest&, const LegacyPredictionResponse&)`。
  后者有效返回空串，无效返回确定性原因 token，异常诊断另存，不授予几何身份。
- [ ] 同一构建运行旧 `spin-faster-newton*`、subspace 回归，保存 hook 尚未接入
  时的科学 trace/event 摘要。检查真 baseline 为当前源码，不借历史二进制。
- [ ] 运行契约 GREEN 与依赖回归；保存差异和结果摘要，不提交整个脏树。

## Task 2：封装原生状态并验证 hook 关闭等价

**Files:** solver 头／源、新 native 测试。

**Interfaces:** 新 `LegacyTransactionResult` 包含原 `SpinNewtonResult solver`、
版本化事务事件和成本附件；同一 native 内核新增可选内部 transaction context。
新 native 内核测试入口：

```cpp
LegacyTransactionResult solveSpinMomentSystemLegacyTransactional(
    const SpinMomentProblem& problem, const SpinSolverSettings& settings,
    const SpinMomentEvaluator& trueEvaluator,
    const std::vector<QuadricCoefficients10>& initialGuesses,
    const std::string& contextIdentity, const Vector10& fixedMomentsTarget,
    const LegacyPredictor& predictor,
    const LegacyTransactionSettings& transactionSettings);
```

该低层入口只用于内核与原生适配；Python 不开放可手填 fixed moments 的入口。

- [ ] 添加 `--disabled` 测试，直接调用同一解析映射，两条路径比较 status、q、
  multipliers、traceHash、eventHash、真实调用次数和顺序。使用现有
  `quadric_hollowing_solver_regression.cpp` 的 identity map 与非终态 initialPoint
  数值；这是合成内核用例，不是合法材料见证。

```cpp
const auto problem = spin::SpinMomentProblem::staticBalance(
    spin::StaticObjective::MaximumMass);
spin::SpinSolverSettings settings;
settings.equationProfile = spin::SpinEquationProfile::LegacyNativeNormalizedV1;
Vector10 p;
p << 1.0, 0.12, -0.08, 0.09, -0.07, 0.05, -0.04, 0.03, -0.02, 0.01;
const std::vector<spin::QuadricCoefficients10> initial{
    spin::QuadricCoefficients10::fromUnnormalized(p)};
const auto map = [](const spin::QuadricCoefficients10& q) {
    spin::SpinMomentMapEvaluation value;
    value.rawMoments = q.values();
    value.momentJacobian.setIdentity();
    value.derivativeDefined = true;
    value.geometryHash = q.hash();
    value.topologyHash = "synthetic-identity";
    return value;
};
const auto base = spin::solveSpinMomentSystem(problem, settings, map, initial);
const auto wrapped = spin::solveSpinMomentSystemLegacyTransactional(
    problem, settings, map, initial, "kernel-fixture-v1",
    Vector10::Zero(), {}, LegacyTransactionSettings{});
require(base.traceHash == wrapped.solver.traceHash, "disabled trace changed");
require(base.eventHash == wrapped.solver.eventHash, "disabled events changed");
require(base.functionEvaluations == wrapped.solver.functionEvaluations,
        "disabled evaluator calls changed");
```

- [ ] 提取内部 `LegacyNumericalSnapshot`，包含当前 q、NewtonState、方向／缩放、
  line-search/iteration/restart 控制游标、缓存世代和认证轨迹前缀。定义独立
  `LegacyCostState`；真实查询、已消费逻辑试探、总重启、失败和墙钟只放后者。
- [ ] 空 predictor 分支直接进入既有行为，不产生兼容附件以外的新原科学事件。
  不改矩阵装配、求解／回退、初值筛选或 stoppingResidual。只共享内核，不能
  再写一份 while Newton 循环。
- [ ] 增加 snapshot 复制后修改 tentative 缓存的测试，认证缓存及原 trace
  不受影响；检查重入回调被拒绝而不是死锁。运行 `--disabled` 和旧 DUAL 回归。

## Task 3：实现单步暂算、真实认证和不可回滚账本

**Files:** solver 源、契约类型／源、native 测试。

**Interfaces:** 附件事件包含序号、parent/call ID、transaction ID、checkpoint、
state identity、角色、结果、原因及异常类型。成本角色：initial、trial、predict、
certify、fallback、final；接受／拒绝是结果字段，几何子操作以同 call ID 关联。

- [ ] 添加 `--single-step` 透传测试与 `--rollback` 故障测试。测试适配器由
  同一个计费真实 evaluator 获取变量矩与操作矩阵，所有 passthrough 真调用
  也被记录，不能只计入 predictor 自报时间。指定故障 nonfinite、exception、
  stale state、true-check-failure 分别覆盖；同时测固定矩非零时只加一次。
- [ ] 在现有 evaluateAttempt 的可选事务分支内构造预测请求，检验响应后构造
  **明确标为 tentative** 的 LEGACY 操作 map。不得向其填入伪造的真实
  geometry/topology/derivativeDefined 认证。原线性代数内部如需操作矩阵，
  使用有区别的来源标识，不能把它输出成已验证数学导数。
- [ ] 按以下协议在原接受位置调用真实校正，使用原有 evaluateState、
  stoppingResidual 和几何失败规则，原子提交只发生在全部检查之后：

```text
snapshot = certified numerical state
proposal = original native trial logic using tentative response
real = counted true evaluator(proposal.q)
if real valid and original acceptance predicate(real, certified):
    commit original numerical state from real; append certified trace
else:
    restore(snapshot) without touching consumed cost/budget
    append rollback; execute original true step if budget remains
```

- [ ] 为附件加入 append-only 宿主 timing 和 parent 子区间；真实调用数不能
  用回调自报值替代。科学摘要排除时间但包括调用角色、失败和接受顺序。
- [ ] 添加 `--budget`：最后一次额度被校正消费、失败回退无额度、初始真查询
  失败、predictor 重复拒绝、非有限 timing；验证计数单调，禁止隐式免费求值。
  明确分别保留 logical trial 与 true call 计数，两者均不超过原配置额度。
- [ ] 运行四个 case 与旧内核回归。透传在充足原预算的夹具上比较接受状态及
  终态；其成本事件／调用数不同于 hook-off，不要求额外开销凭空消失。若原
  预算不够则正确输出 RESOURCE_LIMIT，不能为取得相同终态增加预算。

## Task 4：接通实际请求、固定材料和 Python 入口

**Files:** solver 头／源、pybind、spin.py、新 Python 测试。

**Interfaces:** 新的公开 Python 签名固定如下；evaluator_settings 沿用已有
真实几何设置类型，不能从 predictor 获取真实认证。

```text
solve_spin_it_faster_legacy_transactional(
    request, *, fixed_shell_vertices, fixed_shell_faces, fixed_shell_density,
    evaluator_settings, predictor=None, predictor_identity="",
    maximum_tentative_steps=1, initial_quadrics=(), event_sink=None,
) -> dict
```

上面定义公开调用签名。实现必须依次：验证 LEGACY → 复制
请求及只读数组 → 积分固定网格并绑定其实际坐标/密度 → 原生真实初始化 →
同内核事务 → 原终点/互补材料真实求值 → 返回原结果与版本化附件。
零壳使用显式空顶点／面集合，记录零壳身份；不接受含 NaN 的“空壳”。

- [ ] 添加双层 profile 拒绝与几何驱动测试。Python 真实几何测试保存 request
  后改变原可变数组，原生身份和结果不受影响；密度或目标框架改变必须改变
  context identity。复用 `with_measured_fixed_shell_geometry` 积分语义，原生
  用 `axes.transpose()` 和 `-axes.transpose()*origin` 变换一次。

```python
def test_legacy_hook_keeps_strict_external_guard(self):
    with self.assertRaisesRegex(ValueError, "DUAL_MONOMIAL_RAW_V1"):
        solve_spin_it_faster_with_external_moment_evaluator(
            self.legacy_request, lambda q: {})

def test_cancel_is_not_converted_to_prediction_fallback(self):
    def interrupt(_):
        raise KeyboardInterrupt()
    with self.assertRaises(KeyboardInterrupt):
        self.run_geometry(predictor=interrupt)
```

`setUp` 定义 `legacy_request` 为从实际封闭立方体几何构造的固定配置；
`run_geometry` 调用本 task 新入口并提供相同固定材料几何。普通 ValueError/
RuntimeError 记录原异常类型并真实回退；KeyboardInterrupt/SystemExit 传播。

- [ ] C++ 顶层新增 `solveQuadricHollowingLegacyTransactional`：复用
  `solveQuadricHollowing` 的真实 evaluator、确定性初始化及终点导出流程；
  返回原 `QuadricHollowingResult` 加事务附件，不把核内合成 map 当材料验证。
- [ ] Python/native 对 predictor 重入作明确拒绝；新绑定必须持有／获取回调
  所需 GIL，不从工作线程无 GIL 调用 Python，不把异常字符串当成功响应。
- [ ] 原生重建后验证实际加载路径／hash，与构建回执绑定。运行新 Python
  测试和旧 strict external/源表面双材料受影响回归，不使用旧扩展假装通过。

## Task 5：在已通过单步协议上扩展最多三步

**Files:** solver 源、native/Python 测试。

- [ ] 建立冻结合成非终态夹具的三步故障序列：第 1/2/3 次认证各失败一次，
  检查整块认证 trace 无部分提交，缓存和乘子回到 snapshot。两块失败必须
  禁用 predictor；一次成功块清零连续失败数。
- [ ] 将一个候选改为有界候选链，长度取 settings.maximumTentativeSteps。
  pending convergence、异常、预算边界或不确定分支提前结束块。
  对链上每个状态真实求值，并与前一个**真实**状态执行原接受判定；不能
  比较预测残差，也不能只真检最终状态。

```text
real_chain = []
previous_real = snapshot.state
for candidate in tentative_chain:
    real = counted_true_evaluate(candidate.q)
    if not original_acceptance(previous_real, real): rollback entire block
    real_chain.append(real)
    previous_real = real
commit all real_chain only if every check succeeds
```

- [ ] 拓扑切换测两类：新实际 map 合法且真实残差下降可重新线性化；新 map
  无效或中间状态不能认证则整块失败。反向测试断言旧 Jacobian 未被继续使用。
- [ ] 运行 `--three-step`、单步及全部旧 Newton/subspace 用例。若单步未通过，
  不用增加 block 长度绕开，不启动工程面板。

## Task 6：证据检查、工程额度和实际几何运行器

**Files:** 新 verification 附件模块、Python 测试、新 engineering runner。

**Interfaces:** `verify_legacy_transaction_attachment(record: dict) -> dict`
读取事件和真实材料附件重算，不信 passed；
`reserve_engineering_run(ledger: Path, case_id: str, mode: str, input_hash: str) -> str`
返回一次性 reservation ID，mode 只允许 passthrough/fault；
`finish_engineering_run(ledger: Path, reservation_id: str, result_path: Path) -> None`。
预约使用进程间排他锁，追加带前一记录 hash 的日志并 fsync。尾部不完整拒绝
继续使用，修复需要人工检查，不能丢掉尾部释放预算。

- [ ] 先写负例：重复预约、17 次、并发同键、预约后异常、空统计、丢失回退、
  unknown checkpoint、错误 module hash、只有一次独立进程。使用临时测试
  账本，不写真实共享面板额度。

```python
first = reserve_engineering_run(path, "case-a", "fault", digest)
with self.assertRaises(ValueError):
    reserve_engineering_run(path, "case-a", "fault", digest)
# No finish call: crash must not refund the reservation.
with self.assertRaises(ValueError):
    reserve_engineering_run(path, "case-a", "fault", digest)
```

- [ ] 实现事件状态机检查，验证认证前缀、逐状态 true call、回滚前后状态、
  成本单调、两失败禁用及 final 来自真实材料。互斥成本与父子墙钟分别汇总。
- [ ] engineering runner 支持 `--admit-only`、`--method`、`--panel` 三种模式。
  admit-only 只核对冻结请求映射、源码构建和 checkpoint 可生成条件，不预约。
  每条 panel 只有实际设计域／固定材料／求值器设置可追溯才准入，缓存 q 和
  moments 不足则输出具体缺件。不生成历史 Newton 状态。
- [ ] method 使用独立的明确 SYNTHETIC_REGRESSION 真实几何配置；先封存
  请求、初值协议和身份，再运行 hook-off、passthrough、指定故障。不得读取
  已成功终态生成初值。按真实结果区分成功、数值失败和资源终止。
- [ ] 三条成本链分别输出：hook-off、passthrough、预测／故障。逐状态认证
  次数必须与实际真调用核对，不沿用历史三步一次查询比例；数据不足的加速
  上界返回 NOT_MEASURED，不把代理快速返回或缓存回放算成端到端加速。

## Task 7：集成验收与有限交付

**Files:** 新运行器／测试注册、已有 progress ledger 和对应结果目录。

- [ ] 新建独立结果目录，不复用 run ID。通过现有编译观察器保存新
  compile→link→load 回执；若新增源不在旧 commands 清单中，必须刷新新
  构建清单并核对实际对象，不能重放旧清单遗漏新源。旧 bundle 保留。
- [ ] 在两个独立进程运行同一冻结 method 配置，比较科学事件、真实认证状态
  与材料摘要；计时只列成本，不入确定性摘要。失败必须保留，不能挑成功种子。
- [ ] 仅在所有准入与单步／三步门通过后执行最多 16 次工程 panel。每次先
  持久预约再启动 source job，崩溃和未完成保留；不自动增加预算或重复运行。
- [ ] 运行完整受影响 Newton、strict CGAL/Cherchi 集成、材料与证据测试；
  公共 native 改动所影响的 Spin-It 回归也必须报告。历史 source-cage 失败
  不从统计中删除；专门通过不能写作全项目通过。
- [ ] 保存输入、源码／构建、每次失败、独立检查、材料和 SHA-256 索引。
  更新现有严格目录的扩展范围，不新增原论文通过项。不得用新字段升级旧结果。
- [ ] 提交 scoped 代码审查后再给最终实现结论；不自动 commit、merge、切换
  用户运行时、发送跨线程消息或宣称学习／速度通过。

## 自查与执行交接

覆盖：设计 §1–4 → Task 1/4；§5 → Task 2/3/5；§6 → Task 3/6；
§7 → Task 1–7；§8 → 文件责任、构建约定；§9 → Task 7 的声明限制。
Review Focus 的五项分别在所属 task 中加入测试。所有新实现尚未执行；
本计划的命令是待执行步骤，不是现有通过证据。

建议在当前会话由同一实现者顺序完成以上任务，再做独立审查；native 状态、
预算与几何认证强依赖，不适合并发修改共同内核。书面计划审核与执行方式
确认后开始 Task 1；不再重复设计阶段，也不把这份计划视作论文算法成果。
