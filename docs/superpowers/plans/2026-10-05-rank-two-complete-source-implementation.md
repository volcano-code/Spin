# Faster 完整秩二源构造实施计划

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox syntax. 本计划待用户审阅及执行方式确认；编写本文件不构成产品代码实施。用户保全工作区的要求优先于技能模板中的自动提交步骤：不全量暂存、不提交、不替换共享运行库。

**Goal:** 在现有 Faster 调用链中交付完整秩二源面、正侧封口、真实八维源／材料导数及独立验收，继续服务于 Faster 优先、随后两篇论文完整计算复现的总目标。

**Architecture:** 一个受绑定的精确源对象，经代数分区、有限载体和共同投影闭合，产生现有 `QuadricTriangulation`。同一构造对象进入 source motion 与实际 Cherchi／CGAL 材料微分，不重新生成另一份源面，不建立平行求解器或封存系统。

**Tech Stack:** C++17、Eigen、固定 CGAL/CORE/GMP/MPFR、现有 C++20 Cherchi 探针、Python 3.13、现有 source_job、编译观察器与材料检查器。

**Spec:** `docs/superpowers/specs/2026-10-05-rank-two-complete-source-design.md`，SHA-256 `c9bf8012cfc5cd295f951ed0c845cafd01b36b10ccb5ab3190c41c3dabf47982`。用户在 2026-10-05 明确“批准，继续推进目标”。该批准版本保持不变；其中原有“待审阅”字样由本批准记录承接，不改写历史文件。

## Global Constraints

- `R=/data/xvxinyao/moment-guided-spin-v1/code/make-it-balance-paper-sim`；沿用当前独立脏工作树，保留全部失败和已有结果。
- 新 profile 为 `rank-two-arrangement-source-v1`；它是独立秩二实现，不与论文“先三角化再截断”或旧有限片路径混称。
- 固定法向球、原半径、目标、姿态、密度与科学容差；不以凸包、规则网格、体素、扩大半径或丢分量替代完整构造。
- 精确有理系数是唯一几何权威；数值视图不是第二个多项式。正尺度归一化不得翻转材料侧。
- 源图与材料 Jacobian 保持 `10×8`，顶点速度保持 `3×8`；一般 Eq.14 仍拒绝未经证明的方向。
- 解析质量矩 `≤1e-10`；光滑导数 `≤1e-6`、至少三个步长；最终材料质量矩与约束 `≤1e-6`。
- 现有一般 Faster `α=t_max=0.5`、raw residual `≤1e-12`、每次最多 60 次求值不变。开发种子 0、1、每种子最多八次尝试不变；本批不在低维路径上假造通用 Newton 运行。
- 编译最多四并发；大型任务顺序；单作业 24 GiB、swap 0、60 分钟；数值库单线程。起始空闲保持 100 GiB，运行中低于 20 GiB 安全退出。
- 每个作业根独占且不可重用；只改检查器时重审保存对象，只有有效求解输入改变才重跑优化。
- 本批通过不等于完整 Faster、通用秩变换或原论文案例通过；近锥面 CAD、原始资料、规模实验、壁厚、仿真与 Spin-It 仍保留在总体清单。

## Review Focus

1. 对四次曲线做平方自由化后丢失偶重因子，导致材料符号错误：任务二保留原多项式与重数，并测试平方曲线不翻号。
2. 代数类型只有比较、没有算术，或有限精度符号查询把未知返回为零：任务一、二验证表达式适配和精确／未知区别。
3. 圆环孔中整条射线材料被漏掉，或全在边界的顶点生成零厚度：任务三、四分别检查载体厚度和完整闭合实体。
4. 共享法向边、赤道或三角化分支改变后仍复用旧速度：任务五测试共享身份、赤道定义及失效处理。
5. 新构建的测试实际上加载旧 native／Cherchi 探针，或写到作业根外：任务一、六、八锁定路径和完整依赖；不得以旧路径通过代替当前结果。

---

## 一、文件边界与统一接口

仅新增以下职责明确的单元；以下名称是计划中的新文件／接口，不是当前已存在能力。

| 文件 | 单一职责 |
|---|---|
| `include/spin/RankTwoSource.h` | 精确输入、不可伪造来源附件和预算的公共合同 |
| `src/spin/RankTwoSource.cpp` | 精确系数准入、Schur／正交图、身份与单位视图 |
| `src/spin/detail/RankTwoAlgebraic.h`、`src/spin/RankTwoAlgebraic.cpp` | CGAL 代数数到可算术表达式的认证适配 |
| `src/spin/detail/RankTwoConstruction.h` | 内部分区、弧、点依赖、载体及只读构造计划类型 |
| `src/spin/RankTwoArrangement.cpp` | 原多项式符号、完整曲线分区、规范共享事件 |
| `src/spin/RankTwoCarrier.cpp` | 有界弦化、孔洞／多分量三角化、严格厚度见证 |
| `src/spin/RankTwoClosure.cpp` | 正负材料共同投影提升、球盖隔离、实际网格检查 |
| `src/spin/RankTwoSourceMotion.cpp` | 固定构造分支内的精确图／顶点运动 |
| `tests/rank_two_*_regression.cpp`、`tests/fixtures/rank_two_source_cases.h` | 下述任务的独立可执行正反例 |
| `experiments/spin_it_reproduction/faster/run_rank_two_source_batch.py` | 复用 source_job 的窄范围构建／重载／冻结状态批次，不授予论文完成 |

继续修改现有 `QuadricTriangulator`、`QuadricSourceMotion`、`QuadricDiscreteMomentDifferential`、`QuadricMomentEvaluator`、`bindings/pybind_module.cpp`、严格导数运行器和注册测试。旧分支保持原值与原哈希。

### 对外签名与约定

以下类用私有 `Impl` 保存精确／代数对象；工厂之外不能写入字段。`TriangleMesh`、`QuadricTriangulationSettings`、`QuadricTriangulation`、`QuadricStratumSourceMotion` 沿用现有定义。

```cpp
using RankTwoExactCoefficients = std::array<CGAL::Gmpq, 10>;

struct RankTwoSourceLimits {
    int maximumRefinementDepth = 24;
    int maximumIsolationBits = 512;
    std::size_t maximumArrangementCells = 200000;
    std::size_t maximumVertices = 1000000;
    double maximumDeviationOverInputScale = 0x1p-20;
    double maximumCoordinateErrorOverInputScale = 0x1p-44;
};

class PreparedRankTwoSource {
public:
    static PreparedRankTwoSource fromExact(
        const RankTwoExactCoefficients&, const std::string& frameIdentity);
    static PreparedRankTwoSource fromRepresented(
        const QuadricCoefficients10&, const std::string& frameIdentity);
    const RankTwoExactCoefficients& exactCoefficients() const;
    const QuadricCoefficients10& normalizedView() const;
    const std::string& inputHash() const;
    const std::string& geometryHash() const; // only positive scale canonicalized
private:
    struct Impl;
    explicit PreparedRankTwoSource(std::shared_ptr<const Impl>);
    std::shared_ptr<const Impl> impl_;
};

class RankTwoSourceConstruction; // read-only definition in internal header
class RankTwoArrangement;
class RankTwoCarrier;

RankTwoArrangement analyzeRankTwoArrangement(
    const PreparedRankTwoSource&, const QuadricTriangulationSettings&);
RankTwoCarrier buildRankTwoCarrier(
    const PreparedRankTwoSource&, const RankTwoArrangement&,
    const QuadricTriangulationSettings&);
QuadricTriangulation closeRankTwoCarrier(
    const PreparedRankTwoSource&, const RankTwoCarrier&,
    const TriangleMesh& design, const QuadricTriangulationSettings&);

QuadricTriangulation triangulateQuadric(
    const PreparedRankTwoSource&, const TriangleMesh& design,
    const QuadricTriangulationSettings&);

QuadricStratumSourceMotion evaluateQuadricStratumSourceMotion(
    const PreparedRankTwoSource&, const QuadricTriangulation&);
```

`QuadricTriangulationSettings` 新增 `rankTwoArrangementSource=false` 及 `RankTwoSourceLimits rankTwoLimits`。与旧有限片、projective 或 algebraic-finite profile 同时启用时拒绝；禁止隐式选其一。`QuadricTriangulation` 新增 `shared_ptr<const RankTwoSourceConstruction> rankTwoConstruction`；旧数组仍只表示旧来源。

新 overload 接受完整设计网格，以支持封口包含／距离判断。旧只给设计顶点的入口只可建立严格径向上界；若该证书不足，不能假造完整设计实体或转用顶点采样距离，返回缺失完整设计证据。新应用一律提供实际设计网格。

统一失败使用携带 `code`、`witness` 的异常／结果附件，由现有调用边界映射为 spec 第 10 节状态；不得用空向量表达未知。普通旧接口的异常约定不变。

内部 `RankTwoArrangement` 和 `RankTwoCarrier` 的定义均位于 `src/spin/detail/RankTwoConstruction.h`，对外只读。前者拥有 DCEL 的规范顶点／半边／面记录、各面全部边界环及孤立点；后者拥有点表达式 DAG、约束边、三角形、原区域 ID、上下高度和严格厚度见证。原生 CGAL handle 只用于本次构造，不能写入持久身份。

持久化仍作为现有封存对象的附件，schema 为 `rank-two-source-construction-v1`，不是另一套验收系统：

- 精确有理数统一为互素 `numerator/positive_denominator` 字符串；保存系数顺序、原多项式与正尺度规范形式。
- 代数根保存定义多项式、根序号和 dyadic 隔离区间；表达式 DAG 保存运算、依赖 ID 及分支。禁止只保存打印后的小数坐标。
- 顶点保存 binary64 位模式、误差界和对应表达式 ID；面保存顶点 ID、区域／法向 owner 和 cap 标签。
- 保存 input、normal、settings、chart、combinatorial-plan、source、backend 摘要。隔离区间的加密记录进入证明摘要，不改根 ID；时间／内存不进入科学摘要。
- 重载先校验摘要，再独立验证原多项式符号、全部环／邻接、实际有限面和误差界。摘要相等不能替代几何验证；字段中不存在可直接授予通过的 `passed` 输入。

## 二、运行与构建合同（所有任务共用）

所有下面的构建命令由新批次驱动器放进现有 `run_source_job` 的工作载荷。`J` 为其创建的新作业根，不能预先存在；配置、编译、链接、测试与输出全部在 `J` 内。这是将来执行命令，本计划编写期间不运行。

```bash
/usr/bin/cmake -S "$R" -B "$J/build" -G Ninja \
  -DCMAKE_MAKE_PROGRAM=/home/lin/anaconda3/bin/ninja \
  -DCMAKE_CXX_COMPILER=/usr/bin/c++ -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DLIBIGL_SKIP_DOWNLOAD=ON \
  -DMAKE_IT_STAND_LIBIGL_DIR="$R/libigl" \
  -DMAKE_IT_STAND_BUILD_VIEWER=OFF \
  -DMAKE_IT_STAND_BUILD_REACHABILITY=OFF \
  -DMAKE_IT_STAND_BUILD_SPIN_IT_2014=OFF \
  -DMAKE_IT_STAND_BUILD_SPIN_IT_FASTER=ON \
  -DMAKE_IT_STAND_BUILD_EXACT_GEOMETRY_CHECKS=ON \
  -DMAKE_IT_STAND_BUILD_PYTHON=OFF -DBUILD_TESTING=ON \
  -DPython_EXECUTABLE=/home/lin/anaconda3/bin/python3 \
  -DPython3_EXECUTABLE=/home/lin/anaconda3/bin/python3 \
  -DBoost_DIR=/home/lin/anaconda3/lib/cmake/Boost-1.82.0 \
  -DGMP_INCLUDE_DIR=/usr/include/x86_64-linux-gnu \
  -DGMP_LIBRARIES=/usr/lib/x86_64-linux-gnu/libgmp.so \
  -DMPFR_INCLUDE_DIR=/usr/include \
  -DMPFR_LIBRARIES=/usr/lib/x86_64-linux-gnu/libmpfr.so \
  -DMAKE_IT_STAND_OBSERVE_COMPILERS=ON \
  -DMAKE_IT_STAND_COMPILER_EVIDENCE_DIR="$J/compiler-observation" \
  -DMAKE_IT_STAND_COMPILER_TIMEOUT=600
```

已核对：这里不存在 `/usr/bin/ninja`。现有 pin 为 CGAL `f7c3c8212b56c0d6dae63787efc99093f4383415`、GMP `6.3.0`、MPFR `4.2.1`；CORE 还依赖当前可用的 Boost 库。实际配置继续执行 `cmake/SpinItFasterDependencies.cmake` 的强制核验，不自动安装或下载。

作业驱动器对模式 `backend / native / material / frozen / regression` 使用显式白名单命令；调用例：

```bash
/home/lin/anaconda3/bin/python3 \
  experiments/spin_it_reproduction/faster/run_rank_two_source_batch.py \
  --mode backend --output \
  /data/xvxinyao/moment-guided-spin-v1/code/make-it-balance-paper-sim/.superpowers/sdd/2026-09-29-seven-day-faster-sprint/rank-two-complete-source/backend-red-01
```

驱动器由任务一实现，后续任务扩展白名单。复用 `run_source_job`：backend 600 秒／4 GiB；其他模式最高 3600 秒／24 GiB；磁盘和 swap 门保持全局值。controller 必须是已核对的 Python 3.13，保留 timerfd、user-systemd/cgroup v2 真实检查，不降级为普通 subprocess 冒充受控作业。

输入清单绑定源码／未提交有效差异、工具、依赖、配置、checker；构建后记录真实对象、链接与加载路径。构建观察器本身不等于完整构建证明。后续导数／冻结作业可读取已封存的当前构建包，但须重哈希其全部有效依赖；不为每个小测试重新编译全部 core，不复用摘要不明的共享二进制。

每轮 RED/GREEN 使用独立目录；文件已存在就拒绝，不覆盖。作业返回 `PROCESS_EXITED && returncode==0` 只说明进程完成，科学门由输出重审决定。

## Task 1：精确源工厂与真实代数后端准入

**Files:** 新建 `include/spin/RankTwoSource.h`、`src/spin/RankTwoSource.cpp`、`tests/rank_two_source_input_regression.cpp`、`tests/rank_two_algebraic_backend_probe.cpp`、公共夹具头与上述批次驱动器；修改 `CMakeLists.txt` 的 Faster 条件块。

**Interfaces:** 实现统一签名中的两个工厂和只读 accessor。新 backend 探针直接链接 `CGAL::CGAL` 与 `CGAL::CGAL_Core`，不依赖先编译整个 core。

- [ ] 写失败测试：解析正则秩二、重复非零曲率、精确满秩但最小特征值很小、线性核项为零、错误 frame／附件、正负整体缩放。所有测试使用运行时 `require`，不能依赖 Release 下被移除的 `assert`。

```cpp
inline void require(bool ok, const std::string& why) {
    if (!ok) throw std::runtime_error(why);
}
// q=z-x*x-y*y; rational values, not decimal reconstruction of old results
RankTwoExactCoefficients a{0,0,0,1,-1,-1,0,0,0,0};
auto p=PreparedRankTwoSource::fromExact(a,"target-frame-v1");
auto twice=a; for(auto& x:twice) x*=2;
require(p.geometryHash()==PreparedRankTwoSource::fromExact(
    twice,"target-frame-v1").geometryHash(),"positive scale changed geometry");
auto opposite=a; for(auto& x:opposite) x=-x;
require(p.geometryHash()!=PreparedRankTwoSource::fromExact(
    opposite,"target-frame-v1").geometryHash(),"material sign lost");
```

- [ ] 注册 `rank-two-source-input-regression` / `spin-faster-rank-two-source-input`，先记录 API 缺失或准入错误的真实 RED，不把缺第三方依赖当作算法 RED。
- [ ] 工厂用精确行列式、非零主子式和核线性项验证；选绝对主子式最大的主元，平局按坐标序；保存比较分支。正尺度清除分母和公因子生成 geometry hash，保留原输入 hash 和方向。生成数值视图但不再调用它的梯度进行精确判断。
- [ ] 实现图 `u,E,β,h,G,F` 的精确多项式恒等式测试，`q(X)=0` 和 `uᵀG=β` 必须逐系数成立；加入非零交叉项及有理旋转输入。
- [ ] backend 探针明确实例化 `Arr_algebraic_segment_traits_2<CORE::BigInt>`，插入圆／四次曲线、直线和孤立点，核对真实 arrangement 计数与隐藏闭环。缺 CORE 必须编译／配置失败，不能复制 upstream 的提示后 `return 0`。

```cpp
#if !CGAL_USE_CORE
#error "Rank-two algebraic backend requires the pinned CORE configuration"
#endif
using Traits=CGAL::Arr_algebraic_segment_traits_2<CORE::BigInt>;
using Arrangement=CGAL::Arrangement_2<Traits>;
Traits traits;
auto construct=traits.construct_curve_2_object();
auto x=CGAL::shift(Traits::Polynomial_2(1),1,0);
auto y=CGAL::shift(Traits::Polynomial_2(1),1,1);
Arrangement circles(&traits);
CGAL::insert(circles,construct(x*x+y*y-1));
require(circles.number_of_faces()==2,"closed circle not analyzed");
Arrangement rings(&traits);
const auto s=x*x+y*y;
CGAL::insert(rings,construct(4*s*s-12*s+7));
require(rings.number_of_faces()==3,"two closed contours not analyzed");
```

- [ ] 批次驱动器 backend 模式真实执行新增 probe；随后构建工厂 target 并运行测试。每个 task 结束保存源码 diff、构建／测试回执与限制，不自动 commit。

## Task 2：完整代数分区与算术适配

**Files:** 新建 `src/spin/detail/RankTwoAlgebraic.h`、`RankTwoConstruction.h`、`src/spin/RankTwoAlgebraic.cpp`、`RankTwoArrangement.cpp`、`tests/rank_two_arrangement_regression.cpp`。

**Interfaces:** 内部 `RankTwoArrangement` 是只读 pimpl，提供区域／弧／事件记录及序列化；入口 `analyzeRankTwoArrangement(const PreparedRankTwoSource&, const QuadricTriangulationSettings&)`。记录包括原 `F,h,D`、零集因子及重数、代数分支、法向 owner、内外环与符号，后续不得只消费匿名坐标。

这里 `D` 用作集合时沿用 spec；用于下面的符号判定时明确指多项式 `D(v)=R²−vᵀEᵀEv`，所以内部为 `D>0`，不是相反号约定。

- [ ] 先固定公共夹具：disk 的系数为任务一 `a`；annulus 为 `{2,0,0,1,-1,-1,0,0,0,0}`；two-lobes 为 `{-2,0,0,1,-1,1,0,0,0,0}`；`R=3/2`，二维三角形为 spec 固定的三个点。测试完整面内区域分别为一圆盘、一圆环、两个分量，均不是由边根数量推断。
- [ ] 为代数表达式适配写 RED：`Algebraic_real_1` 只有比较，不能直接计算中点或行列式。用定义多项式与经认证的 dyadic 隔离区间构造 `CORE::rootOf(polynomial, BFInterval)`；比较根、符号、隔离区间和重载身份。包含 `sqrt(2)`、近根、相等代数表达式、无效隔离区间。
- [ ] 禁用有限 `max_prec` 的 `sign_at_2` 零值作为证明：采用无该上限的精确重载，或显式 UNKNOWN。在作业资源封装内运行；被中断不能返回“恰好为零”。
- [ ] 保留原多项式及重数，arrangement 仅使用 square-free 零集因子。把公共因子、重合直线去重为带 owner 集合的曲线，不能丢掉符号来源。

```text
original := {F, h, D}
zero_sets := square_free_factors(original) with owners and multiplicities
arrangement := exact_overlay(zero_sets, fixed_normal_face_halfplanes)
for each open cell:
    evaluate ORIGINAL polynomials at its certified sample
    D>0,F<0             -> CUT
    D>0,F>0,h<0         -> FULL
    D>0,F>0,h>0         -> EMPTY
    zero strata         -> preserve incidence; do not invent 2D area
```

这里 `exact_overlay` 是上面伪代码中的局部步骤，不是现有库 API。实现使用 bundled traits 的 `construct_curve_2_object`、`make_x_monotone_2_object`、`CGAL::insert` 和扩展 DCEL 访问；必须同时处理返回的弧和孤立点。DCEL 的 vertex、halfedge、face 记录分别保存稳定 feature ID、owner 集合、全部边界环与邻接。全局 `F,h,D` 分析缓存一次，各法向面只叠加所需三个半平面，避免无条件把 480 条无限延长线全交成巨大新离散。

- [ ] 固定反例 `F=(s−3/4)^2`：原 F 两侧非负，平方自由化后的圆不能把内部错误标成 CUT。加入共同因子、四交点、法向边反向和 equator 必要切触。
- [ ] 保存／重载事件，检查相邻 normal face 的规范 edge/root IDs 完全一致；与已保存 3,840 边认证交叉检查，不生成新的成功总数宣称。
- [ ] 注册并执行 `spin-faster-rank-two-arrangement`；依次保留 RED、GREEN 与未知／资源失败。未恢复全部区域就不进入网格任务验收。

## Task 3：有限共同载体与严格厚度

**Files:** 新建 `src/spin/RankTwoCarrier.cpp`、`tests/rank_two_carrier_regression.cpp`；补充内部 `RankTwoCarrier` 类型及只读来源记录。

**Interfaces:** `buildRankTwoCarrier(const PreparedRankTwoSource&, const RankTwoArrangement&, const QuadricTriangulationSettings&)` 返回公共二维顶点／边／三角形、区域标签、代数位置表达式、分支计划、曲线到弦误差及严格厚度见证。

- [ ] RED 测试首先复现“只有轮廓顶点的隐藏 disk”“只有赤道顶点的球投影盘”和圆环内部对角线的零厚度。不能只检查正总体积。
- [ ] 单调弧按端点／临界点横坐标的有理二分定义采样；竖直弧用纵坐标。每个新增点绑定曲线、分支序号、两端 feature ID 和分割比例；用精确表达式求坐标。不能使用根区间中点代替曲线点。
- [ ] 用隔离条带和区间导数界证明弦不跨越其他弧及误差界；临界端使用含端点的单调隔离扇区，无法认证则细分或 UNKNOWN，不能凭采样认定拓扑。
- [ ] 将所有共享弧离散先发布到公共登记表；每个 face 消费同一表。对全部外环／孔洞采用确定性约束三角化，约束与对角线按稳定 feature ID 排序。先在任务二验证 `Cartesian<CORE::Expr>` 的精确谓词，再采用 `Constrained_triangulation_2` 的 `No_intersection_tag` 路径；只输入已证明互不相交的约束，不允许库生成没有来源的交点。区域标签来自认证分区，不能用所有约束的奇偶数替代。无法分辨的零行列式不是可用三角形。
- [ ] 内部点从局部垂直分解的左右事件和上下边界构造：横坐标取区间内固定比例，纵坐标取该截面上下弧值的固定比例；保存完整移动表达式。对每条非外边界内部边和每个活跃区域检查厚度。

```text
For each active triangle, let d_i=upper_i-lower_i.
Require d_i >= 0 at every corner.
Require sum(d_i)>0 for its open interior.
For each edge which is not a true material-boundary edge:
    require d_first+d_second>0.
If a condition fails because only seam vertices were used:
    insert the certified moving interior/edge sample,
    retriangulate incident cells through the shared registry,
    recheck all conditions; never delete that region.
```

- [ ] 增加漏孔、重复 seam、反向环、未声明插值点、粗弦交叉、不同 owner 独立舍入等破坏测试。检查转换后的实际有限载体，不把代数 arrangement 的通过直接转授。
- [ ] 注册并执行 `spin-faster-rank-two-carrier`。输出实际有限顶点、面、各类区域与误差界；精度／细分超限保留 `RESOURCE_LIMIT`。

## Task 4：球内正侧闭合与网格出口

**Files:** 新建 `src/spin/RankTwoClosure.cpp`、`tests/rank_two_closure_regression.cpp`；修改 `QuadricTriangulator.h/.cpp`；复用 `ExactBoundaryAudit`、`ExactMaterialRegion` 和内部 `ExactBoundaryDistance`。

**Interfaces:** `closeRankTwoCarrier(const PreparedRankTwoSource&, const RankTwoCarrier&, const TriangleMesh&, const QuadricTriangulationSettings&)` 产生现有 `QuadricTriangulation` 及私有构造附件。`triangulateQuadric` 新 overload 串接任务一至四。

- [ ] RED：annulus 的 full-ray 区域被漏掉；两个 interface lobes 被误当最终材料分量数；所有 cap 顶点在球面而三角形穿过设计域。先固定失败几何，再实现。
- [ ] 按 spec 的三个射线区间提升同一个载体，统一去重 seam。upper 朝 `+u`，lower/interface 朝 `−u`，用载体绕向和统一图的方向符号装配，不逐面猜梯度方向。
- [ ] 从同一载体生成反侧用于互补检查。预先排除任意独立重新三角化；正体积交为空、并集为同一多面体球表示，质量矩相加一致。
- [ ] 对实际 binary64 顶点的精确有理解释检查封口三角形最小半径：

```cpp
using K=CGAL::Simple_cartesian<CGAL::Gmpq>;
// actualCap is a triangle of the binary64 output, interpreted as rationals.
const K::Triangle_3 actualCap(a,b,c);
const auto minimumRadiusSquared=CGAL::squared_distance(
    K::Point_3(0,0,0), actualCap);
const bool radialCertificate =
    minimumRadiusSquared > squaredDesignRadiusWithClearance;
```

上述片段中的 `a,b,c` 从当前面三个实际顶点以 `Gmpq(double)` 转换；`squaredDesignRadiusWithClearance` 是设计域最大半径加既定 clearance 后向外有界取值的保守平方界，不是只看 cap 顶点。`radialCertificate=false` 不能立即伪装成功或强制拒绝：继续用 `verification::detail::exactBoundaryDistance` 的开放非退化三角集接口与完整域分类。不能把开放 cap 伪装成 `ExactMaterialRegion`。距离为正后还须证明 cap 在域外，防止它完全位于域内；两条认证路径都不能证明时返回 `CAP_DOMAIN_INTERFERENCE` 与见证。

- [ ] 检查实际输出闭合、边／点流形、绕向、退化、自交、严格厚度、覆盖和 cap 隔离。用现有 `auditExactTriangleBoundary(..., OutwardPositive)`；失配给出 face／vertex IDs 和来源，不静默修补。
- [ ] 接通新设置、附件和请求哈希；四种旧 flag-off identity 仍逐字节保持。新旧 profile 不共享缓存；负例证明只换附件不能让旧网格获得新身份。
- [ ] 注册并执行 `spin-faster-rank-two-closure` 与旧 `spin-faster-fixed-normal-request-hash`。所有解析正向集合必须获得实体，而不只是“正确拒绝”。

## Task 5：固定分支源运动与真实参数基

**Files:** 新建 `src/spin/RankTwoSourceMotion.cpp`、`tests/rank_two_source_motion_regression.cpp`；修改 `QuadricSourceMotion.h/.cpp` 和附件序列化。

**Interfaces:** 现有返回类型保持不变；新 overload 必须接受任务四同一个 triangulation，校验 input/source/plan hash 后才计算。附件不匹配不能重新构建一份网格再忽略差异。

- [ ] RED：新增点速度全零、误用旧挤出索引、系数基十列填充、相等曲率拒绝、旧图／旧采样计划被接受。
- [ ] 在固定 Schur 主块及齐次 gauge 内计算八个参数的导数；将 `D(a/||a||)/Dθ` 作为系数基，同时对 `u,E,h` 作链式微分。记录 basis 与 parameter order，不能假定新基天然正交。

```text
A = [[C,w],[w^T,w^T C^-1 w]]
dA_kk = 2 dw^T C^-1 w - w^T C^-1 dC C^-1 w
d(a/||a||) = (I-aa^T/||a||^2) da / ||a||
```

- [ ] 法向边球交点用 spec 三方程隐式系统；临界点用其实际两个平面图方程；采样／内点逐依赖微分；赤道用球面、equator 及来源约束，不除以零高度。共享点返回同一个速度记录。
- [ ] 三步长在同一中心 Schur 图内取 `h={1e-5,3e-6,1e-6}`，方向为八个实际参数方向。每个 stencil 都核验图、root、连接和 cap 分支；不一致单列事件，不并入通过均值。

```text
fd_j = (vertices(theta+h*e_j)-vertices(theta-h*e_j))/(2*h)
for each vertex:
    if ||velocity_j|| < 1e-8: require absolute_error <= 1e-10
    else: require relative_error <= 1e-6
require all 8 directions, all 3 steps, and a reported refinement/noise trend
```

- [ ] 检查正系数缩放的零源运动、坐标置换、有理旋转、重复曲率和赤道穿越。真正切触、分区改变和三角化分支改变必须使旧导数失效。
- [ ] 注册并执行 `spin-faster-rank-two-source-motion`；一般 source motion / Eq.14 回归仍必须拒绝把本结果当全维 Jacobian。

## Task 6：同一次实际 Boolean 的材料导数与 Python 接口

**Files:** 修改 `QuadricDiscreteMomentDifferential.h/.cpp`、`QuadricMomentEvaluator.h/.cpp`、`bindings/pybind_module.cpp`、`experiments/spin_it_reproduction/strict/run_faster_cherchi_differential_basis.py`、`python/make_it_stand/verification/cherchi_motion_differential.py`（只为确证缺口）；新增 `tests/python/test_rank_two_complete_material.py`。

**Interfaces:** 新 Python 入口返回现有 source-motion 字段并增加 `exact_source`、`construction_plan`、各摘要及状态附件。输入精确系数使用十个标准有理数字符串，通过 native 工厂构造；不由 Python 重新 normalize 后覆盖。扩展现有运行器 schema `cherchi-rank-two-arrangement-basis-method-v1`，保留旧三类 schema。

绑定一个不可直接实例化／改字段的 native `PreparedRankTwoSource` 类，新增函数名固定为 `prepare_spin_it_faster_rank_two_source(exact_coefficients, frame_identity)` 和 `evaluate_spin_it_faster_rank_two_source(prepared, design_vertices, design_faces, settings)`。前者调用工厂，后者返回带同一私有构造对象的 native 求值及上述只读导出字段。保存／重载通过现有附件入口重新核验，不允许仅传回 Python 字典就跳过工厂。现有 Cherchi 运行器仍负责一次 Boolean 与独立检查；这两个函数不隐式执行优化或构建第二份源。

- [ ] RED：绑定层丢失精确系数／把 Schur 基当正交基、构造摘要不一致、错误密度、错误 probe 路径、10 列冒充8列。
- [ ] 串接新 source motion 与现有 `differentiateRestrictedPolyhedralMaterialMoments`。Cherchi 路径复用当前精确交点来源传递：实际一次 Boolean 的 boundary、birth faces、construction records 同时供积分和微分；禁止为导数再运行一次 Boolean。
- [ ] 现有 strict runner 在新 schema 下以 exact Schur 参数生成 49 份求值（中心1 + 8方向×3步长×2侧），记录所有状态及完整依赖。旧 `_retract_paraboloid` 的谱投影仅留旧 schema，不能用在新证据。
- [ ] 固定小型实际设计盒 `[-1/4,1/4]^3`，`q=z−x²−2y²`、`R=3/2`、法向 level 2、密度 `3/2` 作为独立方法夹具；输入与步骤在首次求值前封存。它不提供优化初值、不计原论文案例、不保证所有分支预先成功。
- [ ] 新隔离 Python 构建启用 `MAKE_IT_STAND_BUILD_PYTHON=ON`；构建 `_make_it_stand` 并明确加载 `$J/build/python`。测试原型：

```python
self.assertEqual(record["derivative_coordinates"], "EXACT_SCHUR_RANK_TWO_8")
self.assertEqual(len(rows), 49)
self.assertEqual(np.asarray(rows[0]["jacobian"]).shape, (10, 8))
self.assertEqual(len({row["label"] for row in rows}), 49)
self.assertTrue(all(row["independent_moment_pass"] for row in rows))
self.assertEqual(record["probe_sha256"], digest_file(explicit_probe))
self.assertFalse(record["optimizer_executed"])
```

`digest_file` 沿用现有运行器摘要函数或在该测试以 `hashlib.sha256(path.read_bytes()).hexdigest()` 明确定义；不读取 caller 的 passed 字段。

- [ ] 运行已有 Cherchi 探针独立构建，并用 `--probe "$J/cherchi-build/cherchi_probe"` 明确传入：

```bash
/usr/bin/cmake -S "$R/experiments/cherchi_reference_probe" \
  -B "$J/cherchi-build" -G Ninja \
  -DCMAKE_MAKE_PROGRAM=/home/lin/anaconda3/bin/ninja \
  -DCMAKE_C_COMPILER=/usr/bin/cc -DCMAKE_CXX_COMPILER=/usr/bin/c++ \
  -DCMAKE_BUILD_TYPE=Release
/usr/bin/cmake --build "$J/cherchi-build" --parallel 4 --target cherchi_probe
```

原 `test_faster_cherchi_basis.py` 与 `test_faster_fixed_normal_runner.py` 有旧 probe 硬编码，增加显式测试配置注入，默认仅为旧兼容测试保留；新必需门无显式新路径就失败。runtime path 使用已有隔离 `compiler-runtime`，不覆盖共享路径。

- [ ] 保存实际输出，再由现有精确材料检查器重载积分与域分区。CGAL 对照使用同一源和设计操作数，分别报告实体、质量矩和导数；不因质量矩相近忽略实体差异。密度只乘一次。
- [ ] 成功方法再由独立进程重复整个49求值链；所有失败也保留。注册 Python 测试与运行器检查，不开启新的未知种子长批次。

## Task 7：原八状态、来源重载和跨进程重复

**Files:** 扩展 `run_rank_two_source_batch.py`；新增 `tests/rank_two_frozen_source_probe.cpp`、`tests/python/test_rank_two_frozen_batch.py`；复用已有根证据和严格检查机制。

**Interfaces:** probe 读取冻结 manifest 中的实际系数／法向球／设计域，输出当前实际 mesh、构造附件和失败见证。batch 保存每状态记录，不接受手填成功。

- [ ] 初始 manifest 精确锁定 `{0,1,2,18,56,57,58,59}`，不是任意八个不同 ID。以 `.superpowers/sdd/2026-09-29-seven-day-faster-sprint/` 为前缀：系数读取 `rank-two-edge-certificate/green-03/bundle/native.tsv` 的 native 表示；法向球读取 `rank-two-normal-boundary/diagnostic-01/bundle/normals.tsv`；设计实体读取 `finite-patch-coverage/seed1-01/output/design.obj`。与 `finite-cap-diagnostic/diagnostic-01/bundle/state-N.txt` 中保存的设计顶点逐个核验，绑定相关旧 config／完成回执；半径固定4，不自行换盒体。系数只从 native.tsv 读取，不能再次调用旧 `fromUnnormalized` 改变其位模式。
- [ ] 注册反向测试：缺 ID、重复 ID、错系数、半径增加、法向 level 变化、漏 full-ray 单元、附件与 mesh 不符、截断输出及只有一次运行冒充重复。

```python
expected = {0,1,2,18,56,57,58,59}
require(len(records)==8 and {r["state_id"] for r in records}==expected,
        "exact frozen state set required")
require(all(r["truncation_radius"]==4 for r in records), "radius changed")
require(all(r["normal_input_sha256"]==manifest["normal_input_sha256"]
            for r in records), "normal discretization changed")
```

- [ ] 每个状态顺序执行 source 构造／封口／实际边界检查；整个八状态批次共享原作业60分钟预算，不为每个状态重新获得60分钟后仍说同一作业。资源不足写明未完成，不跳过困难状态。
- [ ] 全部成功才可以标记完整八状态源构造通过。部分成功只报告对应 ID；未知／正确拒绝不能升级整个门。
- [ ] 在未变源码与配置下独立第二进程重复，比较源网格、实际材料（若已进入该层）、组合事件与来源摘要。计时、内存不进入科学摘要，原始资源记录保留。
- [ ] 批次结束重新读取所有输出，核对来源、分区、mesh 与代数证据，记录新结果与旧有限片失败的差别。不得将新 profile 通过反写为旧 profile 成功。

## Task 8：受影响集成、范围审计与交付

**Files:** `CMakeLists.txt`、现有严格目录／审计注册与进度文件；必要时修正 `tests/cmake/assert_exact_geometry_dependency.cmake` 的隔离运行参数，不改科学门槛。

**Interfaces:** 保持原逐论文审计入口和 `PaperRunEvidenceV2`／`SealedSpinArtifact`，只追加有版本的构造、导数与来源附件。

- [ ] 新 target 与 CTest 对应关系如下，全部实际注册并发现，不能仅编译 executable：

| Target | CTest |
|---|---|
| `rank-two-algebraic-backend-probe` | `spin-faster-rank-two-backend` |
| `rank-two-source-input-regression` | `spin-faster-rank-two-source-input` |
| `rank-two-arrangement-regression` | `spin-faster-rank-two-arrangement` |
| `rank-two-carrier-regression` | `spin-faster-rank-two-carrier` |
| `rank-two-closure-regression` | `spin-faster-rank-two-closure` |
| `rank-two-source-motion-regression` | `spin-faster-rank-two-source-motion` |

Python 新必需测试精确指定 `test_rank_two_complete_material.py` 和 `test_rank_two_frozen_batch.py`；它们分别验证真实执行和保存证据，不以模拟返回值代替所有 native 调用。

- [ ] 构建并运行新测试及既有受影响基线：

```bash
/usr/bin/cmake --build "$J/build" --parallel 4 --target \
  rank-two-algebraic-backend-probe rank-two-source-input-regression \
  rank-two-arrangement-regression rank-two-carrier-regression \
  rank-two-closure-regression rank-two-source-motion-regression \
  truncation-root-certificate-regression quadric-fixed-normal-paraboloid-regression \
  quadric-fixed-normal-request-hash-regression quadric-stratum-tangent-regression \
  quadric-source-motion-regression quadric-discrete-moment-regression
/usr/bin/ctest --test-dir "$J/build" --show-only=json-v1 \
  -R '^spin-faster-(rank-two-.*|truncation-root-certificate|fixed-normal-paraboloid|fixed-normal-request-hash|stratum-tangent|source-motion|discrete-moments)$'
/usr/bin/ctest --test-dir "$J/build" --parallel 1 --output-on-failure \
  -R '^spin-faster-(rank-two-.*|truncation-root-certificate|fixed-normal-paraboloid|fixed-normal-request-hash|stratum-tangent|source-motion|discrete-moments)$'
```

- [ ] 若运行 `exact-geometry-dependency-contract`，先把子构建根和 Ninja 路径从 CTest 参数传入：`-DTEST_OUTPUT_ROOT=$J/dependency-contract -DCMAKE_MAKE_PROGRAM=/home/lin/anaconda3/bin/ninja`。脚本内所有子配置沿用显式参数；不得仍向共享仓库 `build/` 写文件，也不得在 source_job 的最小 PATH 下默默换生成器。最小批次不把这个未修正的旧测试冒充已运行。
- [ ] 受影响 Python 集成覆盖现有 source-motion/basis、Cherchi Newton/Table/Pareto dispatch、材料重载与 frozen batch。新 profile 不向一般 Newton 偷送 reduced Jacobian；旧 full-rank 方法结果仅在依赖有效时复用。
- [ ] 全局数学与接口复核分别对照 spec 第1–13节。报告明确列出：完整源构造、局部八维导数、实际材料调用、通用必要方向、原论文实验；最后两项不能因本计划前几项通过而自动打勾。
- [ ] 保存当前源码 diff、构建／依赖／加载证据、所有失败、两个进程科学摘要、输出与独立检查。保留当前 Spin-It source-cage 失败和近锥面 CAD 未通过状态，不宣称整个项目全绿。

## 三、计划自检与执行交接

覆盖对应：spec 1–3→任务一／八；4–5→任务一／二；6→任务三／四；7→任务二／三／四／六；8→任务五／六；9–10→所有接口及负例；11→任务一至七；12→共同作业合同；13→本次设计批准记录。

只读规划发现的实际陷阱已经纳入任务：代数值缺少算术、bounded sign 的未知零值、square-free 重数丢失、开放 cap 不能伪装闭合材料、Python 测试加载旧 probe、子 CMake 越过作业目录。上述发现不等于已修复或测试通过。

本计划推荐 **Native 执行**：由主执行者顺序完成各任务、每步先 RED 后实现；阶段末由独立审查者检查整体。原因是图坐标、共享来源、载体与导数接口紧密耦合，核心代码拆给多个同时实现者更容易出现不一致；只读数学／构建核查仍可并行。若用户选择逐任务 subagent 实施，则各任务和审查按上述依赖串行，不能并行改公共类型。

自检结论：已逐节对应批准设计，并将 D 符号、代数算术、序列化、实际输入路径、私有对象绑定和 cap 双路径检查补成具体合同。现有主 CMake 的相关 target／CTest 名称及 bundled arrangement 插入 API 已只读核对。未执行任何 probe、编译、算法测试或源网格生成，因此没有新增计算通过状态。

当前只等待本实施计划的审阅及执行方式确认，不再重复要求批准已通过的总体目标或书面设计。批准后从任务一的真实代数后端与精确输入失败测试开始。
