# Spin-It 2014：保持可行性的盒约束 QP 与完整入口重跑

## 本次结果

修复了固定重建茶壶在旧 `FINITE_BOX_QP_V1` 下的工作集数值卡点。
新策略 `FINITE_BOX_PRIMAL_V2` 从同一输入运行两次，三个连续优化阶段均达到
各自转换前的一阶 KKT 门，但最终结果均为 **BINARY_VALIDATION_FAILED**。
没有通过二值材料验收，没有新封存几何或仿真视频，不是原论文完整复现。

项目：`/data/xvxinyao/moment-guided-spin-v1/code/make-it-balance-paper-sim`。
全部本轮证据：`build/spin2014-primal-box-20260915`。
前轮证据保持原样：`build/spin2014-finite-box-native-20260915`。

## 实现与保全

- 新增保持可行性的 primal active-set QP 路径。沿当前可行点到工作集最优点
  的线段移动至首个阻挡面，再检查或释放错误符号的乘子。
- 使用实际 L-BFGS 作用，不在生产求解器展开 N×N Hessian。
- Schur 对角均衡仅改变线性代数，不改变原始等式、驻点和盒约束残差门。
- 非零等式右端的 Phase I 使用独立最小范数可行性问题；Phase I 或主过程
  失败仍为 UNKNOWN，不能据此给出不可行证书，也不保证通用收敛。
- 默认、tangent、有限盒 V1 路径保持独立。新策略配置、trace 和 hash 明确区分。
- Python 独立 replay 对 V2 采用一致的均衡秩检查，但仍按原始量检查完整 KKT、
  边界、符号、互补性及接受后的下一状态。
- 审查发现并修复 pinned Eigen 零行乘法的未定义行为，以及浮点端点微小越界
  导致加入冗余工作面的情况；没有截断自由变量来伪造通过。

`pre-change.tar.gz`、`pre-change-tests.tar.gz` 保留修改前文件；
`implemented-source.tar.gz` 保存本轮实现与测试。未改主工作区、未提交或打 tag，
没有宣称当前脏工作树是可发布版本。旧结果相对新源码属于历史证据，
重新核验会触发 STALE_RESULT；没有改写旧 provenance。

## 新鲜验证

- 固定旧 iteration=108 的 45 变量、5 等式、实际 BFGS 历史作为回归输入。
  旧 V1 再现 RANK_DEFICIENT，新 V2 通过。
- 独立显式 dense BFGS 检查该小型 QP：等式误差 `2.01038e-16`，
  驻点误差 `1.04083e-17`，边界/符号/互补性误差为零。
  与 SciPy SLSQP 的目标差 `5.84385e-17`。dense 路径仅用于独立回归检查。
- `ctest.log`：11 个定向测试组全部通过，无 Not Run。
- `pytest.log`：96 项选定 Python 回归通过。
- `rounding-tests.log`：6 项；`stage-tests.log`：13 项检查器回归通过。
- `ubsan.log`：使用仓库 pinned Eigen 编译的解析 QP 测试通过 UBSan。
- `frozen-v1-compat.log` 与 `new-compat.log` 完全一致：三个小型 default、
  tangent、finite-v1 配置的完整诊断摘要在冻结旧模块与新模块间不变。
- `REVIEW.md` 记录两轮独立只读审查和修复。本轮采用先复现失败测试、
  再修复并独立检查的流程，而非仅凭原生 checked 字段验收。

上述测试有重叠，不将测试组与用例相加，也不冒称全仓库发布门已经通过。

新加载模块 SHA-256：
`56355966dc78e0ef3f5ed0ba0257dee213cc0bd80f6f5659afe22a656eb2b1bd`。

## 两次完整入口运行

`PROTOCOL.md` 预注册协议，`low.json` 相对 V1 唯一改变是 direction_policy。
输入仍为重建 G6 茶壶，γL=0.0001、初值、几何、容差及预算全部保持不变。
每阶段最多 960 次迭代、12 阶段、12000 叶节点；每次 24 GiB 地址空间与
300 秒墙钟上限，OMP/BLAS 单线程，顺序运行。两次正常结束，无资源超限。

| 阶段（从 0 开始） | 填充变量数 | 接受更新 | 转换前驻点残差 | 转换前等式残差 |
|---|---:|---:|---:|---:|
| 0 | 44 | 109 | 4.33681e-19 | 5.91147e-10 |
| 1 | 66 | 4 | 1.08420e-19 | 2.14671e-9 |
| 2 | 94 | 7 | 5.45215e-8 | 3.55555e-15 |

两次结果相同：三个阶段均 CONVERGED，全部 120 次 QP 被检查并接受，
接受步长均为 1，fallback/restoration 均为 0。
独立 replay 逐项检查全部 120 次更新：最大 QP KKT 误差 `4.00722e-15`，
下一状态误差为零。此检查不证明非线性全局最优。

运行器耗时分别 28.7314 / 28.6870 秒，峰值 RSS 337156 / 336744 KiB。
两次独立 run ID、127 个完整科学事件，归档 SHA-256 一致：
`7d0be4df50615f7ff7ebcee1ff2ca4a417414e868201bbb015fee2ad061ae956`。
输出 hash 一致：
`7ae660e22f0ebff8668377fcd9370eaf3d0acfa45bec6d0a52b76aeffa129704`。

`verification.json` 检查源码观察集合、未提交差异、构建配置、实际模块、
输入、完整 trace 和输出绑定。沿用明确的 CURRENT_SOURCE_OBSERVATION，
不把这一观察记录冒称为历史编译的外部认证。

## 当前剩余卡点：二值化

按固定 `PAPER_POST_TRANSITION_V1` 停止策略，第三阶段将 level-5 分数单元
细分到 level 6，终态分区为 252 叶节点，然后舍入。**没有在这 252 变量
分区上再次优化**，因此 `terminal_relaxation_kkt_applies=false`。
三个阶段的 KKT 不能直接挪用成最终分区的 KKT，更不能代表二值材料 KKT。

独立 `stage-audit-1.json` / `stage-audit-2.json` 确认：

- split/merge 的质量矩变化仅约浮点舍入量级，没有本次显著误差来源证据。
- 最终舍入前的归一化等式残差约 `2.37e-17`。
- 按固定 beta>=0.5 规则，186 个单元改变材料贡献。
- 舍入后最大归一化约束残差 `4.64702e-6`，门槛为 `1e-6`，失败项为横向
  一阶矩 X。独立重新计算质量矩误差仅 `2.44e-16`，不是质量积分明显不一致。
- 184 个 level-6 分数单元贡献 X 约 `3.64567e-6`；两个 level-4 单元贡献
  约 `1.00135e-6`。后两者 beta 为 0.0121571、0.0150713，小于既定
  fractional_epsilon=0.05，所以未继续细分，但最终仍会舍入为 0。

这定位了失败，不证明任何输入不可行。也不能断言仅提升最高分辨率一定通过：
粗层小分数单元的舍入贡献同样需要研究，改变细分/舍入规则必须单列扩展 profile。

## 下一工作及边界

1. 保持本轮协议和失败证据冻结，用固定输入研究真实更高采样/优化层级的
   舍入误差，同时记录各层贡献，不使用上采样结果冒充高分辨率优化。
2. 分开比较现有 unchanged-partition-KKT 数值扩展与当前 post-transition
   停止策略，不能悄悄改变论文重建 profile 的停止位置。
3. 若需舍入感知细分或二值约束修正，先核对方法来源，作为显式扩展单列；
   不改变阈值、删除失败或把改良算法冒称作者原算法。
4. 只有二值约束通过后才做材料封存、几何合法性、完整惯量及仿真验证。

本轮提升的是重建模型上的求解与诊断能力。原始模型/姿态/参数、level-9
完整案例、联合 cage、双材料切割、论文实验矩阵及冻结门仍未因此完成。
Spin-It Faster 未在本轮修改或升级完成声明。
