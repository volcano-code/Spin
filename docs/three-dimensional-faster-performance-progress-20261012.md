# 三维 Faster 性能收尾执行记录 — 2026-10-12

## 授权与范围

用户已批准落实“三维陀螺 Faster 正确性与速度收尾”，同时继续带柄心形和密集面条交付。沿用现有目标、材料、独立检查与重复条件；质量最小化、学习代理、浮力不进入当前实施批次。

Ruling: 在现有工作区直接实现，保留全部脏修改、历史输出与失败；不自动提交、合并或删除。执行技能的提交、清理和全项目测试步骤服从用户明确的保留与受影响测试约定。

## 已核验基线

- Flower v2 的四个二维求解调用耗时 0.458606、0.454262、0.466289、0.466672 秒。两个完整求解作业各包含 seeds 0/1，耗时 1.537255、1.556172 秒。不是三维球面 Newton 性能证据。
- Flower 材料准备 10.617567 秒；保存后独立材料检查 224.150972 秒；64 场景 PyBullet 作业 34.301947 秒。开发检查和动力学不计入普通优化时间。
- Noodles stage04 未认证。已有诊断定位：JSON 解析 0.095 秒；有理坐标规范化各约 1.5 秒；原生 `make_exact_material_region_from_rationals` 1708.743 秒。该诊断不授予前驱、整体输入或优化准入。
- 原生 assemble 内含精确边界审计、三角形相交检查、Polyhedron/Nef 构造及合法性检查。尚无各内部阶段的实测占比，不先猜测某一项可以省略。

## 当前实施批次

1. 已新增 `experiments/quadric_physics_demo_20261011/native_geometry_profile.py`：由 GDB 启动独立诊断子进程，每 60 秒采集一次最多 24 个原生栈，保留完整 MI 日志和真实退出码。采样类别不是精确阶段耗时；所有 admission 字段为 false。未修改现有 native 二进制或被冻结的 Heart 求值器。
2. 已新增 `launch_native_geometry_profile.py`：验证调度器命令及进程起始身份后，仅暂停 Heart 调度器，等待当前 point 完成，再运行单个受限诊断；正常或 Python 异常退出均由上下文恢复调度器。不停止 Heart worker，不并行运行两项大型作业。
3. 新增原生诊断测试先观察 3 项缺实现失败，再通过 3/3；串行调度测试先观察 2 项缺实现失败，再通过 2/2。真实 GDB 子进程退出成功/失败、采样及保留错误已测试；额外阶段链测试 5/5 通过。
4. 已启动 launcher，exec session 51779。运行目录 `experiments/quadric_physics_demo_20261011/three-dimensional-faster-performance-v1`；其 events.jsonl 当前记录等待 Heart point-job-23。诊断 job 为 noodles-native-job-01，24 GiB、swap0、60分钟。执行结束后恢复 coordinator PID1766941。
5. Heart point23 已退出0，串行原生诊断已实际启动。首个60秒采样只得到GMP内部调用与未知帧，分类 `EXACT_ARITHMETIC_UNRESOLVED_CALLER`。`nm -C` 证实现有模块无普通符号；不能由此断言是 PairScan 或 Nef 阶段。若余下采样仍不能定位，下一步必须采用独立、具有原生阶段计时或符号的候选构建，不能继续盲目堆叠相同采样。当前诊断的原材料检查继续按原函数执行，科学准入仍为 false。
6. Heartbeat 10:33UTC（18:33北京时间）：诊断仍运行，调度器按预登记保持暂停；没有重复大型作业。第121、181、241秒三个原生栈均出现动态符号 `verification::auditExactTriangleBoundaryRational`，其中一份叶调用为 GMP gcd/mul → CGAL orientationC2，调用链包含多层未知递归帧。这确认这些采样时刻正在精确边界审计，不能转述为该阶段的完整耗时或已精确定位 PairScan。普通符号缺失不代表动态符号全部缺失。等待完整采样与终态后决定候选；Heart point23作业回执退出0、551.504秒，现有队列尚未写其 POINT_COMPLETED 是调度器暂停造成的正常状态。
7. Heartbeat 11:03UTC：原生诊断退出0，1746.230秒；24个采样中23个落在精确边界审计，1个仅可见GMP。原dispatcher已按finally恢复并推进point24。证据表明边界审计是本次主要调查方向，仍不是阶段精确占比或全材料准入。
8. 已登记最小候选：共享索引边的非退化三角形，无需构造完整相交多边形即可用精确共面/半平面判定是否在公共边以外相交。三角形支持平面不同时仅可交于该边；共面时第三顶点位于同侧才存在面积重叠。所有判定仍为Gmpq精确算术，无容差、舍入或删检查。新头文件 `ExactSharedEdgeIntersection.h`；`ExactBoundaryAudit.cpp`增加该分支；`SPIN_DISABLE_EXACT_SHARED_EDGE_SHORTCUT`仅供独立基线对照。未覆盖/重建任何生产native模块。
9. 先新增完整CGAL相交构造oracle测试，并观察缺实现退出1（RED）。候选需通过近共面、同侧/异侧、倒序、一般有理坐标与退化拒绝，以及现有44接触/三组穷举/BVH回归。尚未编译通过，不授予候选通过。
10. 串行候选launcher `launch_shared_edge_candidate.py`已启动，exec session66196，目录 `three-dimensional-faster-performance-v1/edge-predicate-candidate-v1`。当前仅暂停dispatcher、等待point25结束；随后独立构建并运行上述测试，再顺序运行同一冻结面条边界的candidate与baseline审计。每个作业均24GiB/swap0/60分钟；完成或Python异常后恢复dispatcher。优先检查此候选events.jsonl和launcher进程，不能因上一个GDB launcher已完成就独立恢复dispatcher。边界审计加速不能代替Nef材料检查或整体输入准入。
11. 19:23北京时间复核：v1 launcher在构建准入前被 `UNSUPPORTED_JOB_ENVIRONMENT` 拒绝，没有编译或性能结果；finally已恢复Heart，第26点继续运行。确证错误是重复传入OMP/OPENBLAS/MKL变量；作业监督器已固定这些变量为1，调用者只允许PYTHONPATH/LD_LIBRARY_PATH扩展。移除非法扩展，保留监督器默认单线程协议；新增输出目录参数，使修复重跑进入独立v2目录，v1失败记录保留。
12. 修复后的launcher已实际启动，exec session6602；当前应优先检查 `three-dimensional-faster-performance-v1/edge-predicate-candidate-v2/events.jsonl`。它等待Heart point26终态后执行构建和对照，仍由上下文恢复dispatcher。不要把v1失去launcher进程当成v2无需暂停；不要独立恢复正在被v2拥有的dispatcher。
13. Heartbeat 11:33UTC（19:33北京时间）：v2独立构建退出0，15.503秒。共享边候选与完整CGAL相交构造oracle的488次精确对照全部通过；baseline和candidate各通过44个字面接触/绕向案例、三组各2556对的穷举/BVH对照及既有原生边界回归。候选面条实际边界计时已运行约8分钟，尚无终态；baseline实际边界计时尚未开始。没有加速结论或完整材料准入。dispatcher仍由v2拥有而暂停，Heart point26已完成；不得独立恢复或修改被候选作业绑定的源码。

## 跟踪规则

- 先检查 performance-v1/events.jsonl、launcher 进程、diagnostic receipt 与 profile/profile.json。
- 若仍等待当前 Heart 点或诊断运行中，不重复启动作业、不独立恢复调度器、不修改依赖。
- 诊断结束后先确认 Heart coordinator 已恢复，再审查栈分布及原始栈。没有符号或类别 UNRESOLVED 时保留未知，不能据此宣称热点已确认。
- 若 launcher 已终止而 coordinator 仍为 T，核验身份和所有大型工作进程已经结束，才恢复 coordinator，并保留故障记录。
- 科学结果仍按旧计划独立验收。采样诊断成功不把失败 stage04 升级为前驱，不意味着新的性能策略或面条输入已经通过。

## 后续依赖

- 根据实际栈确认热点，登记最小的性能修改假设；候选构建使用独立路径，不覆盖 Heart 当前模块。
- 在冻结状态上比较修改前后实体、质量矩、判定和成本；保留保存重载独立检查。
- 性能修复确证后继续面条阶段排除及完整输入准入。
- Heart 完成 81 点及分支汇总后运行真实优化；两个案例各需优化、独立材料检查、重复、合格 OBJ 及实际动力学记录。
- 正常求解、输入预处理、开发验证、终态验收分别计时；与同输入、目标和最终精度的体素路径比较。未测到加速前不授予性能完成。
