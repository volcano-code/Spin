# Joint material batch execution ledger

Plan: user-approved Spin-It / Spin-It Faster shortest completion path, batch A.

Existing worktree and all unrelated dirty changes are retained. No commit or
integration is performed while the current uncommitted work remains unaudited.

The frozen asymmetric fixture uses its existing geometry, 98 cage handles,
INDEPENDENT_NORMALIZED actual BBW, density 1.3, objective and five weights,
frame (.5,.5,.5), beta linspace [.2,.65], and phi .17. The old four-iteration
regression remains unchanged. A separate experiment protocol uses the existing
solver defaults (120 iterations, 24 line-search steps), four adaptive passes,
the existing resource caps, and unchanged scientific tolerances. These are
implementation defaults, not recovered author parameters. Each new run starts
from the initial state and writes input and controls before invoking a solver.

Pre-flight finding: the existing binary material admission requires a privately
constructed completed adaptive endpoint. A nonconverged run cannot be promoted
to a candidate or sealed artifact. Preserve this gate and report the failure.

## 本批结论

状态仍为 `PARTIAL_BLOCKED`。本批没有修改求解器目标、容差、几何防护或
舍入规则，也没有获得新的作者原始案例完成记录。

| 工作 | 新增的真实证据 | 尚未通过的门 |
| --- | --- | --- |
| A：完整 cage 联合形变 | 从初始状态运行完整尝试，保存每个试探及几何拒绝原因；两个独立进程逐文件一致 | `RESTORATION_STALLED`，没有完成连续自适应流程，因此二值材料 `NOT_ADMITTED` |
| B：双材料几何 | 实际壳和内件分别导出，Python 独立积分；两次导出逐文件一致 | 重建 Sumone 内件为空；没有获得非平凡 Sumone 优化成功或合格封存 |
| C：Cuboid 舍入诊断 | G4/G5 逐单元、逐层级十项质量矩误差和约束误差分解 | 两组二值残差仍不通过，未启动 G8/G9 |

测试通过只说明实现与上述预期失败的记录机制正确，不能把失败案例算作
成功的论文复现。A 的连续量积分不是二值材料几何验收；B 的独立积分
不是完整材料重叠、自交或壁厚检查。

## A：联合形变完整尝试

复用 `SourceCageAdaptiveSolver` 和 `SourceCageBinaryMaterial`；新增的
`source_cage_joint_batch.h` 是原回归程序的实验适配层，不是替代求解器。
`run_source_cage_joint_batch.py` 绑定实际二进制、扩展、有效源码、依赖、
初始输入、输出和检查器。无最终准入凭证时，不强制导出二值对象。

当前有效重复：

- `build/source-cage-joint-material-20260923-f/`
- `build/source-cage-joint-material-20260923-g/`

两次均从初值开始，运行约 15.8 秒，生成时新鲜度 `CURRENT`。
11 个 native 科学文件逐字节一致，汇总 SHA-256：

`99613001d49a22f4dc792336d113b401c26641993ba9cf702f26dea2d177afd7`

| 指标 | 初始 | 终态 / 变化 |
| --- | ---: | ---: |
| 总目标 | 2.093136365485859 | 0.524952409318367 |
| 最大绝对约束（原函数尺度） | 4.8654358349533395e-4 | 1.6527401991271541e-4 |
| 填充变量 L2 变化 | — | 0.10162001769335599 |
| 全 handles 参数 L2 变化 | — | 2.0939044669844167 |
| 连续材料独立积分误差 | — | 1.2810265668751805e-16 |
| 外表面顶点 RMS 位移 | — | 0.2967419551397852 |
| 外表面最大位移 | — | 0.35313535168533977 |
| Hausdorff 下界 / 上界 | — | 0.28944878963497084 / 0.35313535168533977 |

Hausdorff 下界来自顶点到对侧三角面采样；上界来自相同三角拓扑的
对应顶点位移。没有把采样值声称为精确双向 Hausdorff 距离。

第一自适应 pass 的 14 条迭代记录包含 13 个正接受步，有实际填充和
handle 更新，但最终
KKT stationarity 为 `0.01667531630691544`，未满足收敛条件。186 次
guard 调用中，拒绝阶段记录包含：

- 163 次 `primitive not subset of actual outer material`；
- 6 次材料质量非正；
- 1 次 handle 尺度非正；
- 1 次外表面精确边界审计失败。

163 次包含关系失败的首个见证均是 primitive 0。检查器遇首个失败即
停止，所以这不证明后续 primitives 合法。最后三个正接受步从
`9.5367431640625e-7` 缩至 `1.1920928955078125e-7`，已经接近当前
线搜索最小试探；下一诊断应针对该见证及其变形外表面，而非继续盲跑。

最终目标线搜索与恢复线搜索各 24 次试探均无效，最小试探步长
`1.1920928955078125e-7`。最后保留的合法状态与试探失败必须区分。
这是当前离散表示与几何防护下的数值停滞，不是不可行证书，也未证明
存在积分实现错误。所有 callback、guard 对应关系、分区、rebind 和
二值准入原因均保留，后续可以精确重放问题状态。

输出 `last_accepted_outer_DIAGNOSTIC.obj` 仅是最后合法外表面；
不存在 `binary_material.obj`。因此本批没有合格的变形后二值材料，
不得启动其正式物理仿真或用它替代论文 Dancing Couple 结果。

同输入 hollowing-only 对照为 `FIXED_PARTITION_DIAGNOSTIC`：固定初始
分区与 handles，仅优化 beta 和角变量，收敛到目标
`2.090088159334087`，趋向恢复实体材料。它不是完整自适应 hollowing
对照，也不是非平凡挖空成功案例；这一比较要求仍未完成。

## B：两个材料区域的真实网格与独立积分

复用原双材料回归，新增 `--export-directory`；严格运行器新增
`--export-geometry`。没有修改原测试的密度、输入或停止条件。
输出依实际 label mask 提取每种材料的定向边界，空材料用 `mesh:null`
明确表示，不制造空 OBJ 或把期望标签计数抄入结果。

当前重复：

- `build/spin2014-dual-geometry-20260923-b/`
- `build/spin2014-dual-geometry-20260923-c/`

六个材料文件（两个 manifest、三个实际 OBJ、一个预登记协议）以及
`independent_geometry.json` 两次逐字节一致，生成时新鲜度 `CURRENT`。

| 输入 | 实际壳 / 内件 | 独立组合质量矩最大误差 | 判定 |
| --- | --- | ---: | --- |
| 4³ dense regression | 56 / 4 单元，密度 1.17 / 8.1 | 1.7763568394002505e-15 | `TWO_MATERIAL_GEOMETRY_INTEGRATED`，仅方法几何验证 |
| 16³ reconstructed Sumone | 192 / 0 单元，密度 1.17 / 8.1 | 1.1102230246251565e-15 | `NO_RETAINED_DENSE_INSERT`，保留失败 |

关键几何摘要：

- dense 壳：`f0afd6ddb8f3c4539e2990e92dbf890320bc2be0e82c221f36ead9bd5ed4c0a1`
- dense 内件：`92f198ec52f52bd398e4d5126d7367fc287e060b989df07e4186f8c4af765cd6`
- Sumone 壳：`f4aa51ee703410a8622f9d16dc087f85973308ba488721c8d00a0576d67feb9b`

独立路径读取导出的 OBJ，并重新计算质量、一阶矩、二阶矩、质心和
完整惯量。改变壳网格顶点的反向测试会失败。实际空内件的 native
减法残差约 1e-15 被原样保留，按既定 1e-10 门槛检查，不强行改成零。

重要限制：Sumone 继承的旧二值门槛为 `0.005`，已显式写入 manifest；
其旧 `CONVERGED` 不能升级为当前严格 `1e-6` 通过。dense fixture 仍是
合成几何，新增的是导出和独立积分证据，不是非对称原案例成功。
完整重叠检查、约束独立复算及制造封存仍需补齐。

## C：历史 Cuboid G4/G5 舍入误差

新增 `rounding_moment_breakdown.py`。限定读取不超过 64 MiB 的现有
receipt，不展开大型轨迹，不重新运行优化，也不宣称独立重算历史几何。
按 `beta >= 0.5` 原规则，对每个叶节点计算：

`delta_s = -density * (beta_binary - beta_relaxed) * leaf_geometric_moments`。

将十项增量变换到原登记目标坐标系，并逐层、逐单元分解所有 yo-yo
约束，包含 Givens 角。边界 beta=0.5、NaN 和不符合舍入规则均有测试。

| 历史输入 | 二值归一化残差 | 重组质量矩误差 | 重组约束误差 |
| --- | ---: | ---: | ---: |
| G4 | 2.558653046689352e-5 | 4.135903062765137e-20 | 5.066460261753961e-19 |
| G5 | 1.4472403638176351e-5 | 4.963083675318165e-19 | 8.918040979087328e-19 |

结果：`build/cuboid-g4-rounding-breakdown-20260923.json` 与
`build/cuboid-g5-rounding-breakdown-20260923.json`，均有原 receipt SHA-256。
两组主要约束偏差与舍入十项增量一致，差分恒等式误差不超过 3.5e-18。
这一检查排除了这些保存数据中的质量矩求和或目标坐标复算不一致；
它不证明整个自适应算法没有缺陷，更不是连续/离散不可行证书。
没有修改容差、执行后处理补救或重复相同失败优化。

## 测试、审查与证据边界

相关五项 CTest 实际执行通过（18.92 秒）；包含新增三项注册测试，
覆盖协议初值锁定、禁止覆盖已有运行、虚假准入拒绝、从几何重算轴约束、
截断/空文件保留、完整失败元数据、真实双材料积分和舍入端点。

独立代码审查后的已落实修正：

1. 空 stdout / 截断 JSON 不再阻断原始证据归档；归档精确保留原字节。
2. pass 的 rebind、终止及 binary 原因进入证据，不只保存终态 token。
3. 实际执行的 run-evidence / canonical helper 纳入有效输入摘要。
4. 双材料测试显式使用项目根目录，避免 CTest cwd 导致找不到既有资产。

完整回归 **281/281 通过**，exit 0，903.24 秒；独立解析日志得到 281 个
不重复的已执行测试 ID，均为 Passed，未发现 `Not Run` 或 unittest
跳过记录。本次构建启用精确几何、2014、Faster 和 Python，viewer 构建
关闭；不把这个测试集合声称为尚未执行的其他构建配置也通过。

最终日志：`build/spin-joint-dual-batch-ctest-20260923-b.log`。
SHA-256：`fd15ded230a2964f8daf75018bce4c642d9db9c59252c12540a1037b3f185cf3`。
本次套件包含 40 项 Faster 标签测试，均通过。systemd 报告整套任务
约 1.9 GiB memory peak、swap peak 0；不是高分辨率优化资源基准。

首轮在 systemd 环境中
因 PATH 不含 Anaconda 的 Ninja 导致依赖配置测试失败，约 12 秒时主动
停止并保留 CTest 临时日志
`build/spin-joint-dual-batch-ctest-20260923.log.tmpef476`（被停止的 CTest
未将临时文件重命名为最终日志）。确认 Ninja
实际位于 `/home/lin/anaconda3/bin/ninja` 后，显式固定 PATH 重新运行；
未修改测试预期或源码以绕过该错误。最终四份实验 evidence 在实际加载
声明的 `_make_it_stand` 扩展后重新检查均为 `CURRENT`，差异列表为空。

成功的进程退出只代表捕获完成，不能授予科学成功。当前记录绑定实际
加载文件，但 `build_source_attestation` 仍非已验证的可复现构建证明。
已有脏工作区保留，没有最终 commit、tag 或发布。更早 a–e 运行目录
原样保留，不将其旧新鲜度升级为当前版本证据。

资源策略：24 GiB、swap 0、3600 秒 systemd 上限，数值库单线程，实验
依次执行。快速单作业的 systemd memory-peak 部分记录显著偏低，不把它们
当作实际 native 峰值；因此也不能据此推导 G8/G9 内存需求。

## 下一可执行工作与未完成要求

1. A：利用已保存试探状态，检查 PrimitiveSubset 的具体失效单元及外面
   几何关系；只有发现实现不一致才修复。补全真正的自适应 hollowing
   比较和成功后的二值变形材料独立检查；不靠换权重或缩减 handles。
2. B：补材料间重叠/域归属及独立约束检查；Sumone 内件全挖空维持失败，
   不能因此调密度或目标。已完成的等密度退化检查继续作为回归。
3. C：登记源表面 Heart G8 资源预估，再决定运行；利用实测估算 G9，
   未运行不得获得 Level-9 优化能力声明。
4. Armadillo 设计域、Faster 秩亏类型及完整应用集成、作者原始资产/
   参数、正式论文矩阵和最终封存/仿真仍未完成；本批没有提升这些状态。
