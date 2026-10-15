# 已有工作获取与首次接入验证（2026-09-17）

## 结论与边界

已获取七项公开参考源码，并真正使用 ProxQP 和 Cherchi Boolean 做了独立验证。
这不是两篇论文官方完整代码，也不是完整复现；项目仍为 `PARTIAL_BLOCKED`。
未修改原生优化器、未替换论文 profile、未创建最终 commit/tag，未覆盖主工作区和旧实验。

开发目录：`/data/xvxinyao/moment-guided-spin-v1/code/make-it-balance-paper-sim`。
下载目录：`external/reference-sources-20260917`，约 867 MiB。
原始压缩包、解压源码、许可证和来源均保留；压缩包/源码缓存单独 gitignore，
避免把其他许可证的完整上游工程批量纳入本项目 MPL 源码提交。

## 获取与使用状态

精确 URL、完整 commit、压缩包 SHA-256、许可证文件在下载目录的 `DOWNLOADS.json`。
其中 fetched/built/used 是**获取时快照**；本报告记录后续使用。

| 工作 | 锁定版本/提交前缀 | 本轮实际使用 |
| --- | --- | --- |
| ProxSuite / ProxQP | v0.7.3 / `b93d7778` | 项目专用 Python 环境安装同版本上游 wheel；真实 QP 两次独立求解 |
| qpOASES | 3.2.2 / `680f18c8` | 已获取，未编译/接入 |
| CppAD | 20260000.0 / `5d51b2aa` | 已获取，未编译/接入 |
| Ipopt | 3.14.20 / `1e71ba4e` | 已获取，未编译/接入 |
| Cherchi Boolean | `7bd6c2697695f555f9a15d7e154e084d841316f0` | 已编译独立 headless 适配器；真实交集/相减测试 |
| OCCT | 7.9.3 / `a016080bf6738d6aeae020badee4e888ad1540a5` | 已获取，未编译/CAD 接入 |
| Mirtich volumeIntegration | 上游 tar 的 SHA-256 锁定 | 已获取，未编译/接入 |

已有 libigl、CGAL 只记录版本，没有重复下载或修改。
Mirtich 的 `volInt.c` 源文件头声明 public domain/free use，并有版权署名；均原样保留。
只检查 README 会遗漏该声明，初始获取说明已据独立审查纠正。

独立复核及本地检查器确认：七个压缩包 SHA-256 全部匹配，44,641 个原始普通文件
与压缩包逐一字节一致，构建未修改这些上游文件。证据：
`build/reference-qp-20260917/verified-downloads.json`。

## Spin-It：真实 QP 的外部交叉验证

输入不是人工构造的已知解，而是已有封存 G7O7 重建茶壶轨迹中的 iteration 57
子问题：9,304 个填充变量及 1 个角变量，1 对有效 BFGS 历史。
它属于重建模型诊断，不属于作者原始案例。

输入文件：`build/spin2014-g7-repair-20260916/qp-512-repeat0.json`，SHA-256：
`0d7b8b36eb053cf28024a5d1bad84f2919a3de901d37b737265f159a68cbb67b`。
保留相同的梯度、线性化等式、全部有限盒约束及 BFGS Hessian；
不把保存的原生方向用作初值。ProxQP 从 `NO_INITIAL_GUESS` 开始。

| 指标 | 当前源码重复 0 | 当前源码重复 1 |
| --- | ---: | ---: |
| ProxQP 迭代数 | 9 | 9 |
| Hessian 组装 + 初始化 + 求解 / 秒 | 8.5760 | 8.5317 |
| 其中求解 / 秒 | 4.2053 | 4.1988 |
| 进程峰值 RSS / KiB | 3,923,320 | 3,925,284 |
| 独立 stationarity 残差 | 6.12e-17 | 6.12e-17 |
| 独立等式残差 | 3.05e-16 | 3.05e-16 |
| 最大盒约束违反 | 6.96e-13 | 6.96e-13 |
| 最大互补残差 | 6.47e-15 | 6.47e-15 |

两次科学输出 SHA-256 均为：
`ac1d293f8a598e4d0468a42a608e82d6dd6861c871951b1bb0d0e38663515365`。
与原生方向的最大绝对差为 `2.0581e-11`，原生方向亦重新通过独立 KKT 检查。
未放宽原 `1e-9` KKT / `1e-10` 盒约束与自由变量乘子门槛。

这确认了**这个固定子问题**的局部最优解一致性，不证明后续非线性优化、二值化、
几何封存或物理稳定性。包含稠密组装的外部总耗时高于此前原生回放约 5.58 秒；
原生耗时是历史参考、并非本轮同环境新鲜 A/B 基准，不能宣称获得加速。
更不能只比较 ProxQP 求解部分而忽略初始化和矩阵展开。

实现：`experiments/teapot_spin2014/external_qp_replay.py`。
当前证据：`build/reference-qp-20260917/proxqp-final{0,1}.json`、`verified-pair.json`。
初始两次运行因随后补强来源记录而被检查器识别为 `STALE_RESULT`，仍原样保留。
来源记录绑定输入、实际 Python、求解器扩展、53 个实际映射的共享库（含 MKL/OpenMP）、
源码、版本和线程环境；执行上限 300 秒、24 GiB 虚拟地址空间，实际内存另记 RSS。
预编译 wheel 的上游编译选项未自行重建，不声称与作者的软件环境相同。

## Faster：公开 Boolean 组件的实际构建和验证

适配器：`experiments/cherchi_reference_probe`，与现有优化器解耦。
调用下载的真实 `booleanPipeline`，使用其附带依赖，不使用 voxel/凸包回退。
固定 TBB 最大并行度为 1，不修改上游源码、不对输出进行排序或后处理来制造一致性。
适配器编译使用 GCC 13.3、C++20、Release、最多 4 个编译任务。

非对称相交盒体为 `SYNTHETIC_REGRESSION`：

- 交集体积 3，一阶矩 `(6, 3.75, 2.25)`。
- 相减体积 6，一阶矩 `(7.5, 5.25, 4.5)`。
- 全部二阶矩（含非零交叉项）与解析值及独立有向四面体积分一致。
- 交集两次原始 OBJ 哈希相同，相减两次原始 OBJ 哈希相同。
- 控制器在另一个新目录重新运行，结果及哈希再次一致。
- 3 项保护测试验证：已有输出不覆盖、已有实验目录不复用、未知操作拒绝执行。

初始多线程版本出现不同的原始网格哈希，历史结果和旧二进制全部保留于
`build/cherchi-reference-20260917`。当前单线程版证据位于
`build/cherchi-reference-singlethread-20260917`，其中 `REPORT.md` 有完整命令和哈希。

该适配器只用于受控夹具验证，尚非可接收任意模型的生产 Boolean 服务。
本次确定性只覆盖该夹具/构建；上游编译警告仍存在。没有验证复杂曲面类型、
Boolean 拓扑切换、离散交点导数、论文完整 Newton 求解器或作者定制版本。

### 同一离散二次曲面输入的 CGAL—Cherchi 交叉检查（2026-09-22）

此前“盒体通过”不能证明 Faster 实际使用的二次曲面离散与现有 CGAL Boolean
结果相符。为此新增受证据绑定的运行器
`experiments/spin_it_reproduction/strict/run_faster_reference_boolean.py`：它从当前
原生 Faster evaluator 生成五个**实际穿过设计域**的离散曲面（椭球、单双叶双曲面、
平面和抛物面），直接导出同一份设计网格和二次曲面网格，分别作为 CGAL 与固定
Cherchi 组件的输入。每个操作由原生 `quadricSolidIsPositiveRegion` 的封闭截断体
符号决定，不从参考结果倒推交集或相减。没有重采样、
体素化、凸包替代或几何后处理；OBJ 写出仅将数值相等的 `-0` 规范为 `+0`，不会改变
顶点数值、连通性或曲面离散。

两份独立运行证据位于：

- `build/faster-reference-boolean-evidence-20260922c-run-1/evidence.json`
- `build/faster-reference-boolean-evidence-20260922c-run-2/evidence.json`

各夹具两次 Cherchi 原始 OBJ 哈希一致；使用项目独立质量矩积分重新测量十项 raw
moments，相对当前 CGAL 输出的最大绝对差为 `3.5527136788005009e-15`（预登记门槛
`1e-12`）。严格审计器只据两份当前证据授予
`faster.geometry.boolean_independent` 这一**合成方法**项；它不授予任何作者案例、
CAD、完整曲面类型、拓扑切换或论文级几何结论。

离散椭球导数回归还新增了真实设计域顶点接触的拓扑临界测试：临界参数处
`derivative=false`，两侧的 topology/scientific identity 不同。该标记已写入两份
`build/faster-discrete-moment-evidence-20260922c-run-{1,2}/evidence.json` 的必需原生
转录中。因此，当前结果只在固定来源特征拓扑的开区间上主张导数有效；临界点和非椭球
封口不会被错误地重用为该导数。

### Faster：十项原始质量矩与系数对偶基础（2026-09-22）

新增 `experiments/spin_it_reproduction/strict/run_faster_moment_foundations.py`，在两个
独立的新目录中绑定当前原生回归程序。它核验十项 raw mass moments 中体积、一阶矩和
二阶（含交叉）矩的约定，以及二次多项式中 `xy`、`xz`、`yz` 系数与相应原始质量矩的
字面对偶；同时对由该系数到实际几何的局部微分走一条独立的原生回归路径。

两份当前证据为：

- `build/faster-moment-foundations-evidence-20260922e-run-1/evidence.json`
- `build/faster-moment-foundations-evidence-20260922e-run-2/evidence.json`

严格审计器只据二者授予 `faster.moments.raw` 和
`faster.moments.coefficient_duality` 两项**合成方法基础**要求。运行器明确拒绝授予
Faster Theorem 1／2、全局最优、布尔几何、作者案例、Table 或 Figure 的结论；这些都
需要各自的条件检查和原始资产证据。

### Faster：Equation 7 的最终材料数值条件审计（2026-09-22）

新增 `spin::auditFasterTheoremOneNumericalConditions`。该检查器从最终几何二次曲面
恢复论文的 literal monomial 坐标（交叉项先乘二再归一化），独立检查：单位长度、等式
约束残差、严格正的 \(\mu\)，以及
\(\nabla_r f=-\mu p+\operatorname{Jac}_r(g)^T\lambda\) 的有限精度残差。非有限或
维度不一致输入返回 `UNKNOWN`；有限但违例的输入返回 `VIOLATED`，不会静默重求乘子。

**canonical FullRank 严格路径**把该审计结果封入输出身份，并在
`finalGeometryVerified` 前强制要求 `SATISFIED_NUMERICALLY`。单元夹具、真实离散材料
integration 与重复 replay 均已通过。两次重新绑定的支持性证据为：

- `build/faster-theorem-numerical-conditions-evidence-20260922d-run-1/evidence.json`
- `build/faster-theorem-numerical-conditions-evidence-20260922d-run-2/evidence.json`

它们分别运行条件拒绝夹具和 canonical 非平凡材料 replay；审计器将其列在
`supporting_evidence`，其效果固定为
`SUPPORTING_ONLY_NO_CATALOGUE_REQUIREMENT_VERIFIED`。因此它不会把 Theorem 1／2 的目录项
或任意作者案例改成 `VERIFIED`。

这个模块的范围刻意固定为
`NUMERICAL_EQ7_AND_EQUALITY_CONDITIONS_ONLY_NOT_A_THEOREM_PROOF`：它不证明特征函数和
连续正区域在测度零集之外的一致性，不证明无限维局部极小值，也不处理 Theorem 2 的
\(p=0\) 退化分支。因此严格目录中的 Theorem 1／2 需求不会仅因该数值审计而被授予通过。

同日的重复运行诊断发现，未 canonicalize 的旧 FullRank 兼容路径会直接以
`mesh_boolean` 返回的顶点／面顺序累积导数；即使值几何不变，该顺序也会造成约
`1.65e-14` 的 Jacobian 差异并漂移 Newton 的事件哈希。浮点舍入模式在调用前后不变，
因此不是原因。已有 canonical 路径以不焊接、不舍弃面、不取整的确定性重排材料网格和
birth-face 对应关系后再微分；同一进程五组完整 replay 均通过。严格证据只能使用后者；
旧路径保留作兼容诊断，不能用于确定性发布门。

### Faster Equation 14 的原生 Newton 协议证据（2026-09-22）

新增 `experiments/spin_it_reproduction/strict/run_faster_eq14_newton.py`，两次独立
执行当前原生 `quadric-hollowing-solver-regression` 的四个保护模式：默认完整球面
Newton 回归、`--raw-stopping`、`--dual-pairing` 与 `--dual-coordinates`。它们共同
检查 Equation 14 的原始残差停止语义、双单项式坐标配对、每步乘子、球面测地线步长
上限 `0.5`、回溯因子 `0.5`、每次重启至多 60 次函数评估，以及失败首次尝试仍要保留
并继续重启的事件计数。

初始证据保留在 `20260922a` 目录。由于上述对偶回归的二进制摘要已更新，当前源码的
重新绑定证据为：

- `build/faster-eq14-newton-evidence-20260922e-run-1/evidence.json`
- `build/faster-eq14-newton-evidence-20260922e-run-2/evidence.json`

严格审计器重新读取二者、绑定二进制与 CMake/源代码摘要，并比较四个模式的科学输出。
它只授予 `faster.eq14.newton`、`faster.eq14.geodesic`、`faster.eq14.multiplier` 和
`faster.eq14.restart` 四项**合成方法**要求。篡改任一模式转录会得到
`STALE_RESULT`。这既不证明全局最优，也不证明任一作者输入、Table 1、Figure 或 CAD
案例已经复现；作者资产与关键参数仍是独立阻塞项。

### Faster：实际封闭网格浮力几何（2026-09-22）

新增 `spin::HydrostaticGeometry` 与 Python 的 `compute_hydrostatic_geometry`：对严格封闭、
外向的材料网格，以实际水面负侧做 CGAL 网格交集，并从所得封闭边界积分排水体积、
浮心、水线面积、水线质心和水线中心二阶矩。它拒绝没有真实交线或源三角面与水面共面的
退化输入，避免把手填水线数据伪装成几何测量。

随后新增 `find_hydrostatic_equilibrium`：它只在调用者明确给出**闭合排水边界**、固定方向与
目标排水体积时，用实际网格切割的体积残差求一个静态水面。该接口明确不推断空腔是否进水、
不求物体姿态，也不证明浮力稳定；这些物理语义仍必须由正式案例的原始组件定义。

立方体的水平、平移和旋转水线解析夹具、目标排水量为 1/2 与 3/4 的静态水面、以及无水线／
全体积目标拒绝路径均通过。两份新格式严格证据：

- `build/faster-hydrostatic-geometry-evidence-20260922b-run-1/evidence.json`
- `build/faster-hydrostatic-geometry-evidence-20260922b-run-2/evidence.json`

审计器只授予 `faster.hydrostatics.actual` 这一合成**方法**项。该结果不是 Floating Island
或 Tea Infuser 原实验：原始几何、可优化域、密闭/进水部件解释、姿态和参数仍未取得，
因此对应案例与最终复现声明继续保持阻塞。

## 实际测试与复核

最新控制器测试：**83 个 pytest 测试通过**，另有 **3 个 Boolean 安全测试通过**，
以及上述两种 Boolean 运算的重复解析检查。无 skip。
这是本轮相关测试组，不是全仓库测试或完整论文验收。
先记录失败测试，再实现/修正，再执行通过测试；失败日志没有删除。
独立只读审查已核对数值映射、来源、获取文件及修正，不用手填布尔值授予论文完成。

外部 QP 依赖组位于 `tests/reference_qp`，不会强制核心运行环境安装 ProxSuite。
wheel 锁定文件：`experiments/teapot_spin2014/reference-qp-requirements.txt`。
测试命令（从开发目录执行）：

```bash
PYTHONPATH=python PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
build/reference-qp-20260917/venv/bin/python -m pytest -q \
  tests/reference_qp/test_external_qp_replay.py \
  tests/python/test_native_qp_replay_audit.py \
  tests/python/test_slqp_direction_replay.py \
  tests/python/test_independent_mass_moments.py \
  tests/python/test_mass_moment_slivers.py
```

Boolean 重跑时必须提供新的、尚不存在的结果目录；不要复用已封存目录。

## 下一优先级与仍未解除的门槛

1. Spin-It：剖析更多真实 QP/SLQP 步骤的组装和线性代数成本；利用稀疏盒约束、
   低秩 BFGS 做数学等价的性能优化，以本轮外部证书作交叉检查。之后再跑完整 G7，
   仍需连续收敛、自适应终止、二值验证、几何封存和仿真，不能只提高预算。
2. Faster：已完成一个真实裁剪离散椭球上的 Cherchi—CGAL 同输入交叉检查，下一步
   仍需覆盖其他曲面类型与拓扑切换，并检查完整离散几何—质量矩映射的导数；一个
   合成夹具不能代替这些项。
3. 已下载的 Mirtich、CppAD、qpOASES、Ipopt、OCCT 可逐项接入；本轮不把“下载完成”
   写成“全部模块已实现”。
4. 两篇论文的作者原始资产、关键姿态/参数、完整案例及图表等缺口没有因依赖下载而消失。
   现有 G7 重建茶壶仍未得到新的收敛模型；没有新的论文原案例仿真视频。

没有发布 `VERIFIED_COMPUTATIONAL_REPRODUCTION`，没有降低原严格复现标准。
