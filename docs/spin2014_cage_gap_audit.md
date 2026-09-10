# Spin-It 2014 cage 路径：本轮源码核查与后续门槛

核查基线：`13afcf0`。依据原始 SIGGRAPH2014 论文第6节、Figure4、第6.1节（本机 `/home/lin/下载/SpinItPaper.pdf`，已核查完整第6页）；[作者官方项目页](https://igl.ethz.ch/projects/spin-it/)用于核对版本。本文不授予任何完成状态。

## 已定位的差距

| 环节 | 当前源码证据 | 严格路径需要补齐 |
|---|---|---|
| 完整 cage | `AutomaticCageBinding.cpp` 调用 farthestPointSample；C++默认16个 handles，Python入口默认12个 | 直接从原始表面确定完整单元边界，所有边界角点成为handles；不得抽样 |
| 双层网格 | `GridDeformationDomain::fromVoxelGrid` 读取完整corners，分配完整cell表，再三重循环根格点 | 明确较粗cage层与较细BBW层；稀疏索引只存实际域，不能由低分辨率occupancy反推原表面 |
| BBW离散 | `BoundedBiharmonicWeights.cpp` 使用无权图Laplacian的LᵀL，逐handle解后截断并归一化行 | 与论文引用的BBW离散及边界条件做独立数值核对；现有方法不能未经证明直接标为论文等价 |
| 权重检查 | 当前diagnostics检查bounds、行和及Dirichlet约束；没有归一化后的独立最优性证据 | 保存原始解与归一化变化；验证实际最终权重对应的能量/约束/KKT，不把可行性当最优性 |
| 联合变量 | `solve_spin_it_2014_cage` 明确拒绝PAPER_GIVENS；旧CageSpinDeformationFunctional只有beta与handle参数 | 接入显式phi与完整相应等式约束、对所有变量验证导数 |
| 旧Python桥 | `bindings/pybind_module.cpp` 将已有cavity作为单个voidCell，fillLaplacian为1×1零矩阵，三个正则系数固定1e-3，minimumShellThickness为0 | 新严格桥必须保留实际优化叶节点及其各自beta、真实邻接/独立权重和输入壳厚，不能直接复用该旧入口 |
| 联合优化终止 | 旧CageSpinDeformationSolver使用自己的投影梯度/活动集路径 | 与已审计的SLQP可行性恢复、KKT及历史更新契约统一验证，不能继承其他求解器的测试结论 |
| 形状与几何 | 旧functional存在双拉普拉斯/位移能量和顶点质量矩导数 | 在完整cage与新单元上验证整个W映射、几何更新、翻面/自交/壳厚以及封存链路 |

这些是源码现状和验证缺口，不是对“任何图Laplacian都不可能近似BBW”的证明。论文未公开的离散控制参数仍须标注来源，不能推测成作者数值。

## 当前增量的边界

Task5g只补完整源表面cage几何。它为后续权重和优化提供完整角点集合，但不自行完成较细网格、BBW、Givens联合优化或Dancing Couple实验。旧API保留兼容，不获得新路径的验证身份。

## Task5g 已执行结果（2026-09-10）

实现提交：`b61993772bd9429b08425353af61a969aba6ce71`。新增原生 `SourceSurfaceCage`：从原始表面分类选定层级的相交单元，去除内部共面边界，保留全部边界角点。支持空腔、穿孔及多分量；奇异边界明确拒绝；预算耗尽不降低层级或抽样handles。旧抽样API保持不变。

- 全目标构建及当前已注册测试：200/200 通过，492.41秒。
- Faster关闭的独立配置：新增测试1/1通过，确认没有隐藏Faster依赖。
- 独立控制器导出：216个单元、218个handles、432个三角面；两次JSON逐字节相同。
- 独立Python有向四面体积分与固定解析包络对照：体积 `27/64`，质量矩/惯量检查最大绝对误差 `1.113567084231137e-19`。
- 稀疏level9/20测试验证真实层级的几何提取，不等于完成对应层级的优化。

原始日志和独立导出保留在 `build/task-5g-evidence-vL6KVn/`；源码、二进制及日志摘要见 `.superpowers/sdd/2026-09-05-spin-it-strict-computational-v2/task-5g-controller-verification.md`。这些解析夹具属于 `SYNTHETIC_REGRESSION`，不能计入论文原始案例完成数。完整复现状态仍为 `PARTIAL_BLOCKED`。

网格对齐时的闭单元接触采用明确的`IMPLEMENTATION_DEFAULT`。网格形成奇异边界时应拒绝，不应默默删handles、补孔或改分辨率。cage包络的体积和惯量只用于几何检查，**不是原模型的材料质量矩**。

## 后续可验收顺序

1. 基于完整cage构造较细的源表面权重域；记录两个层级、实际cell/node/handle数量、资源和哈希。验证新生成叶节点的三线性插值与共享边界一致性。
2. 核对并实现严格BBW离散和求解；解析/非对称测试、边界one-hot、partition-of-unity、能量及独立最优性检查；不静默以归一化掩盖求解失败。
3. 显式beta/phi/全部handle参数的联合functional和SLQP；检查链式导数、KKT、活动集变更、几何拒绝与回滚记录。
4. 完整变形材料重新积分和封存，然后接入展示及独立仿真。取得原始Dancing Couple姿态、域和权重前，原论文对照实验继续阻塞。

以上各项仍属于已批准的严格方案；没有移除缺失案例，也不以解析立方体代替原论文模型验收。

## Task5h：细网格与插值已实现（2026-09-10）

提交`4c6f89b`新增`SourceCageGrid`，修正提交`da1bd34`补齐浮点精度边界防护。新路径将实际粗cage单元均匀细分，使用稀疏索引和三线性场传递，不删除粗控制点。孔洞和多分量沿用已验证的cage域；外部点严格拒绝，不吸附到边界。此域是计算域，不是材料模型，也不是由低分辨率结果上采样冒充高分辨率优化。

实际独立检查：1,728个细单元、2,197个节点、6,084条边，保留全部218个粗handles；修正后二次导出逐字节一致且与普通坐标下旧输出相同，解析场误差最大`1.3877787807814457e-17`。最终输出在`build/task-5h-controller-fix1-GSRMx9/`；最初201/201与Faster关闭1/1记录在`build/task-5h-evidence-20260910/`，明确属于修正前结果。修正后两种配置各4/4相关回归通过，包含12组细网格测试；最终全量结果单独记录，不沿用旧源码通过状态。

独立审查发现并复现了大坐标情况下相邻细网格面重合的反例，随后复审确认已修正。当前要求细步长严格大于根区间最大ULP，因此也会保守拒绝某些原本可用的等ULP对齐输入；这是显式数值适用范围，不是模型物理不可行。不会自动改变坐标或分辨率。上游源八叉树的类似精度判断仍需另行统一审计，不能把本次防护推广为所有几何路径都已修复。

这完成了上述第1步的域/一般插值基础，但尚未给出实际BBW解。下一实质工作是权重的离散能量、完整cage面约束及求解器验收，再推进联合变形。细网格上的解析函数测试不计入论文原始模型成功数。

补充参考：[BBW作者官方项目页](https://igl.ethz.ch/projects/bbw/)明确将该方法描述为带界约束的Laplacian能量最小化。旧代码的bounds/行和检查只能建立权重可行性，不能独立证明它达到了这个优化问题的最优点。正式离散选择仍需核对原2011技术细节与2014网格使用方式；本次没有改动或替换旧BBW实现。

2026-09-10 后续核查补充：2011原文也公开了逐handle求解后归一化的加速变体，因此不能把“存在归一化”本身判为不忠实。必须分别验证对应原始问题及后处理结果，不能混用最优性声明；还需要cage面上的线性约束。具体源码证据及下一求解器验收条件见 [BBW契约核查](spin2014_bbw_contract_audit.md)。

最终修正后验收：`da1bd34` 的全部目标重新构建成功，已注册CTest **201/201通过，482.35秒**，测试ID无重复或遗漏，无NotRun/Skipped记录。原始输出保留在 `build/task-5h-final-AkTS5P/`；源码、实际二进制、独立导出和日志SHA见 `.superpowers/sdd/2026-09-05-spin-it-strict-computational-v2/task-5h-controller-verification.md`。独立复审已通过。本轮没有新增原论文案例或OpenGL结果；总体状态仍为 `PARTIAL_BLOCKED`。

## Task5j：FEM能量与完整节点约束已验收（2026-09-10）

新增`SourceCageFem`，不改旧BBW接口：完整细网格上的体积FEM、显式两种质量离散、`KᵀM⁻¹K`的矩阵无显式Q应用，以及所有边界细节点的约束。标准夹具保留218个粗handles，实际约束866个细边界节点。独立有理数对照、重复输出、FasterOFF3/3和全量202/202通过；独立审查无必须修复问题。详见 [FEM验收记录](spin2014_source_cage_fem.md)。

这补齐了上述第2步的能量/边界前提，**尚未完成实际BBW权重求解及KKT验收**。beta/phi/完整handles联合变形、最终变形材料封存、原始Dancing Couple等实验仍保持原有门槛，不因新增测试而升级。
