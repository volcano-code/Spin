# 三维 Faster 核心最终交付 Implementation Plan

> **For agentic workers:** 使用 superpowers:executing-plans 逐项执行，沿用本聊天直接实施方式；本文件不授权恢复案例或自动提交。

**Goal:** 在登记支持范围内完成正式 Cherchi 核心的正确性、工程链路和性能分别验收。

**Architecture:** 复用正式 Cherchi、原生 Newton 和同次材料封存链。先闭合当前候选证据，再冻结统一发布候选，最后执行三规模十种子与性能验收。已有不受影响的证据直接复用。

**Tech Stack:** C++、Python、Cherchi、CGAL 显式对照、精确几何与独立四面体积分。

**Spec:** 用户批准的“先正式 Cherchi，再修复几何”计划；执行事实见 docs/three-dimensional-faster-evidence-closeout-v2-progress.md。

## Global Constraints

- 心形、面条完整案例和自动跟踪保持暂停，保存失败仅用于回归。
- 不增加浮力、制造、视频、质量最小化、学习加速或 2014 五目标完成声明。
- 大型作业串行，24 GiB、swap 0、60 分钟，编译最多四并发，数值库单线程。
- 保留磁盘启动 100 GiB／运行 20 GiB 门槛；不删除历史、不提交、不合并、不覆盖旧模块。
- 原始残差 ≤1e-12，alpha=t_max=0.5；种子 0–9，每种子八次尝试，每次最多 60 次求值。
- 解析质量矩 ≤1e-10；适用生产导数三尺度 ≤1e-6；终态独立质量矩和约束 ≤1e-6。
- 不换失败种子、不答案初始化、不拼接超时结果、不静默换后端。

## 已核实基线

正式 Cherchi 小型闭环已有 17 次接受更新、终态保存和记录核验；尚非整体准入。三个 legacy 固定状态各 81 点数值与有限分支认证通过。新 projective selector 的 81 点作业正常结束，30 项方向／尺度检查通过，独立积分最大差 1.734723475976807e-18；material_branch_certification=false、core_admission=false。单次导数复用只有候选代码，未构建验证。30 种子和最终性能尚未执行。

## Review Focus

- 多分量内部见证冲突或零裕度：必须拒绝，不能多数投票。
- 坐标别名、合法极小面、真正退化：分别判定，不删面绕过。
- 构建宏与实际对象不匹配：禁止启用同实现导数复用。
- 后端、材料、Jacobian、日志身份失配：终态封存必须拒绝。
- 细分导致实体变化：不得进入“同实体三规模”比较。

## Task 1：关闭当前几何修复的有限证据

**Files:** tests/quadric_component_selector_regression.cpp；src/spin/QuadricMomentEvaluator.cpp；experiments/faster_core_repair_20261012/launch_component_selector_test.py、launch_component_selector_branch.py。

- [ ] 读取 component-selector-derivatives-job-01 的 result.json 和 source_side/report.json，记录数值通过但分支未通过的事实，不重跑 81 点。
- [ ] 串行执行 `python experiments/faster_core_repair_20261012/launch_component_selector_test.py`；一致符号通过，冲突与零裕度拒绝。
- [ ] 串行执行 `python experiments/faster_core_repair_20261012/launch_component_selector_branch.py`；从保存材料精确重放构造与有向平面链，逐模式分类。
- [ ] 失败只保存并诊断首个不能解释的模式；不得以质量矩相近代替实体／连接证据。
- [ ] 复核登记的封口、方向、细小面、真正退化与来源失败；逐项列修复或明确不支持，不能把后者记作成功。

## Task 2：消除已定位的重复导数

**Files:** src/spin/FullRankMomentEvaluator.cpp；experiments/faster_core_repair_20261012/build_single_transport.py、test_single_derivative_transport.py、compare_single_transport.py。

- [ ] 保留已观测 RED：精确材料配置中 transport 发生两次。
- [ ] 执行 `python experiments/faster_core_repair_20261012/launch_single_transport_build.py`，校验源 hash、编译参数、对象 hash 和真实链接对象一致后才启用宏。
- [ ] 使用新模块执行 `pytest experiments/faster_core_repair_20261012/test_single_derivative_transport.py -q`，要求调用数为一。
- [ ] 执行 `python experiments/faster_core_repair_20261012/launch_single_transport_comparison.py`；三状态各冷一次／同进程三次，实际材料、质量矩、Jacobian 与判定一致，调用数 2→1。
- [ ] 单独登记新构建科学身份；该复用只适用已证明同实现的 overlay，不转授标准 CMake 排序实现。

## Task 3：冻结统一正式候选并完成闭环回归

**Files:** python/make_it_stand/verification/cherchi_inprocess.py、cherchi_run_evidence.py、faster_material.py；现有正式入口与对应测试。

- [ ] 将通过的 selector 和复用构建冻结为同一个候选，记录实际加载路径、对象与模块摘要、请求设置、支持范围。
- [ ] 运行小型真实种子完整入口，要求实际接受更新、同次 Boolean 附件封存、保存后独立重载，不在终态重建另一份材料。
- [ ] 执行现有来源错误、旧 J、密度重复、跨输入缓存、后端失配、材料篡改、缺日志及无接受更新测试；任何漏检先补 RED 再最小修复。
- [ ] 对实际改变映射的路径补对应证据；数值不变的复用以同状态严格对照佐证，未受影响的旧批次不重跑。
- [ ] 最终证据必须指向统一候选，禁止用不同版本的局部通过拼成完整验收。

## Task 4：三规模与三十种子

- [ ] 从已有非零固定材料小型陀螺生成共享边细分的约万面／十万面版本，独立核验相同分片平面实体、封闭性、体积与质量矩，冻结输入摘要。
- [ ] 先运行每规模 seed 0 的正式链路冒烟检查；失败先定位，不批量启动其余种子。该次完整合规运行可纳入正式账本，不重复计算。
- [ ] 在同一候选上依次完成每规模 seed 0–9，保留全部尝试、回溯、重启、失败与资源终止。
- [ ] 每规模十种子全部取得合格终态才通过；每规模最小编号成功种子独立进程重复，比较科学事件、材料和摘要。
- [ ] 单独核验外表面身份、固定／保留／空腔完整合法分区、独立质量矩和约束；不以数值收敛替代材料验收。

## Task 5：性能与交付审查

- [ ] 分开统计输入准备、正常求解、开发验证、终态交付；正常求解含全部实际回溯与重启。
- [ ] 固定状态冷一次／同进程三次，同硬件、同输入、同检查范围比较；报告全部样本、中位数和最慢值。
- [ ] 约万面核心求值／完整求解中位数目标 ≤1／10 秒；约十万面 ≤5／60 秒。未达标明确标记性能未完成。
- [ ] 未达标时仅优化分段计时证明的热点；改变数值映射则使对应证据失效，禁止以减少检查范围冒充加速。
- [ ] 最后独立代码审查，覆盖上述五类 Review Focus，修正阻断项后仅复测受影响内容。
- [ ] 交付正式命令、版本／输入清单、论文步骤与实现映射、失败关闭记录、30 种子账本、重复及终态材料、正确性／工程／性能三张结论表。

## 预算与停止规则

证据关闭和复用验证 0.5–1.5 工作日；统一正式候选 0.5–1 日；规模基准、重复与性能 1–2 日；独立审查 0.5–1 日。剩余合计约 2.5–5.5 工作日，仅在未暴露新算法缺陷时成立。首两任务完成后用实际成本重估，不以到期为通过理由。

任何失败先记录阶段、输入和首个见证；无修复或新诊断假设不重跑同一长作业。核心验收报告提交后，心形与面条仍暂停，等待用户明确恢复。
