# 三案例核心算法优先 — 2026-10-12

## 最新用户指令优先于原深度验收调度

用户要求去掉不必要额外检查，先使用现有 Spin-It Faster 核心算法取得心形、面条结果。开发导数81点、逐点独立Nef重建、有理材料整套认证与性能benchmark不再阻塞核心求解。Newton raw residual 1e-12、level2、Faster目标、输入模型、单位密度、seed0/1、8次尝试/60次求值仍保留。

历史所有输出保留，未独立验收的核心结果不授予完整材料、制造或论文复现状态。

## 已执行的调度切换

- 已停止旧 Heart dispatcher PID1766941（先前完整完成points00–26，81点集合不完整）；不要恢复或重启旧队列。
- 已停止额外 `edge-predicate-candidate-v2/candidate-audit-job-01` 单元 `spin2014-job-d32e28c1f1b84ee5accf3a66b3a7ae93.service`。该操作是用户改变优先级后的中止，不是数值算法失败。无性能比值结果。所有编译、回归和部分计时输出保留。
- 新核心入口 `core_faster_heart.py` 复用 `rounded_solver_journal.py` 与现有外部实际质量矩映射 Newton，保留耐久初值/失败/迭代记录。显式使用binary64源面、CGAL与实际离散Jacobian；不再生成高精度有理源面及Nef权威材料。原真实固定材料积分与冻结请求继续复用。

## 心形：实际快速运行及确证修复

- `rounded-clearance-v3-core-faster-v1`两个seed实际完成。seed0为129.417秒、seed1为186.783秒；各8次初值求值，0接受更新，状态GEOMETRY_EVALUATION_FAILED。不能作为优化成功。
- 首个失败状态进行一次已登记假设核对：`rounded-clearance-v3-core-fixed-state-job-01`仅改用现有精确binary64相交谓词（不是Nef材料认证）。69.406秒后实际材料和生产Jacobian均返回，derivative_defined=true。证明该状态被旧近似相交查询误报，不能泛化为所有状态都无错误。
- 核心配置因此改用 `exact_binary64_self_intersection=True`，保留闭合、正体积和导数有效性等必要拒绝；仍关闭rational_source_coordinates、compact_rational_source_coordinates、exact_binary64_material_boundary。策略回归观察RED→GREEN，2/2通过（另一次从错误工作目录执行的discover命令仅发生导入失败，未计为通过）。
- 修复后完整核心运行已启动：`experiments/faster_heart_contact_20261011/rounded-clearance-v3-core-faster-v2/events.jsonl`，launcher exec session53844。不要修改绑定的worker/journal/请求/native等源文件，不重复启动大型作业。
- 数值成功才生成实际终态内件、固定件和装配预览OBJ；装配预览是同坐标组件组合，不是已合并/独立审核的打印实体。预览basic检查和完整审计状态均如实保存。

## 面条：核心预处理与表示问题

- `noodles-dense-v5-core-preprocessing-v1`沿同一inner-ellipsoid与固定81部件顺序，CGAL Boolean后仅进行快速有限坐标、闭合、绕向、正体积、零面积检查，不逐步重建Nef。
- 约11秒完成22部件；102.778秒完成72部件。第73个（tube-070-fixed）在105.848秒被4个binary64零面积面拒绝。输入构造显著变快，但全设计域未完成，不能开始优化。
- 新增只合并完全相同binary64坐标、移除合并后重复索引面、核验受影响顶点link的表示候选。不移动/近似吸附顶点，不删除三个不同顶点的共线面。3项真实mesh回归观察RED→GREEN并通过。
- `noodles-dense-v5-core-preprocessing-v2`从保存的72部件核心前缀重放同一第73次Boolean；该表示候选仍被CORE_BOUNDARY_NOT_CLOSED_ORIENTED_POSITIVE拒绝。候选未获采用，不把删除零面后的开口网格继续送入优化。
- 下一步仅定位这个固定失败Boolean，保存原有理输出及binary64派生物，检查坐标别名、相反绕向重复面、局部边/link；经零位移等价或已登记有界表示证明并通过闭合检查后继续剩余9部件。不得改变管束外形、轴、权重或删除固定部件来获得成功。

## 下一次跟踪应做的事

### core-v2终态与身份连接修复（12:44 UTC跟踪）

- core-v2两seed已结束：519.792秒与576.425秒，各8次求值、0接受更新，GEOMETRY_EVALUATION_FAILED。进程正常退出不代表数值成功。
- 实际保存的evaluation-000001返回合法材料和有效Jacobian。Newton回调拒绝原native后端字符串中的分号，错误为external evaluator map/backend/derivative identities must be nonempty and delimiter-safe。这是连接缺陷，不是可行性判定；另有一状态存在细小面退化拒绝，继续保留。
- rounded_solver_journal现在将完整原后端来源UTF-8字符串SHA-256编码为安全身份；原来源仍保存于每份实际evaluation，不更换后端或改变科学数值。新增身份测试观察失败后通过；保存的实际第一求值在同seed、一次求值诊断预算下通过真实native Newton解析，无需再做Boolean。4项相关测试通过（中间测试曾因回调tuple无copy导致失败，修正测试类型处理后通过）。此重放不是完整求解或初值答案。
- 已启动rounded-clearance-v3-core-faster-v3，同原请求、seed0/1、原预算，仅包含上述连接修复。运行中不得改journal/worker/绑定源文件，不重复启动大型作业。下一次优先检查v3而非已终止v2；旧v2结果保持原样。

1. 优先检查Heart core-v3事件、数值结果、进程。它正在实际求解，不是81点开发验证；正常进度保持安静。
2. Heart作业运行时只准备面条定向修复，所有大型实验串行。Heart终态后执行面条固定失败重放与最小表示修复，完成全域、真实固定材料和原目标请求后调用同一核心求解器。
3. 任何核心成功可先交付真实材料视图/OBJ及数值摘要，并标明独立审计待完成。物理接触输入仍需基本合法性和对应真实质量、质心、完整惯量，不能把无效mesh强行用于仿真。
4. 只通知实际结果、明确失败或需用户选择；不因普通检查计数增加而通知。不要再启动旧81点、Nef性能长对照或恢复已取消dispatcher。

### 13:14 UTC跟踪

- core-v3 seed0仍在运行，已记录12个接受事件、25次求值；最新试探被原下降规则拒绝，raw residual约0.62187，未收敛。身份连接修复已使实际Newton进入优化，不能把正常接受更新当作完整成功。
- 已准备capture_noodles_failed_boolean.py，保存第73次实际Boolean的原有理输出、binary64派生物、坐标别名、重复三角形方向及异常边。仅检查Python语法；尚未执行，不宣称发现了重复面或修复成功。Heart占用串行作业槽时不启动该Boolean重放。

### 13:44 UTC跟踪

- Heart core-v3 seed0仍在原60分钟作业内；不修改冻结求值器或扩大预算。
- 已启动launch_noodles_failed_boolean_capture.py，exec session8688。它按PID872658及进程起始身份等待现有Heart launcher退出，然后只执行一次已登记的第73次失败Boolean采集；目录noodles-dense-v5-core-failed-boolean-capture-v1。当前WAIT_EXISTING_HEART_CORE_V3，不存在并行大型作业。采集仅保存实际几何与局部证据，不授予域准入，不主动修复或删除面。
- 下一轮同时检查Heart v3终态和上述capture的事件、diagnosis.json及原始输出；不重复启动采集。

## 新批准核心修正优先

用户已批准docs/superpowers/plans/2026-10-12-three-dimensional-faster-core-repair.md；优先读docs/three-dimensional-faster-core-repair-progress-20261012.md。Heart v3 seed0已按60分钟RESOURCE_LIMIT_WALL终止，不续拼、不无修改重启；Noodles采集因ndarray序列化失败，旧有理文件是部分输出，不能使用。采集序列化已修正，待核心串行空位登记修复重跑。

核心新候选已对固定Heart新参数试探取得72.608→8.005秒（约9.070倍）且完整生产输出逐JSON相等。它是验证/边界证明复用修复，尚未完成收敛修复、后端实体对照或完整优化。候选位于experiments/faster_core_repair_20261012，旧模块保持原样。
