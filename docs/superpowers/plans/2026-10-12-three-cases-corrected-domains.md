# 三案例：纠正设计域、真实优化与论文式交付 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement task-by-task. Checkboxes record evidence, not intended success. Do not automatically commit, merge, delete, or overwrite historical results.

**Goal:** 花盘、带柄心形和面条各取得真实优化、独立检查、同版本重复、合格 OBJ 和同对象物理对比，并采用用户参考图的表达形式展示。

**Architecture:** 输入构造、求解、材料验收、展示和物理分别设门。花盘采用现有二维圆边界/dogleg 后三维挤出；心形和面条采用 Faster 三维二次曲面求解及其惯量—高度目标。三例不借用2014五目标线的成功记录。

**Tech Stack:** 现有 Python/native、CGAL、GMP/MPFR、精确材料与独立积分、PyBullet、Blender。

**Spec:** 当前用户要求“三例均优化成功，呈现论文类似效果”；沿用已批准三案例计划，纠正2026-10-12实际OBJ核查暴露的花盘固定中间层问题。参考为用户提供的三张视频截图，不视为作者CAD或完整参数。

## Global Constraints

- RECONSTRUCTED_METHOD_VALIDATION；不声称原作者输入复现，不保证任意输入成功。
- 成功必须是真实接受更新、目标及材料变化、独立约束通过、材料保存重载和同版本重复；漂亮渲染不是成功。
- 不硬编码右下切口、预制双曲面空腔、旋转结果冒充期望方位或按优化成败调形状/权重。
- 三维显式CGAL、level 2，原单位/密度/目标权重，种子0/1，各最多8次尝试、每次60次求值；alpha=t_max=0.5，raw residual≤1e-12。二维保留自身dogleg协议。
- 解析矩≤1e-10；适用导数三尺度≤1e-6；独立终态质量矩及约束≤1e-6。
- OBJ内部边界误差≤1e-6L、对称差体积≤1e-6L³，L为原输入bbox对角线；自身分区、拓扑和包含合法。制造资格不自动通过。
- 每作业24GiB/swap0/60分钟，编译≤4、数值库单线程，大型作业串行；磁盘低于20GiB安全停止。
- 三例全部通过才关闭本计划；失败保持具体阻塞，不改分母、不替换种子、不降低门槛。

## Review Focus

1. 花盘中间孔边被额外固定会阻止真正贯通切除；输入测试必须检查固定材料的空间角色，不能只核验标签。
2. 把 cavity/design/original 与 final OBJ 一起装配会填回缺口；最终包必须与诊断包隔离。
3. 心形局部pin证书不能替代整模型连接、间距及分区；整体准入检查必须保留。
4. 面条增加采样不能消除中心线曲率导致的实体折返；单管及并集均须检查。
5. 动画、物理、惯量与OBJ必须绑定同一终态，不能沿用旧花盘物理或将夹持组件计入转子。

## Task 1：建立新版花盘输入，而非事后挖空

**Files:** 新建 `experiments/quadric_physics_demo_20261011/paper_flower_fixture_v2.py`、`test_paper_flower_fixture_v2.py`；复用 `paper_flower_fixture.py`、`paper_flower_capacity.py` 和现有二维入口。旧v1不修改。

**Consumes:** 冻结未优化轮廓、48孔及轴承装配。

**Produces:** 独立v2输入、角色清单、准入记录；不携带旧切割圆。

- [ ] 增加失败测试：中间原始花形实体应全部归入设计域；不得自动沿全部孔边侵蚀产生全高度固定骨架。上下固定层仍使用原轮廓，轴承座/外圈单独积分，夹持件不计入转子。
- [ ] 新候选明确采用“上下固定盖层＋可优化中间花形层＋独立中心轴承装配”。取消旧 `original.buffer(-RIM_WIDTH)` 引入的中间孔边固定网；允许中间固定网为空，不输出伪造非空网格。
- [ ] 使用初始几何核验中心连接、完整分区、层间接触、轴承座及质量矩；不把此项称为已验证打印强度。若机械连接需要额外固定结构，先说明具体几何依据，另登记输入，不根据解的位置添加。
- [ ] 生成优化前正面、侧面、三层展开和角色叠加图；固定+x向右、+y向上，说明新旧输入差异并确认结构后冻结。
- [ ] 重新计算容量界和物理量，重新执行两个种子及独立重复；旧圆心、残差、物理轨迹不能转授新版。

**Tests:** `python -m unittest discover -s experiments/quadric_physics_demo_20261011 -p test_paper_flower_fixture_v2.py`。检查 design/fixed 空间交并关系及不读取旧解，而不只检查字段。

**Exit:** 一个合法、冻结的论文式分层重建输入。切口位置和形状由新优化给出；若并非右下，不改解，报告输入/目标/坐标差异。

## Task 2：Heart整体准入；Noodles修复扫掠后冻结

**Files:** Heart复用 `experiments/faster_heart_contact_20261011/admit_bounded_request.py`、`launch_bounded_request_admission.py`、`bounded_pin_attachment.py` 及 `test_bounded_request_admission.py`；Noodles新增 `experiments/quadric_physics_demo_20261011/noodles_dense_v5_input.py` 和对应测试，复用 `diagnose_noodles_tube.py`、`noodles_dense_v3_admission.py`。

**Consumes:** Heart权威间距结果和有界pin附件；Noodles tube003自交见证及现有79管布局。

**Produces:** 两套整体输入准入、四视图、固定/设计域分区及冻结请求。

- [ ] Heart先运行有界pin与整请求反向测试：旧默认严格路径不变，仅显式策略采用经过认证的派生pin；错误摘要和缺附件必须拒绝。
- [ ] 在新独立作业核验整输入连接、pin与固定件对应、0.375mm间距、包含和完整分区。局部PIN_CORRESPONDENCE_CERTIFIED不升级为整输入通过。
- [ ] Noodles保留密集横向细管、内部连接、清晰握柄和底尖，拒绝五根粗管造型作为展示交付。保留旧tube003负例，修复中心线局部曲率与扫掠折返，不仅增加三角形数量。
- [ ] 以原布局、管径和未优化参考形态登记一次确定性曲线修订；逐管检查，再做真实并集、内部椭球与固定管束扣除、实体连接及分区。
- [ ] 两例各输出正面、侧面、透视和剖面，显示实际参与求解的全部初始材料；隐藏部分标重建假设。造型确认后冻结，不用求解结果反向改形。
- [ ] 容量界证明不可行则保留证书；未排除仅UNKNOWN，不阻断合规求解但不宣称已证明可达。

**Tests:** 运行 `test_bounded_pin_attachment.py`、`test_bounded_request_admission.py`、`test_noodles_tube_diagnostic.py` 及新增v5单管/装配测试；逐项保存真实native准入结果。

**Exit:** 两个合法可求解输入，或明确见证；渲染完成不计为优化结果。

## Task 3：三例的真实优化和独立材料门

**Files:** 复用现有运行器、精度策略、事件日志与材料接口；新配置/证据进入各自版本目录。不建立演示专用求解器。

**Consumes:** Task1/2冻结输入及适用生产导数证据。

**Produces:** 各例成功记录、独立终态检查和同版本重复。

- [ ] 花盘保持二维应用目标；心形和面条显式绑定Faster惯量—高度目标。请求包含单位、权重、密度、后端、离散与精度策略、实际加载二进制。
- [ ] 只补未覆盖生产映射的导数；三维沿用冻结81唯一点/三尺度协议，不能用Cuboid或诊断探针替代。修改有效映射才使证据失效。
- [ ] 冻结构建后按花盘、Heart、Noodles顺序串行运行；初始化不读答案，全部接受/拒绝/失败/资源终止及时落盘。
- [ ] 重载保存材料独立计算质量、质心、完整惯量、原目标与约束；检查固定材料、保留和互补空腔分区、连接与包含。
- [ ] 每例至少一个完整成功由独立进程从相同初态重复；两个预登记种子的全部记录保留。
- [ ] 失败只针对首个确证原因修复；输入修订需单独新版本和理由，不偷偷换权重/姿态/设计域。完整作业超限不通过拼接续跑授予成功。

**Exit:** 三例逐一满足真实优化、独立材料及重复；任何一例未过，总交付未完成。

## Task 4：分离可装配OBJ与诊断对象

**Files:** 修改 `paper_flower_output.py`、`paper_flower_saved_check.py` 并保留v1兼容；扩展 `test_paper_flower_output.py`、`test_paper_flower_saved_check.py`；三例复用同一交付清单约定。

**Consumes:** Task3同源材料，不读取预览模型替代。

**Produces:** 每例 final-assembly、before、diagnostics 三个清楚隔离的目录及可下载包。

- [ ] 增加反向测试：final-assembly只能包含最终实体和装配件；cavity/design/original出现即拒绝。转子、固定夹持组件分别列出。
- [ ] 统一边界构造/浮点转换，处理新版空中间固定网；保留所有旧失败OBJ不覆盖。
- [ ] OBJ重载验证自身拓扑、包含、无正体积重叠与完整分区，按既定几何/质量误差认证与权威实体的关系。
- [ ] 提供中间层单独OBJ及完整装配入口；另用诊断对象显示被移除材料，默认隐藏，不能装配填回空腔。
- [ ] 在中间层正面诊断图上标出实际移除区域；用OBJ几何而非重画理想圆证明缺口。摘要绑定图、网格及物性JSON。

**Exit:** 用户只导入final-assembly即可看到真实最终结构；不是要求其手动从15个混合文件猜出哪些要隐藏。

## Task 5：论文式表达与同对象真实物理对比

**Files:** 复用 `flower_paper_display.py`、`reference_studio.py`、`reference_studio_blender.py`、`flower_run_physics.py`、`flower_render_physics.py` 及现有Heart接触仿真入口。

**Consumes:** Task4合格OBJ、Task3物性和真实轨迹/事件。

**Produces:** 三例各一组输入/优化剖视/完整装配图、二维说明、交互三维、OBJ、物性JSON、正常/慢放视频及全场景统计。

- [ ] 固定灰、设计粉、保留黄、空腔透明；心形和面条必须显示真实握柄与底尖。展开/剖切只改变展示变换，不改变物理对象。
- [ ] 花盘新版必须重算理想轴承与柔性支撑矩阵，不能沿用v1的64场景。10Hz/.1支撑由新版原始移动质量确定，优化前后共用；不将轴承模型称桌面自由旋转。
- [ ] 心形/面条先通过无损刚体和完整惯量坐标变换基准，再执行35Hz、8秒、八方向2°、两时间步和两代理，各64场景。无世界固定支点/持续驱动，摩擦.8及既定零阻尼协议不改。
- [ ] 视频用首个登记相位/方向，并附全部统计；真实动力学轨迹而非相机动画。比较结论随代理/步长反转须标阻塞，不调参数制造稳定。
- [ ] 最终逐例核验优化、OBJ、惯量、图和物理摘要一致；独立报告数值成功、表示合法、物理比较结论及未标定制造状态。

## Scheduling and completion

输入与展示代码可以独立推进；正式实验依赖准入及适用导数，大型任务串行。优先新花盘形成第一个完整闭环，同时处理Heart准入、Noodles扫掠。不再安排仅用于美化未优化候选的新宣传渲染。

条件式工程估计：输入修订/准入2–4工作日；适用导数、求解、材料和重复3–6工作日；导出、物理与论文式交付2–3工作日。部分代码工作重叠，整体暂按7–12工作日预算，并非成功保证。前48小时以新花盘初次求解、Heart整体准入和Noodles单管/并集成本更新估计；遇到新算法或不可达阻塞单独报告。

结束条件是三例全部优化和交付门通过；不以到期、相似截图或两个成功案例替代第三个。2014五目标功能线继续独立记账，不被本计划替代或借用成功。
