# 三个真实优化结果 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Checkboxes record actual completion, not intended success.

**Goal:** 花盘、带柄心形、面条各有真实 Faster 优化终态、独立材料检查、可用 OBJ、重复证据及绑定同一对象的展示和物理对比。

**Architecture:** 花盘复用现有二维圆边界/dogleg及已合格三维挤出结果；心形与面条采用现有 Faster 三维二次曲面生产求解器。输入准入、导数适用证据、终态材料和动力学分别设门，不用渲染或历史其他案例证据替代。

**Tech Stack:** 现有 Python/native 求解器、CGAL、GMP/MPFR、有理材料、独立积分、PyBullet、Blender。

**Spec:** 当前聊天已批准的《三案例收尾规划》及三张参考图的展示修订；最新明确任务是取得三个优化后结果。参考图约束显示与重建目的，不提供优化答案。

## Global Constraints

- 所有案例为 RECONSTRUCTED_METHOD_VALIDATION，不称作者原始输入复现。
- 只保留固定外形功能；不恢复体素、形变、全论文统计、作者资产搜寻。
- 花盘使用其二维原目标；心形、面条使用 Faster 惯量比加质心高度目标，不混入2014平方高度目标。
- 三维显式 CGAL、level 2、冻结单位/密度/权重、种子0和1；每种子最多8次尝试，每次60次求值；alpha=t_max=0.5，raw residual≤1e-12。
- 解析矩≤1e-10；支持方向导数三尺度≤1e-6；独立最终矩及约束≤1e-6。
- OBJ内部边界偏差≤1e-6L，对称差体积≤1e-6L³，L为原输入bbox对角线；自身分区、包含、拓扑合法。
- 两输入线可独立开发；大型科学作业串行、24GiB/swap0/60分钟；编译≤4，数值库单线程。原磁盘准入要求不降低。
- 原结果、失败、所有修改保留；不提交、合并、删除或覆盖旧版本。

## Review Focus

1. 有效pin附件不能自动授予完整Heart输入通过，必须重算整个请求关系。
2. Noodles造型渲染不能掩盖管身自交；修补仅增加采样不等于修复扫掠折返。
3. 不得以既有Cuboid导数通过转授Heart/Noodles，也不得在检查清单冻结后换native模块。
4. 静态或初始化满足约束不能冒充优化；必须有实际接受、目标及材料变化。
5. OBJ、图像、惯量和物理轨迹不得来自不同对象；非对角惯量不能丢失。

## Task 1：花盘直接交付已有真实优化结果

**Files:** `experiments/quadric_physics_demo_20261011/paper-flower-method-v1/`，`reference-studio-v2/`；原结果不修改。

**Consumes:** solve-job-01/02、prepare-job-03、saved-independent-job-03、physics-job-02及qualified-flower-obj-v1.zip。

**Produces:** 已有真实解的来源索引与显示/物理下载入口。不是新优化次数。

- [ ] 重核两次求解、独立材料、64场景及所有导出SHA；只重新审查，不安排第三次相同求解。
- [ ] 明确二维圆边界/dogleg而非通用三维Newton；解释右下切除来自实际解。
- [ ] 区分立体静态渲染、真实物理视频和展示性展开。
- [ ] 运行已有来源/物理准入反向测试：

```bash
/home/lin/anaconda3/bin/python3 -m unittest discover -s experiments/quadric_physics_demo_20261011 -p test_flower_result_figure.py
/home/lin/anaconda3/bin/python3 -m unittest discover -s experiments/quadric_physics_demo_20261011 -p test_flower_physics_admission.py
```

**Exit:** 花盘标为已有优化、OBJ及工程轴承物理交付；不等待其他两例。

## Task 2：Heart完成派生输入整体准入

**Files:** `experiments/faster_heart_contact_20261011/admit_bounded_request.py`、`launch_bounded_request_admission.py`、`bounded_pin_attachment.py`及已有测试。

**Consumes:** authority admission job01、bounded request job02、PIN_CORRESPONDENCE_CERTIFIED companion job01；后者明确geometry_admission_granted=false。

**Produces:** 新独立整体准入结果及显式pin-policy身份，不更改外形或权威pin。

- [ ] 运行附件与整体策略反向测试，确保失配来源、缺附件、原策略误用和修改文件被拒绝。
- [ ] 原严格pin包含默认不变；仅显式策略采用与同一固定实体相交得到的有界派生pin。
- [ ] 在新job03重算附件与实际pin交集、完整分区、设计域包含、轴杆排除和0.375mm间距；记录原pin严格关系失败，不改写旧记录。

```bash
/home/lin/anaconda3/bin/python3 -B experiments/faster_heart_contact_20261011/launch_bounded_request_admission.py \
 --admission experiments/faster_heart_contact_20261011/rounded-clearance-v3-admission-job-01/output \
 --candidate experiments/faster_heart_contact_20261011/rounded-clearance-v3-bounded-request-job-02/output \
 --source-binding experiments/faster_heart_contact_20261011/rounded-clearance-v3-admission-job-01/source-binding-v1.json \
 --pin experiments/faster_heart_contact_20261011/rounded-clearance-v3-job-01/output/pin.obj \
 --job experiments/faster_heart_contact_20261011/rounded-clearance-v3-bounded-admission-job-03 \
 --native build/spin2014-objectives-quadric-20261007/python/_make_it_stand.cpython-313-x86_64-linux-gnu.so \
 --expected-native-sha256 ae617f1cf909be67611d0de7f51a3f47e8eb52ba4e5a841ac2623a6c666aab30 \
 --pin-policy BOUNDED_PIN_INTERSECTION_V1 \
 --pin-companion experiments/faster_heart_contact_20261011/rounded-clearance-v3-bounded-pin-job-01/output/terminal.json
```

- [ ] native发生行为修改时另行登记新哈希与证据，不修改以上冻结作业。
- [ ] 保留可达性UNKNOWN；容量界未排除不是可行证明。

**Exit:** 整体输入通过才能接入求解；失败固定具体见证，不盲目重启。

## Task 3：Noodles只修几何折返，再冻结合法输入

**Files:** 新建`noodles_dense_v5_input.py`及相应测试/准入适配，复用`noodles_dense_v3_input.py`、`diagnose_noodles_tube.py`、`noodles_admission.py`；不改冻结v3/v4。

**Consumes:** tube003独立精确交点见证、实际中心线与半径；参考图未优化结构。

**Produces:** 新候选预览、冻结manifest、全管准入和实际Boolean分区；不是仅视觉合格。

- [ ] 固定tube003为负例，先运行`test_noodles_tube_diagnostic.py`；27 native面交与独立21对严格交叉均保留。
- [ ] 保留79根横向管的布局、管径、握柄、底尖及完整包络设计意图；只修改造成高曲率折返的中心线扫掠。预先登记一次确定性几何修订，不按求解成败扫描参数。
- [ ] 为修订中心线登记局部曲率/半径约束；局部曲率检查仅为生成条件，不能代替完整单管自交判定。
- [ ] 每根管分别检查封闭、方向、退化、自交；再进行真实实体并集、内椭球扣除固定管束、连接及完整分区检查。
- [ ] 前后完整输入预览都展示实际初始材料，粉色候选域不把灰色固定管束重新标为可优化。
- [ ] 若完整包络容量界证明不可达，保存证书并提出单独输入修订；不为了优化成功偷偷缩小内部域。

**Exit:** 合法、冻结的新输入或明确几何阻塞；不能以渲染漂亮授予准入。

## Task 4：Heart/Noodles真实优化与独立终态验收

**Files:** 复用现有请求、生产求值器、材料保存、独立检查及受控作业系统；新增案例配置仅进入独立证据目录。

**Consumes:** Task2/3整体准入、请求与固定几何积分、相应生产导数证据。

**Produces:** 两个案例各自的实际成功与同版本重复，或明确失败。

- [ ] 审查生产导数适用范围；未覆盖路径按已有冻结81点/三尺度协议建立自己的记录，不能借Cuboid证明。
- [ ] 冻结Faster目标、归一化/单位、权威系数策略、CGAL、level2及实际加载模块；整个作业映射固定。
- [ ] 从请求/种子生成初值，不读取成功终态。串行执行seed0、seed1，保留所有筛选、失败和资源截断。
- [ ] 实际有接受更新、材料变化，目标和约束独立重算；raw residual≤1e-12，不把未收敛叫不可行。
- [ ] 保存固定、保留、互补空腔及系数；保存后独立积分、包含和完整分区，生成合格OBJ。派生修复不改变求解结果。
- [ ] 成功配置由独立进程从初始完整重复一次，比较科学事件和材料摘要；失败种子不替换。

**Exit:** 每案例至少一组完整成功及其独立重复；未通过者不显示黄色优化终态成功标签。

## Task 5：三案例绑定对象的最终交付

**Files:** 现有展示/OBJ/仿真入口；`reference_studio.py`、`reference_studio_blender.py`只读真实终态，保持原输入预览。

**Consumes:** Task1真实花盘结果、Task4合格Heart/Noodles终态及物理属性。

**Produces:** 每个案例的输入/优化/空腔/装配、OBJ与JSON、真实物理视频及场景统计。

- [ ] 三维剖视显示黄色保留、透明空腔、灰固定；不得画预期双曲面或强制右下切除。
- [ ] 心形/面条优化前后同外形、同碰撞代理，使用各自完整惯量；正确惯性坐标系变换。
- [ ] 无损刚体基准先通过，再运行35Hz、8秒、八方向2°、两时间步及两代理，每案例64场景；无持续驱动/固定世界支点，摩擦0.8、其余既定阻尼摩擦0。
- [ ] 花盘复用已有理想轴承/柔性支撑；不得描述成桌面自由接触或已标定人手。
- [ ] 视频采用首个登记方向/相位，附所有场景；结论随步长/代理反转则保留阻塞。
- [ ] 审计优化、OBJ、惯量与轨迹对象摘要一致，禁止无效OBJ进入接触。

## Scheduling and completion

先交付已有花盘；Heart整体准入与Noodles几何修订可独立开发，大型科学任务串行。输入通过后安排生产导数与优化，不做新宣传渲染占主线。

条件式预算：Heart准入/适配1–2工作日；Noodles修订/准入2–4工作日；适用导数、两个新案例优化/保存/重复1–3工作日；最终物理与包1–2工作日。存在重叠，顺利目标5–10工作日；新算法或资源阻塞另列，不能承诺期限必成功。首个生产点与完整作业实测后更新预算。

完成判定：三个独立案例全部有实际优化、权威终态、合格OBJ、独立检查、重复及绑定的物理记录。渲染完成、代码存在、已有其他功能成功均不代替该判定。五目标A线单独记账。

执行者不得自动提交、删除历史、延长60分钟或拼接被截断优化来声称一次完整成功。
