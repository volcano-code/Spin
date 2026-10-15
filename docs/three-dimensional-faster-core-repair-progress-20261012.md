# 三维Faster核心修正执行记录 — 2026-10-12

批准范围：同输入Cherchi/CGAL对照、实际求值计时、几何与Newton修复、性能收尾。论文§5.1/5.2与补充材料为数学基准，现有成功问题和保存失败状态为回归输入。

Ruling: 直接在用户批准目录新增隔离候选与证据目录，保留全部未提交修改、历史与旧模块；不提交、合并或删除。执行技能的提交、清理流程不适用用户约束。大型实验串行，24GiB/swap0/60分钟，编译≤4，数值单线程。

任务顺序：冻结诊断集合和分段计时 → 同操作数Boolean与导数对照 → 去重与几何修复 → Eq14独立复核 → 核心Cherchi接口 → 完整求解及性能验收。

初始事实：Heart core-v3 seed0在60分钟上限RESOURCE_LIMIT_WALL终止；未取得完整成功，保留耐久日志。Noodles第73次采集因ndarray JSON序列化失败，已保存文件是部分文件，不作为完整有理输出使用。

当前源码：value侧material audit、differential侧source/design/material audit存在重复调用；须对实际加载模块计时确认。适配器的拓扑哈希混入source_scientific_hash；标签不影响当前残差下降分支，不能据此断言它是收敛根因。

第一批使用隔离诊断代码，先生成实际固定几何积分的小型箱体请求，固定三种二次曲面状态；四次求值（冷启动一次，同进程三次）与Cherchi同操作数对照。小型状态仅作诊断，不作为模型优化成功。

第一批实测：same-input-job-03在26.577秒内完成。三状态各执行四次CGAL值求值、四次Cherchi同操作数Boolean；另一次明确计数的fullrank调用提供生产Jacobian。质量矩最大误差分别7.11e-15、3.55e-15、2.66e-15；Jacobian最大误差分别1.38e-14、6.39e-14、7.11e-15。尚未认证实体对称差，因此不是完整后端等同通过。

CGAL值求值中位数约0.057、1.284、1.210秒；Cherchi进程+Boolean+导出约0.0058、0.0323、0.0232秒，两者工作范围不同不能直接报告加速比。Cherchi首次来源运动与导数约0.132、0.238、0.161秒。双曲面CGAL总时间约1秒未被现有observer覆盖，增加canonical导数路径阶段计时。

计时候选01独立构建/加载/三状态检查完成，生产fullrank完整输出与旧模块逐JSON相等，但旧canonical导数对象仍未传内部observer。候选02编译链接成功，加载失败：原二进制的canonical导数对象缺少新增observer重载；保留失败记录，未进入生产。候选03补独立canonical导数对象，所有函数别名单独命名，正在构建；成功加载后重新比较数值及阶段。

额外诊断入口失败记录：same-input-job-01输入faces为NumPy整型导致Python入口拒绝，改为标准整数列表；job02来源设置类名错误，纠正为实际SpinTriangulationSettings。两次都未授予科学通过。

Noodles采集的ndarray序列化已修正，保留旧部分输出。大型构建期间不启动采集重放；在空位进行一次明确修复目的的重跑。

计时候选03已独立编译、链接、实际加载并完成三类曲面对照；verify_timing_candidate.py确认三类完整生产输出逐JSON相等，canonical源点运动/材料transport阶段均完整，每类一次。单/双叶源运动重建约1.13秒；材料边界在value和transport各一次审计，source/design也单独审计。候选未被生产采用；实体等同、差分、完整Newton及性能门仍未关闭。

已启动heart-state-timing-job-01（exec session55203），对core-v3第一尝试最后实际更新的保存状态执行一次实际生产映射。只为复杂状态分项计时，不是重新优化或续跑；运行期间不修改计时候选和绑定依赖。

Heart单状态诊断返回真实映射，72.691秒：Boolean0.622秒，源面0.780秒，canonical源运动0.819秒，值侧材料审计3.638秒，导数中固定design审计41.093秒，重复material审计3.687秒。其余未细分的约21.8秒仍需继续调查，不能转述为Boolean耗时。

第一项性能修复：仅缓存成功的固定design正向闭合边界证明，thread_local单条容量，完整hex binary64坐标、面索引、形状和精确谓词策略逐字绑定，不缓存移动source/material。新增真实native回归先观察第二次无cache-hit而失败；同一输入应复用、顶点/策略改变必须失效。候选fixed-design-cache-build-job-01正在独立构建（exec3012）。源码修改只影响候选，旧运行模块未覆盖。

第一项修复已独立加载、3项相关回归及三类生产输出一致检查通过。Heart同进程同状态冷调用72.755秒，随后32.120秒；固定design审计41.180秒→0.00663秒。这是约2.27倍单状态暖调用提速，不代表完整优化或论文性能通过。

第二项修复定位computeSpinDesignRequestHash也再次调用validateSpinDesignRequest，每次outer/design固定输入审计重复。候选仅在validateMesh内部缓存成功的完整索引几何和orientation/exact策略，容量两条，所有请求标量、轴、密度和目标检查仍执行。fixed-input-cache-build-job-01正在构建（exec21209）；需要真实暖输入无效几何/密度/轴回归、三类科学输出一致及Heart实测后才能采用。

第二项已完成独立构建、加载、4项相关回归及三类完整生产输出一致检查。Heart冷调用58.914秒，随后11.491秒；两次科学求值与旧模块完全相同。当前Heart outer4788面、design45208面，不能按作者10k面输入直接比较论文秒数。

第三项修复：auditQuadricMaterialBoundary以完整索引几何和精确策略绑定两条成功边界证明，供同次value与transport复用，未缓存Jacobian。boundary-proof-reuse-build-job-01独立加载并通过上述回归与三类输出逐JSON一致。Heart同状态冷55.038秒，重复4.148秒；同状态重复包含边界缓存，因此不当作普通新Newton步成本。

预登记新参数检查：只将保存状态的constant加0.001并归一化，作为诊断试探而非优化初值。候选heart-distinct-state-cache-timing-job-01第二次实际构造返回，8.005秒，新的材料value边界检查3.790秒仍实际执行，transport复用0.00287秒。旧模块在heart-distinct-state-original-baseline-job-01以完全相同系数/输入重放：72.608秒；两份完整生产evaluation逐JSON相等，系数逐值相等。该固定试探单次成本约9.070倍改善，冷输入准备仍约55秒，尚非完整优化提速或性能门通过。

本批实际交付：分段计时候选、三个同操作数Boolean/Jacobian组件对照、保存失败集合、固定输入/设计域与同次边界证明复用、真实无效输入/策略失效回归、同一Heart试探的前后科学一致与计时。

仍未关闭：后端精确实体等同；Cherchi进程内生产接口；双曲面封口及合法细小面/真退化区分；Newton装配独立复核与适用三步长导数；小型完整求解/独立重复/复杂度分档性能。优先继续这些，不无修改重跑Heart或恢复旧81点批次。

## 2026-10-13：进程内后端、小型完整求解与独立复核

上述未关闭列表为当时状态，以下为新增证据，不覆盖旧记录。

- `exact-entities-job-01/output/report.json`：椭球、单叶、双叶同操作数CGAL/Cherchi有理实体均完成精确集合对照，正向闭合且对称差体积为0/1。限于这三个冻结状态。
- `cherchi-inprocess-build-job-01`：固定版本Cherchi实际Boolean及构造来源已在独立pybind模块内运行，不逐次启动程序或写读操作数。`cherchi_core_map.py`使用同次材料计算质量矩、来源运动及Jacobian；当前是数值诊断适配器，尚未获得正式封存/通用生产准入。
- `small-top-inprocess-solver-sequence-01`：真实非零固定壳、下方支撑原点的Faster惯量＋线性质心高度新方法夹具，答案无关seed0/1。seed0为87.688秒、133次求值、尝试[60,60,13]、raw残差1.908e-15；seed1为19.996秒、23次求值、raw残差3.678e-16。全部失败、拒绝与重启保留。不是论文Heart输入，也不是复杂模型完成。
- `small-top-inprocess-solver-repeat-01`：两个种子各独立进程完整重复，约88.644/20.048秒。`small-top-repeat-comparison-01.json`确认配置、所有非耗时科学事件、完整数值结果和终态材料逐值一致。原event_hash包含耗时而不同，未篡改；比较仅剔除明列的三个耗时字段。
- `small-top-derivative-probe-job-02/output/report.json`：规则双叶终态十方向、三个预登记尺度、81唯一生产求值点全部完成，逐项最大误差2.090e-9。保存材料以独立GMP精确四面体公式重新积分，与生产矩最大差1.776e-15。源面分支保持稳定，但材料连接有13种身份；未完成同实体共面重三角化分类，`branch_certification=false`，不据此授予完整导数门通过。
- `independent-eq14-pilot-01.json`：独立SymPy公式重算两种子有payload的63/23次求值，目标、乘子、raw残差最大差5.330e-15；首个接受Newton球面更新差≤8.882e-16。离线重建矩阵与方向的作用得到验证，但成功事件未保存native矩阵，尚非逐元素原生矩阵比较。
- YOYO小型对称壳夹具曾数值收敛14.149/24.486秒、11/29次求值，但终态趋近a=b=0的圆锥退化极限（符合论文§4.5对称情形）。只保留为退化观察记录，不能代替普通满秩三维成功/导数证据。

新反向/手算检查：重复比较必须检出目标、失败事件与材料身份变化；独立Eq14各向同性手算目标和交叉项坐标约定；独立精确四面体平移解析盒基准。新增测试先观察缺实现失败再实现并通过。

当前剩余关键门：材料连接变化的精确分类；实际失败封口/细小面回归与修复；原生Newton矩阵捕获逐项核验；Cherchi核心正式生产接口及完整终态准入；万面/十万面分档计时及复杂输入收敛。约20秒的小型seed1已进入工程量级，但seed0约88秒含重启、Heart暖新试探8.005秒、冷输入约55秒，均不能授予整体性能通过。Heart/Noodles仍无新的完整交付，不重启旧作业或恢复旧审计队列。

## 2026-10-13 00:40后：原生矩阵准入、连接取证与论文源面缺口

本轮未发现上述大作业仍运行；未恢复旧作业，磁盘仍有约409GiB可用。

### 原生Newton矩阵捕获已补齐（限小型前缀）

新真实回归`test_successful_newton_capture.py`先在旧候选观察0份成功矩阵而失败。定位为`QuadricHollowingSolver.cpp`只在失败或行列均衡时调用已有captureNewtonSystem。最小修复仅在显式captureEvaluationPayload时也发出普通成功线性系统事件，不改变矩阵、方向、求值或回溯。

`newton-capture-build-job-01`独立编译链接并实际加载，不覆盖旧模块。真实3次求值测试证明没有新增材料求值，原生15×15矩阵按请求交付。8项相关native/适配器/缓存/几何回归通过；10项独立公式/比较/精确积分/计时/平面证书测试通过。当前不是完整里程碑，不跑全项目长测试，也不声称全项目通过。

`native-newton-matrix-audit-job-01`因诊断配置含dataclass不能JSON序列化而失败，未进入材料求值；改成asdict后job02保存实际3次求值、3份矩阵，但检查器错误预期只有2份而报告失败。代码事实是预算耗尽后的下一轮在回溯之前仍装配最后接受状态矩阵，额外矩阵没有增加一次材料求值。Ruling: 保留全部失败，修正诊断计数为3并只离线重审保存事件，不为报告错误重跑优化。

`native-newton-matrix-recheck-job-01/output/report.json`：0次新增映射求值，3份原生矩阵和右端项分别与独立SymPy Eq14重建逐项比较通过；矩阵差≤8.882e-16，右端差≤3.470e-16，乘子差≤6.939e-16。3次求值的系数、目标、raw残差、乘子、质量矩和Jacobian与旧小型seed1前缀逐值一致。只关闭这3个规则状态的装配/观察缺口，不代替复杂Heart或完整新版本重复。

### 81点连接变化尚未认证，不升级导数门

新增精确有理共面边证书（不是近似采样）：拒绝真正零面积、非闭合、错误绕向；1e-60高的非零面保持合法；对角线翻转只有精确共面才合并，1e-30非共面不能通过。真实`small-top-planar-patch-classification-job-01`只读取已有81份材料，不调用Boolean。

所有点76顶点/148面，但索引面有9种；精确共面合并后仍有61–65个patch，索引边界不全部一致，故`full_branch_certification=false`。初步空间匹配诊断发现基点→point002顶点17/18/19循环重编号（最大位移约5.70e-4）；空间近邻仅用于取证，不能当作来源对应证明，也不能断言发生真实实体拓扑改变。下一步需从构造来源建立跨状态顶点身份，区分重编号、共面重三角化和真实几何变化，不改变基础步长或重复整套81点。

### 新确认的论文源面支持阻塞

`saved-source-profile-comparison-job-01`预先冻结两份Heart失败系数，分别只调用旧法向构造和显式论文projective构造；共4次源面求值，无Boolean、无Newton、无输入修改。Heart冻结配置实际为`paper_projective_hyperboloid=false`，因此既有Heart轨迹不能称为论文projective源面路径通过。

- 原封口失败：旧路径0.003秒重现顶点8有2出边/0入边；论文路径0.038秒拒绝`UNSUPPORTED_PAPER_IDEAL_FACE: multiple crossings at face 0`。
- 原细小面失败：旧源面本身0.747秒返回，细小面拒绝发生在随后材料Boolean/审计；论文路径0.044秒拒绝同类多重交点face125。

## 2026-10-13：批准的新修复批次启动

执行用户《三维Faster修复与心形、面条交付计划》。任务顺序为固定失败取证/几何修复、正式生产链与适用验证、完整运行及交付。Ruling: 延续当前独立候选目录和证据账本，不建立新的提交/清理流程；用户明确要求保留脏工作树、不提交、不合并、不删除。当前没有模型或候选构建进程在运行，磁盘约409GiB可用。

首个串行作业是面条第73次Boolean完整有理采集（capture-job-02）；仅重放已修正序列化的失败采集，不恢复旧优化或81点队列。源码调查同步审查Heart论文projective多重理想交点支持。几何修复之前，两个模型仍无完整优化结果。

### 本批实际修复和验证结果

面条`capture-job-02`完整采集成功，原有理输出的4个面均为真实非零面积，转binary64后塌缩；3个坐标别名还产生两对反向重复面。新增`bounded_rational_core_mesh`先验证有界误差，再抵消反向面和合并相同binary64坐标；同向重复面、无边界见证、超误差或非流形顶点均拒绝。相关3项回归通过。第73步的边界距离L1保守上界2.767e-10模型单位，扫掠体积上界1.101e-9，远小于预登记误差阈值。该证书不代替独立实体有效性检查。

`noodles-dense-v5-core-preprocessing-v3/domain-job-01`以有界派生边界反馈继续第73、74步，随后第75步出现4面共用短边并终止；保留该实际有理输出。定向收缩该边的诊断虽闭合却产生不连通顶点link，故没有采用，也不删除部件。

Ruling: 第73步之后改用已有`libigl_exact_boolean_construction_from_rationals`传递实际有理操作数，只在最后导出binary64，不逐步反馈舍入网格；登记为`CGAL_CORE_RATIONAL_SUFFIX_V1`新输入构造profile，不升级历史binary64结果。保留相同72步前缀和全部81个固定部件、几何、坐标、密度和轴。`noodles-dense-v5-core-rational-suffix-v1/domain-job-01`实际完成全部81部件，749.848秒（约12.50分钟）为后缀输入准备时间，不是优化时间。第75步的非流形边在该有理路径中未重现；最终167170面，导出`design-core.obj`和`domain-081.json`。派生边界偏差上界2.767e-10、扫掠体积界1.276e-9；输入整体、独立边界及实际固定材料尚待准入，不能称为优化成功。

心形`normal-face-inspection-job-01`使用实际固定球面和inverse-Gauss根，确认两个保存状态各有4个多重交点面，包含单边两个理想根。新增有限正法向带裁剪，按正区间及负锥来源保留全部有限边界分量。候选01错误假设每个负锥只产生两个根，候选02修正为全部正边界路径，但逐弦面梯度定向仍破坏共享边；候选03改为由精确Gaussian曲率符号决定固定法向图整体方向，源面闭合且可微，但第一状态的均值扇形远端封口与interface有8组真实相交。扩大半径或把均值径向推到球面的诊断没有消除这8组相交，未采用。

Ruling: 对多重交点法向带使用显式远端平面延伸封口。由同一截断环Newell方向建立远端平面，原环和新环均严格位于设计AABB同一支持平面之外，精确投影ear裁剪闭合后执行整源面精确边界检查。新增顶点沿用已有ray来源槽编码`{-2, interface-ring-vertex, loop, ordinal}`，几何hash包含`remote-planar-cap-extrusion-v1`；未改变接口数据布局。`QuadricSourceMotion`同步传递环方向及投影构造的真实导数，不填零cap Jacobian。普通既有源面及旧接口继续保留。

`projective-belt-build-job-04`独立编译Triangulator和SourceMotion两个对象，链接保留已验证缓存/Newton捕获对象并实际加载。`test_saved_heart_projective_source.py`在旧模块重现两个UNSUPPORTED错误，在新模块通过两状态源面、来源、motion、selector及精确闭合/定向/流形/无自交检查。16项本目录相关回归实际通过；不是全项目测试通过。

`repaired-heart-fixed-states-job-01`进一步调用真正CGAL全映射：原封口失败状态返回5697顶点/11390面、有效Jacobian，56.213秒包含冷输入验证；原细小面失败返回11978顶点/23956面、有效Jacobian，10.658秒。两份材料、质量矩、Jacobian和分段耗时已保存。不是完整优化或最终收敛，配置显式`paper_projective_hyperboloid=true`，不沿用旧false profile成功记录。

`repaired-heart-source-derivative-job-01`预登记h=3e-6、h/2、h/4，两个保存状态各81个源面求值点，十方向五点差分通过；最大逐顶点相对误差2.618e-7和1.144e-8，interface/点来源/radius及cap连接稳定。该证据仅为source坐标导数，不代替生产质量矩Jacobian验收。

### 当前正在运行的串行队列

`launch_repaired_verification_sequence.py`启动`repaired-production-verification-sequence-01`。先执行`repaired-heart-production-moments-job-01`（原封口失败状态）的81点真实CGAL生产映射及保存后独立GMP积分；数值门通过后执行job02（原细小面失败状态）。两个作业分别24GiB/swap0/60分钟，十方向、三尺度均固定，不事后缩步。只有数值通过不能升级未完成的材料连接认证。

随后队列执行`noodles-dense-v5-core-design-audit-v1/audit-job-01`，一次精确binary64 indexed-boundary与独立质量矩检查，不加入逐点Nef/制造审计。设计域通过后仍需实际outer并集、fixed材料、完整分区和容量界；未齐备时禁止面条进入优化。

队列运行期间不修改其绑定脚本、模块和输入，不另起大型作业，不恢复Heart旧dispatcher、旧81点或旧Nef队列。下一任务是读取本队列结果，关闭适用材料分支认证及实际Noodles固定几何，再将进程内Cherchi正式材料入口接入并执行seed0/1完整优化；成功后独立重复、合法OBJ、MuJoCo+Blender同材料仿真。Heart/Noodles目前仍未有完整优化交付。

### 串行验证终态与新的定向修复

`repaired-heart-production-moments-job-01`实际81点完成，约760.34秒；十方向、三尺度全部数值通过，最大误差1.515e-8，保存后独立GMP积分最大差2.554e-15。只关闭该保存状态的数值Jacobian检查，不授予尚缺材料连接认证的完整导数门。

job02实际返回点0–3并保存，第4点因旧几何面积阈值而拒绝：真实浮点面squared double area=4.08584e-24，闭合、edge/vertex manifold、绕向与正体积均通过，但旧阈值记为degenerate并在执行精确相交检查前退出。该失败保留，不缩步，不续拼整套通过。`heart-exact-material-gate-diagnosis-job-01`对同一manifest第4点使用现有显式`exact_binary64_material_boundary=true`，实际材料及有效Jacobian返回；87.054秒，其中Nef材料证明25.337秒，冷design约41.967秒，native及canonical transport各一遍约3.6秒。证实应使用准确几何政策，而不是删除小面或放宽科学门槛。

新增显式`paper_core_settings`与Heart运行器`--geometry-profile paper-exact`选择`FASTER_PAPER_CORE_EXACT_BOUNDARY_V2`，旧默认入口仍保留。完整优化尚未启动，新的production profile需先针对失败点、相关回归及受影响证据验证。

面条独立设计域检查明确失败：`noodles-dense-v5-core-design-audit-v1/audit-job-01/output/boundary-audit.json`显示闭合且edge manifold，但顶点20166非流形，78组精确自交，12499604次精确pair tests。因此`design-core.obj`不是合格输入/仿真网格，`geometry_admission=false`保持。此前的有界距离和扫掠体积证书没有被用来覆盖这些失败。

当前实际运行：`noodles-rational-authority-diagnosis-job-01`独立CLI对原有理domain-081执行indexed-simplex精确边界诊断，无Nef、无转换和修复，用于区分72步binary64前缀继承缺陷与最终派生舍入缺陷。期间不重跑预处理、优化或另启大型构建。

待上述作业终态后，已启动的`core-material-followup-sequence-01`按串行规则构建`exact-core-material-build-job-01`，并重测第4失败点。该候选对严格精确嵌入、闭合、可定向、只有一个正体积边界分量的材料，使用Jordan-Brouwer区域证明复用indexed-simplex事实，避免再次构造Nef；多个分量保留原有signed-shell/Nef嵌套验证。改变只在显式exact material政策生效，身份标记v2，几何、面积非零、流形和精确自交门不放宽。候选未验证，不能登记为性能成功。单作业24GiB/swap0/60分钟，源码和旧构建全部保留，不提交、不合并、不删除。

相关结果：核心目录16项回归、Boolean表示6项回归实际通过；这些是局部检查，不是全项目或两个案例完成。下一步必须从有理诊断结果决定完整输入构造修复或有界派生重三角化，不能让失败OBJ进入求解/动力学。

### 2026-10-13：有理边界终态、连接见证与精确材料候选复核

上述两个作业与followup均已终态，没有重复启动。`noodles-rational-authority-diagnosis-job-01`原有理边界检查约886.43秒：闭合、edge manifold、无零面积面、无精确自交，但vertex20170非流形。局部`vertex-link-witness.json`证明link由两个独立环组成（8和3个邻点），不能靠顶点编号拆分就声称合法嵌入。浮点派生物另外78组相交仍是独立出口缺口。

`prefix-link-trace.json`只读历史前缀，发现相同binary64坐标在domain-072顶点17665已经有两个link环（8和5个邻点）；此前匹配到的70、71前缀没有该坏star。这是局部继承缺陷证据，不是整个71前缀已认证。下一定向取证应重放第72步的原有理Boolean并检查此star，以区分该步真实点接触与binary64反馈破坏；不盲目重跑全部81部件，不删除部件、不收缩短边。整体outer/fixed/容量界仍未准入，面条无优化结果。

`exact-core-material-build-job-01`已实际编译链接并加载。`heart-exact-material-gate-diagnosis-job-02`同一个失败point4返回真实材料及有效Jacobian，29.715秒。与job01保存evaluation逐字段比较：材料坐标、面、几何身份、质量矩、Jacobian及判定完全相等；仅scientific_hash因显式几何政策v2身份不同。该单状态墙钟对照不能冒称完整优化或总体性能达标。新候选核心目录16项回归实际通过8.214秒，顶点link取证3项测试通过。

当前串行大型作业为`repaired-heart-exact-production-moments-job-01`：明确选择exact material政策，从同一保存细小面失败状态重新登记完整81唯一点，十方向、三尺度和原基础步长不变。它是确证面积误判修复后的新profile验证，不拼接旧job02点0–3、不缩步；全部材料保存后独立积分。运行期间不修改绑定依赖。数值通过仍需适用材料连接认证，不能立即转授完整优化准入。Heart与面条均仍无完整优化交付。

### 2026-10-13 04点后：第二份Heart数值门完成，面条舍入反馈根因对照

`repaired-heart-exact-production-moments-job-01`实际完成81点，约1486.31秒（24.77分钟），最大三尺度五点误差7.894e-8、保存后独立积分最大差1.199e-14，均达标。`full_material_branch_certification=false`与`optimization_admission=false`保持，不能把数值门升级为完整优化或连接认证。

`noodles-original-boolean72-capture-job-01`对原71浮点前缀实际重放tube-069-fixed，原有理输出同坐标vertex17665仍为两个link环。随后预登记的`noodles-rational71-72-diagnosis-job-01`只改变第71与72之间的表示传递：从相同保存70前缀构造71原有理输出，直接送给72，无浮点反馈。原坏坐标精确不变，71/72两个star均成为单个闭合环。31.78秒为该定向两步诊断，不是完整准备/优化。

Ruling: 此对照支持“第71步舍入边界反馈破坏第72步局部star”根因，采用`CGAL_CORE_RATIONAL_FROM71_V1`继续相同81部件，仅提前保留有理操作数；不改几何、管束、轴、密度或目标。局部修复不代替整体准入；仍可能有其他前缀或最终派生缺陷。

当前串行大型作业`noodles-repaired-rational-suffix-job-01`从新实际72有理输出继续73–81并逐步封存。等待队列`noodles-repaired-input-sequence-01`在其终态且完整输出存在后，执行`noodles-repaired-rational-audit-job-01`，复用已绑定精确有理indexed-boundary CLI，不再次编译，不同时启动大型任务。队列不运行优化；通过完整权威边界之后仍需派生出口、outer/fixed实际积分、分区和容量界。Heart与面条目前均无新的完整优化交付。

### 2026-10-13 04:43后：新面条后缀完成、核心离散身份修正

新面条有理后缀已实际完成全部81固定部件，940.758秒为新后缀准备时间。`noodles-repaired-rational-audit-job-01`正在完整独立精确边界检查，运行中的输入/CLI/worker未修改。新增`noodles-repaired-export-sequence-01`只等待该审计终态；仅权威边界通过后才执行新有界派生表示和独立binary64边界检查，失败则不启动。输出名字`design-candidate.obj`且保持geometry/optimization admission false，不把未检查的文件当作可用OBJ交付。

核心接口复核确认诊断Cherchi适配器旧topology hash含连续parameter_t及输出顶点编号。新增显式`SOURCE_FEATURE_BRANCH_V2`，按同一次构造的唯一input-object/input-vertex特征代替输出编号，保持有向三角形，剔除连续插值参数；歧义或重复特征明确拒绝。source连接/类型/封口信息仍参与身份。历史默认legacy原topology算法保留；新策略进入上下文身份。不自动转授共面重三角化或完整邻域可微证明。

新增4个身份行为测试先观察缺实现失败，随后实现并通过；另有歧义/重复特征拒绝测试。真实椭球、单双叶三状态显式新策略与legacy的实际质量矩、Jacobian及材料几何逐值一致，身份不同。核心目录27项局部回归实际通过15.186秒，未执行全项目里程碑测试，不宣称全项目绿或案例完成。

### 2026-10-13 05:13后：面条新整体失败与输入构造架构纠正

`noodles-repaired-rational-audit-job-01`完整精确有理边界检查1065.986秒后仍失败：closed/edge manifold/无精确自交/无退化面均成立，但vertex20149有两个link环（3和4个邻点）。它与此前vertex20170不是同一精确坐标，距离邻点多为1e-11–1e-10模型单位。只读edge连接取证表明两个star属于同一个167134面组件，并不是可以独立移除的小碎片。`noodles-repaired-export-sequence-01`明确停止，没有生成新合格OBJ、启动优化或仿真。

Ruling: 局部抵消别名、保留73以后有理、提前到71的三条尝试均未得到整体合法输入，不继续修补带binary64历史前缀。构造架构应在第一个Boolean起保留同一有理边界；本轮登记`CGAL_CORE_RATIONAL_ALL_INPUTS_V1`一次原始输入对照，所有81固定部件和原outer-ellipsoid不变。原始binary64坐标精确编码成numerator/positive denominator，未改模型、轴、密度、权重或删部件。该对照仍可能揭示原输入真实接触，不能保证成功；通过全域检查之前禁止优化。

`noodles-all-original-rational-job-01`在首个Boolean参数解析前失败：GMP对零/整数的str形式缺少/1，native要求显式分子/正分母。保留该失败，未发生几何运算。新增规范binary64坐标编码器以as_integer_ratio保留原值；零/整数、非十进制改写、NaN/Inf三项测试先RED后GREEN，旧表示6项回归仍通过。`job-02`已使用确证修复从原始输入开始执行全链构造。单作业24GiB/swap0/60分钟，不续拼超时作业。

当前大型作业`noodles-all-original-rational-job-02`；串行等待`noodles-all-original-input-sequence-01`仅在81部件全部实际返回后启动`noodles-all-original-rational-audit-job-01`，复用现有CLI。旧失败输出/退出队列全部保留。Heart已完成的两份保存状态数值门不重复；适用材料连接证据与正式Cherchi准入仍未齐备，两个案例仍未优化完成。

### 2026-10-13 05:43后：运行中的原始构造与零Boolean来源试点

全原始有理构造仍运行，约29.71分钟完成39/81部件，正在40；没有启动其他大型构建/实验或改变该作业绑定依赖。它比有舍入前缀的构造慢，仍由原60分钟上限控制，不因进展慢增加预算或续拼成功。

利用已有small-top点000/002保存的有理材料，用同一原冻结native模块重新生成2份源面（0次材料Boolean），实施精确三角形内包含、精确共面patch与精确共线边界细分见证。5个真实算术测试先观察缺实现失败后通过；1e-60非共面/非共线偏差不会被近似吞掉。仅用于证书视图，没有删改实际材料顶点和面。

`saved-support-incidence-pilot-01`支持平面身份不足：每点有8个rank2细分点且存在重复键。`pilot-02`以精确共线见证在证书视图消去这些冗余细分，得到每点68个唯一rank3核心角点，但跨状态精确patch来源标签仍不匹配，`matched_patch_boundaries=false`、`full_branch_certification=false`保持。不能用近邻对应、误差小或测试数覆盖该缺口；下一取证应保留同次Boolean的原始顶点类型/操作数来源，区分原输入顶点与交点，再处理偶然共面标签变化，而不是重复整套81点。

该取证暴露出Cherchi诊断适配器在返回后丢弃原始construction records。新真实测试重现last中0记录、实际170顶点而失败；最小修复在同次调用的last附件保留原始精确操作数记录和已解析来源，不额外调用Boolean、改坐标、积分或Jacobian。历史默认topology算法不变，附件为`same-call-cherchi-construction-v1`；不能靠附件存在授予完整邻域或封存准入。修复前核心32项局部回归通过15.612秒；附件修复后的全部本目录回归正在执行。

附件修复后核心目录33项局部回归实际通过16.672秒，包括真实同次170顶点逐索引construction记录、三状态科学输出及相关失败回归；未执行全项目测试，不授予正式核心或案例完成。全原始面条构造约35分钟完成41/81，正在42，串行审计仍等待；作业参数与上限未修改。

### 2026-10-13 06:13后：全有理构造资源终止、必要来源取证

`noodles-all-original-rational-job-02`已`RESOURCE_LIMIT_WALL`终止（3600.11秒）：实际完成50/81部件，正在51时截断。MemoryPeak=1010872320 bytes，swap peak0；不是内存耗尽。`noodles-all-original-input-sequence-01`明确记录不完整构造、不启动后续审计。所有完整前缀/失败和输入身份保留，不增加预算、不续拼成功、不无修改重启该整批。

Ruling: 原始输入全程有理构造尚不能在原预算内交付，下一步应对保存第51次操作做预登记分段性能定位，并审查不损精度的构造组织方式；不能再以相同顺序一小时长跑替代修复。输入并集/固定材料/容量界及优化均未完成。

`saved-small-provenance-capture-01`只补两个必要取证点000/002，共2次小型材料Boolean，未重跑81点：与原保存的有理有向三角集合精确相同，各76顶点均保存原始类型和精确操作数来源；构建、适配器、bridge及保存材料摘要记录。不是导数整批重跑/新优化。

新增明确的约束身份正规化：LPI的线两端与平面三顶点置换不改变定义约束；TPI三平面及其角点置换同理。3项真实行为测试先RED后GREEN。仅供来源证书分析，旧legacy与SOURCE_FEATURE_BRANCH_V2求解身份未改，不重新排列实际导数操作数。两取证点得到完全相同的76个唯一约束键集合，从而闭合这对状态的顶点重编号对应缺口。

`saved-small-source-patch-classification-01/02`仍未通过完整patch比较，失败保留。进一步`saved-small-fixed-geometry-retriangulation-01`将两份连接模式分别应用于各自同一固定有理几何：point002精确共面边界相等；point000的另一模式产生精确零面积三角形而拒绝。该零面积只在证书用的假设重三角化中出现，未修改原实际合法材料；不能把两份不同参数的实体直接比较为相等，也不能静默删除证书退化面授予通过。完整81点分支认证仍false；下一步需针对这个固定证明状态做合法有向平面链证书，而非重复数值差分。

核心目录36项局部回归实际通过17.186秒；不是全项目/完整核心/两个案例完成。当前原始输入构造与等待队列均终止，先核实无残留重作业，再执行有明确目的的性能/来源证明工作。Heart旧优化及旧dispatcher仍不恢复。

### 2026-10-13 06:43后：性能分段、81保存点来源链认证及输入更正

`noodles-rational51-phase-diagnosis-job-01`编译完成但链接缺Boost TLS符号而失败，未运行Boolean。job02验证原探针源码摘要后复用已编译对象，补齐实际native构建已有的Boost静态链接输入；旧模块/对象不覆盖。实际固定操作的分段为parse-operands0.181680秒、mesh-boolean128.465286秒、serialize-result0.064秒，正式native同操作133.911150秒，两者有理有向三角集合精确相同。没有继续第52步或恢复整批构造。

**必须纠正输入范围：**本轮核对manifest与旧core入口发现，我新增的`construct_noodles_original_rational.py` job01/02误将outer-ellipsoid作为设计域种子，冻结公式实际为`DIFFERENCE(inner-ellipsoid-candidate, UNION(ALL fixed_parts))`。因此上述全链50步资源结果和第51次分段只能作为误配外域的诊断，不是冻结面条设计域准入/性能证据，不能用于正式案例预算或成功声明。原始记录不修改、不删除；前面从正确domain-070/071/072前缀取得的失败见证不受此错误影响。始终无优化或制造/仿真通过声明。

已明确向用户更正该错误。新增manifest-owned设计域选择器，3项测试先RED后GREEN，检测内域与外表面偷换、缺少种子和未知公式；旧表示6项、坐标编码3项回归仍通过。全链worker改为`CGAL_CORE_RATIONAL_ALL_INPUTS_V2_MANIFEST_DOMAIN`，显式保存种子名/摘要并绑定inner输入；旧顺序整批launcher暂时拒绝启动，性能构造修复前不安排又一小时重跑。下一构造应先在小型几何上核验原manifest的集合语义，并评估保持精确坐标的平衡固定件并集再单次内域切除，不改变81部件或冻结几何。真实fixed仍是outer减design，不能把固定件并集直接当作全部固定材料而漏掉外壳。

针对small-top证明视图的零面积假设面，新增精确有向平面2-chain证书：按每个精确支持平面的有向边界1-chain进行线段端点分割、整数重数求和与共线合并；精确零2-chain仅在证书视图记录，不删除任何实际材料。4个测试先RED后GREEN，覆盖对角线翻转、证书空三角形、反向面与1e-60非共面差异。`saved-small-fixed-geometry-retriangulation-02`在point000和002各自固定有理几何上证明两模式代表相同有向表面链；这不是无效假设网格变成可制造网格。

随后`saved-small-recipe-chain-certification-job-01`只读已有81份有理材料，使用同一冻结native重放81个源面（**0次材料Boolean，未重跑数值差分**）。以基点同次EXPLICIT3D/LPI来源独立重建每个有理顶点，81点均建立精确双射；每个参数下自身连接与基点连接均给出相同有向平面链。`complete_frozen_sample_source_recipe_chain_certification=true`。独立LPI重放的手算、奇异拒绝2项测试先RED后GREEN；本目录42项局部回归实际通过15.395秒。

该新附件只关闭旧small-top有限81点中的重编号/共面重三角化来源认证缺口，原数值差分与原失败报告原样保留。`interval_neighborhood_certification=false`、`optimization_admission=false`仍保持；不转授Heart、面条、新构造策略或完整正式Cherchi准入。当前上述作业与离线队列均已终态；无新的完整Heart/面条优化交付。

这是明确的实际源面覆盖缺口，不能凭切换paper标志宣称修复。下一优先项是结合固定法向面和理想边根来源审查多重交点裁剪/连接，继续保留旧失败和显式profile，不静默切换、不调轴/权重/离散等级绕开。

### 2026-10-13 07点后：正确集合语义的小范围构造验证

首先复核队列均终态，无残留大型进程，未恢复旧Heart或Noodles顺序长跑。按最新更正，原manifest明确选inner-ellipsoid-candidate。Ruling: 在独立输入构造CLI中调用相同Epeck/libigl UNION，先验证平衡并集与单次内域切除；不修改求解器后端、不改变部件、坐标或冻结请求。这个组织方式只在全部集合/几何准入重新通过后才可用于正式输入，不能凭并集代码存在授予准入。

`exact-input-union-build-job-01`已独立编译、链接并实际运行。两个真实解析箱体测试先在旧差集探针缺少union入口时失败，候选通过；包含精确体积12、52、216、164以及双向差集为空的手算/集合对照，验证fixed=outer-design保留外壳。新增平衡四部件测试先观察缺实现失败后通过，实际并集与解析长方体双向差集为空、精确体积20。新CLI和旧探针不同文件，旧冻结模块/探针未覆盖。三个局部真实测试通过；不是全项目或案例完成。

当前串行作业`noodles-balanced-eight-part-pilot-job-01`：只取原manifest固定顺序的前8部件，使用正确inner，保持原坐标精确有理表示，比较平衡并集后单次差集与原顺序差集的实际集合，再构造outer并集及fixed=outer-design。所有操作/结果逐次封存，24GiB/swap0/60分钟。它是预登记小范围构造诊断，不是完整81部件、输入准入或优化结果。运行期间不修改其绑定依赖、不启动其他大型作业；完整面条构造仍未重启。

该8部件试点已终态：平衡并集加单次内域切除1.483799秒、相同输入顺序切除4.125141秒，两份实际CGAL结果双向差集为空；outer/fixed/design相交为空。97.496秒总作业还包括两次精确集合比较等诊断，不能用总时间当作正常构造成本。`noodles-balanced-eight-part-audit-job-01`复用已冻结独立有理边界CLI：7092面设计域closed、edge/vertex manifold、无精确自交或零面积，66319次精确pair tests。只适用于8部件，不转授完整81部件或最终物理材料。

本目录回归第一次运行误选旧newton-capture native，45测试中的Heart两份projective回归重现旧多重交点拒绝；这是加载路径配置错误，输出保留，不修改测试。选择已验证exact-core-material候选后，45项局部测试实际通过15.538秒（含三个新真实集合语义测试）。未执行全项目里程碑测试，不宣称项目全部通过。

Ruling: 小型解析集合和正确内域8部件独立边界门已通过，允许一次完整81部件平衡构造来验证该明确组织修复；不是恢复旧顺序长跑或拼接超时前缀。全部叶部件顺序/覆盖显式检查，所有Boolean从原binary64输入精确编码起，不经过浮点反馈；设计域种子来自manifest。fixed仍通过outer-design构造，未漏外壳。

当前大型作业`noodles-balanced-full-input-job-01`独立执行80次平衡固定件并集、一次inner差集、一次outer并集及一次fixed差集（总83次）；没有运行Newton。等待队列`noodles-balanced-input-verification-sequence-01`仅在完整构造成功后依次检查design、outer、fixed有理边界，任一失败立即停止，不启动优化。每作业保持24GiB/swap0/60分钟，运行中不修改绑定依赖。即使三边界通过，整体分区、固定材料质量矩、容量界、适用生产导数和正式求解仍需分别完成。Heart材料连接认证与正式Cherchi封存入口仍待办，两个案例仍无完整优化交付。

本轮末最新实测：`noodles-balanced-full-input-job-01`已实际完成全部81部件及83次Boolean，正确inner种子，完整构造55.931941秒；design/outer/fixed/fixed-parts-union原有理输出和摘要已持久化。只是输入构造成本，独立边界审查尚未完成，不是面条优化、最终OBJ或性能总门通过；不得与误配outer顺序长跑比较为正式加速比。目前唯一大型作业为`noodles-balanced-design-audit-job-01`，串行队列随后仅在通过时检查outer和fixed；不重复构造、不启动模型优化。下一次跟踪优先读取该队列和三份审查终态，以实际合法性结果决定分区/质量矩与核心准入，不按旧heartbeat重放第73次或旧有理后缀。

### 2026-10-13 07:46后：面条三份权威边界通过，推进实际输入质量矩

`noodles-balanced-input-verification-sequence-01`已全部终态：design、outer、fixed三份原有理边界均closed、edge/vertex manifold、无精确自交、无零面积。独立检查分别146.987、20.788、168.582秒；实际pair tests为817601、275932、1111362。此前正确内域顺序/舍入前缀产生的非流形与派生失败全部保留，不把新结果写回旧报告。新结果只关闭这套完整原始输入的三个权威边界门，尚非完整分区、派生网格、优化或碰撞准入。

新增`capacity_from_fixed_moments`复用独立积分的固定材料，不再为容量计算重复积分fixed。保持原有精确半空间积分公式；先在真实箱体观察缺适配实现的两项失败，再通过手算±2容量及非零原点平移(-24,0)测试。仍只给必要连续子集容量界；未排除不可行记UNKNOWN，不授予可行证书。三维轴不为z时当前输入worker明确拒绝。

`noodles-balanced-input-moments-job-01`当前运行：绑定三份材料与已完成的独立审查摘要，核验outer实际只有一个边界连接分量；原inner减原outer实际为空；对design、outer、fixed各执行一次独立GMP精确十质量矩积分，随后核对分区十项精确相等、真实SI质量/质心/完整惯量及预登记目标原点下的容量界。首次launcher预检路径误指contact_matrix.py所在目录，身份取证前拒绝，未创建作业/发生几何运算；定位实际Heart目录后启动同一个新作业，未改科学配置。当前阶段design exact moments；大型作业仍只有这一项，不重复构造或启动优化。GDB只读附加因ptrace权限失败，没有取得栈，不猜测内部热点或据此改算法。

为后续统一派生装配准备`assemble_fixed_boundary`：由同一outer和反向design界面组成fixed，避免各自舍入内部边界。实际解析盒测试先失败后通过，外实体64、设计域8、固定56且native精确分区通过。该函数只是表示装配，不自动判定包含/距离或派生准入；大型Noodles尚未调用它，不称已有合格OBJ。新增三项局部真实测试实际通过；没有运行完整项目回归。

下一步先读上述质量矩作业终态；若成功且容量未排除，再完成同边界派生出口及独立包含/分区检查、适用生产映射与正式入口。Heart材料连接认证和正式进程内Cherchi材料封存仍待办，不转授small-top有限81点证明。不恢复旧Heart优化/dispatcher或旧面条第73次、全有理顺序长跑。两个案例仍无完整优化交付。

### 2026-10-13 08:16后：精确积分进展与共享界面派生出口排队

读取最新账本和实际进程，原质量矩作业仍在原预算内执行，没有新失败。design精确十质量矩在1058.196秒完成，对应实际输入可优化域质量31.710383g；outer在1140.235秒完成，对应44.348119g；fixed开始于1140.280秒。这些是原始输入物理量，不是优化后质量，也不是论文作者模型参数。尚未运行Newton。所有已有构造/边界审查均不重跑。

已准备`export_balanced_shared_input.py`和等待队列`noodles-balanced-shared-export-sequence-01`：只在当前作业真实成功、报告存在且容量未排除后启动一次独立bounded SI派生出口作业；容量排除、积分失败、超时都不启动。它分别有界转换outer/design，并从同一份设计域反向界面组装fixed，保持公共坐标和三角形完全共用；fixed几何/体积界取两边界的组合保守上界。25mm到米按声明的1/40有理单位变换，再派生binary64，L采用权威outer包围盒对角线的保守下界，不放宽误差。

出口仍要求OBJ17位重载身份一致、实际native精确区域/包含分区、独立质量矩≤1e-6；所有失败记录单独写入且停止，不把candidate文件称为可用网格。成功清单也只登记合法派生输入，solver/collision admission仍false，不能冒充优化/仿真或制造结果。使用已通过手算/精确分区测试的共享装配函数；两个新编排脚本仅完成语法验证，完整案例出口尚未执行，不宣布通过。出口作业真实启动时冻结worker、表示/积分实现、native、三材料和当前审查/积分记录；不修改当前运行作业依赖。

当前唯一大型作业仍是`noodles-balanced-input-moments-job-01`，等待出口队列不做几何运算。下一轮优先读取此作业及共享出口队列终态；源输入、独立开发验证、求解和交付耗时继续分栏。Heart和正式Cherchi入口待办不变，两个案例仍无完整优化交付。

### 2026-10-13 08:46后：质量矩完成，派生出口的过保守见证修复

`noodles-balanced-input-moments-job-01`在约2921.86秒内真实完成，未超60分钟：fixed独立积分完成于2841.959秒，质量12.637736g；outer44.348119g、design31.710383g。十项原有理质量矩均满足outer-fixed-design精确为零；容量为UNKNOWN_NOT_EXCLUDED_BY_FIRST_MOMENT_BOUNDS。质量单位/密度只施加一次，完整惯量已保存。约48.7分钟是独立精确输入验证成本，不是Faster普通求解计时。没有Newton或优化通过声明。

等待出口的job01实际失败且保留：outer的`BOUNDED_CORE_REPRESENTATION_ERROR_EXCEEDED`，尚未写OBJ或进入native区域检查。只读`diagnose_shared_erased_patch.py`重建同一SI转换，发现12个collapsed面、0个opposite面；旧点见证把最长约5.824mm的细长面计成边界误差，而同一保留网格中真实存活边的精确投影L1上界约7.302e-19m，原允许距离8.381e-8m。它说明点见证不适用这批薄面，不说明派生网格整体已合法。

Ruling: 新增显式`SURVIVING_EDGE_V2`，旧默认SURVIVING_VERTEX_V1和旧输出行为保留。只接受确实存在于保留边界中的边，精确计算每个原面顶点到该段的有理投影L1上界；距离到凸段的凸性覆盖整片原三角形，反向段由原面对应端点及原坐标舍入界覆盖。额外胶囊体积保守界8*r²*(L1边长+2r)计入原扫掠体积界，不放宽10^-6L/10^-6L³；拓扑、包含、分区、质量矩及独立准入仍不可省略。策略进入派生表示身份，不改变权威材料或求解映射。

手工闭合箱体的真实细长三角分割测试先重现缺少新政策，修复后通过；旧点政策仍拒绝此长薄面，极小长度预算仍被拒绝。首次测试fixture的负坐标在现有GMP→binary64转换下未实际发生别名，纠正为正坐标平移后才观察到目标失败，没有改验收期望。原6项表示回归及新2项测试实际通过；旧第73次记录不升级为合法输入。

`noodles-balanced-shared-export-job-02`明确选择新政策，重用已通过的原三材料、质量矩和容量证据，不重做Boolean、输入构造或精确质量矩批次。outer/design有界转换及共享fixed装配均实际返回，outer OBJ17位写入/重载身份、native精确区域及独立质量矩比较通过；当前运行design的native精确区域检查。唯一大型作业为此job02，24GiB/swap0/60分钟；其依赖冻结，不另启大型作业。job01及旧等待队列保持失败终态。job02尚未通过完整派生准入，不把candidate OBJ交付为合格网格；成功以后仍需适用生产映射和正式求解，两个案例仍未完整优化。

### 2026-10-13 09:16后：面条完整输入OBJ通过，首次实际生产状态暴露新源面分支

`noodles-balanced-shared-export-job-02`已真实终态成功，约884.16秒总作业/880.39秒worker：outer、design、共享fixed三份米制OBJ17位重载身份相等，三实际精确区域均通过；独立质量矩最大归一化误差约1.27e-17，native完整分区通过。MemoryPeak约8.20GiB、swap0。`manifest.json`绑定三份文件、权威材料及输入质量矩摘要，明确INPUT_NOT_OPTIMIZED、solver/collision admission=false、制造资格未验证。文件名candidate保留旧历史命名，不假装已有优化结果。这一输入出口门实际关闭，不重做约48.7分钟的原有理积分。

新增归一化适配复用现有bbox中心/max span、R=density约定；真实固定网格积分提供归一化固定质量矩，不能读取手填壳矩。平移解析箱体测试先失败后通过，检查目标点、M=.875、三个二阶矩31/384、密度只施加一次及单位回转；只验证适配公式，不授予整个生产映射。`noodles-balanced-request-freeze-job-01`实际完成约2.81秒，冻结两种子原八次尝试/60次求值协议、Faster惯量＋线性质心高度、归一化权重1、level2、显式CGAL及paper-exact source/material政策。当前是已批准core-first的BINARY64_CORE_V1，不冒称A线的MPFR1024已贯通。

为固定状态取证而不是启动完整长跑，`noodles-bound-seed0-production-probe-job-01`只允许运行seed0第一初始化的一次真实生产求值，诊断预算1×1单独登记；实际随机初始化器只依赖mode/seed，首个系数没有从旧成功/答案读取。约28.70秒总作业后返回GEOMETRY_EVALUATION_FAILED，1次求值，0次实际优化更新。失败发生在源面构造、尚未Boolean：`UNSUPPORTED_PAPER_IDEAL_FACE: unanchored cone crossings at face 1`。该状态、系数、阶段、请求、native身份及失败日志均保存；进程/等待队列均终态，当前无大型作业。

`test_noodles_seed_source.py`从这份真实记录重现同一源面失败，作为尚未修复的RED回归；不改目标/轴/网格/seed、不用legacy回退或另一成功初值绕过。源码1121–1200附近显示：现有正Lorentz锥路径只接受单根连接一内一外角点，或两个内角点属于相反锥的双根；没有覆盖两个外角点之间存在有限正区间的双根，也要求每个裁剪patch有原有限角点。下一固定修复应补齐这种根归属/有限patch，并同步记录必要有限锚点的构造来源与运动导数，不凭近邻配对、放宽阈值或仅改变显示网格。新映射需重新验证受影响源面/质量矩导数；现有Heart两份成功数值门和small-top有限81证明不转授。

这是一个已登记初始参数状态的曲面支持失败，不是输入几何非法、整个面条不可行、Boolean后端失败或完整优化失败。现不扫描其它方向/步长或启动一小时相同求解，优先定点修复该新实际源面回归及正式核心入口。Heart生产分支认证与完整求解待办仍保留。两个案例仍无完整优化结果、独立重复或真实动力学交付。

### 2026-10-13 09:46后：面条有限区间源面候选与未通过的导数门

针对同一seed0初始状态，固定法向面检查定位四个双根面：两个外角点之间存在有限正区间，另两个面没有原始有限角点。候选补齐该区间的根归属，在同一原球面边上建立共享有限中点锚点；锚点的法向随系数运动，使用显式`-3`边来源并将隐式根导数、法向运动、逆Gauss映射和远端封口运动传入source Jacobian。未改变目标、轴、seed、level2或输入材料，不读取成功答案。该构造改变了有效映射，必须生成新证据；有理源面政策对尚无区间包围实现的新锚点明确拒绝，不宣称MPFR/rational覆盖。

`projective-interval-build-job-01`独立构建后仍失败：两根patch闭合弧重复插入同一锚点，造成边关联错误。该失败及旧二进制保留。修复每个patch的单次插入后，`projective-interval-build-job-02`独立构建约27.80秒，实际首个冻结seed源面及运动回归通过；核心目录50项局部测试实际通过15.723秒。旧冻结模块未覆盖，这不是全项目或优化完成。

`noodles-interval-source-derivative-job-01`实际完成81个源面求值，十方向、三个原登记尺度3e-6/1.5e-6/0.75e-6，保存全部源坐标、运动和构造来源；约20.20秒。连接、锚点/射线来源及截断半径在81点一致，基点精确边界闭合、定向、边/顶点流形、无自交或零面积。但30项源坐标差分中两项未达标：direction1和9在0.75e-6尺度分别1.40658e-6、1.03526e-6，门槛仍1e-6，`source_checks_passed=false`。进程正常退出不等于科学检查通过；production_moment_admission及optimization_admission均false。

只读保存点诊断将两项最大误差同时定位到vertex93、射线来源`[217,40,135,0]`，它不是新增有限中点自身；坐标约(0.7224,-20.9814,36.1935)。原三尺度误差不随缩步单调下降。对原保存坐标先减基点并用long-double求和，误差仍1.40487e-6/1.03505e-6，说明不能仅靠差分累加精度解释或修改报告为通过。该离线检查没有新源面或Boolean求值，没有改步长，全部原失败保留。

当前上述构建与源差分作业均终态，无大型作业运行。下一项应固定该射线的根、归一化方向与远端平面交点进行算术/解析运动对照，判别source构造舍入和导数问题；仅有确证修复或预登记定向诊断才重测。源门未通过前不安排该候选的完整Newton长跑。Heart旧材料分支认证、正式Cherchi入口及两个模型完整优化/重复/OBJ/动力学交付仍未关闭，不恢复旧队列。

### 2026-10-13 10:16后：有限锚点算术定位、显式代数源策略及首次材料检查

`noodles-interval-ray-arithmetic-diagnosis-01`只读取原81点中direction1/9的三尺度位置，使用同一C++固定法向生成器的edge40/135，独立100/180十进制位重建两根、根中点有限锚点、理想方向与球面截断。没有新源面或材料Boolean。高精度三尺度射线导数与基点生产Jacobian相对差约3.80e-12/3.87e-12；100/180位基点射线差约8.25e-100。基点有限锚点坐标与旧构造差最大8.58e-12，射线最大2.39e-13。这个独立模型支持有限点/法向计算舍入问题，不能作为生产映射差分通过证据。需注意ray93确实依赖新增moving-normal锚点92（normal id217），不是原固定法向有限点。

现有`algebraic_finite_map`此前在该冻结请求为false，明确选择true进行九点定向对照，不静默更改原请求。`noodles-interval-algebraic-arithmetic-job-01`以相同失败位置、原最小尺度完成9个真实源面求值：两方向全部源点最大误差1.009e-7/1.097e-7，约2.53秒。这是已有代数有限点入口对同一构造的算术修正，不改初值、法向level、截断预算、轴或目标；旧false profile失败原样保存。

在该明确算术修正基础上，新登记`PAPER_INTERVAL_ALGEBRAIC_FINITE_SOURCE_DIAGNOSIS`，`noodles-interval-algebraic-source-derivative-job-01`重新执行完整81个真实源面求值，原十方向/三尺度/五点及存储参数位置保持，未事后缩步。约22.84秒，30项全源点误差最大2.511e-7，精确源边界检查通过，连接/构造来源/半径稳定，source_checks_passed=true。只关闭该新profile、固定Noodles状态的源坐标数值门；production_moment_admission及optimization_admission仍false，不继承旧profile或授予MPFR/rational支持。

当前唯一大型作业`noodles-algebraic-first-material-job-01`：使用同一合法输入和冻结seed0首次系数，显式代数有限点政策及projective-interval-build-job-02，对真正CGAL生产材料执行一次映射，保存质量矩/Jacobian及耐久阶段事件。诊断预算1×1明确标记，不是完整求解、不要求一次求值收敛；原八次尝试/60求值的完整协议未改变。仅source门通过后允许此定向材料取证，后续仍需适用生产质量矩差分、材料来源/分支及完整运行验证。作业24GiB/swap0/60分钟，运行期间不得修改其绑定依赖或另启大型作业。Heart原材料分支认证待办不变，两个案例仍无完整优化、独立重复或真实动力学结果。

### 2026-10-13 10:46后：面条首次真实材料返回，生产质量矩检查启动

`noodles-algebraic-first-material-job-01`真实终态：约63.14秒，MemoryPeak332378112 bytes、swap0，1次CGAL材料映射返回2877顶点/5870面及有效生产Jacobian。原第一个seed0系数保持，源构造、实际Boolean、材料边界、transport及积分均实际执行，全部阶段和数值附件保存。Newton事件包含initial_evaluation_succeeded和正规系统payload；最终RESOURCE_LIMIT_REACHED来自预登记的1次求值诊断预算，不是60分钟超时或完整优化不收敛结论。0次实际接受更新，optimization_admission=false；这是一份诊断材料，不能作为优化终态、最终OBJ或仿真输入交付。

已登记并启动`noodles-algebraic-production-moments-job-01`，仍选同一projective-interval-build-job-02及显式algebraic finite政策。清单直接绑定已冻结source81点的相同系数、十方向、原三尺度和五点位置；重新执行81次真正生产CGAL映射，不把source坐标门转授质量矩门，不拼接单点诊断。每点保存实际材料/Jacobian，重新读取文件以独立GMP四面体积分；固定材料只独立积分一次、核验目标坐标/密度后使用，不能读手填质量矩。逐阶段、开始/失败/返回事件立即持久化。全部数值达标后仍需适用材料连接/来源认证，报告预置full_material_branch_certification=false、optimization_admission=false。

当前唯一大型作业为该81点生产检查，输入与worker/native/积分器/差分权重及全部源点摘要冻结，24GiB/swap0/60分钟。不能在运行中改依赖或启动另一个大型作业。没有启动完整心形或面条优化、重复或动力学；Heart已有数值通过记录不重跑，正式Cherchi生产封存和Heart材料连接门仍待办。

### 2026-10-13 11:16后：面条生产数值门完成、同次CGAL来源附件

`noodles-algebraic-production-moments-job-01`实际完成81点、约1236.87秒（20.61分钟），MemoryPeak435429376 bytes、swap0。十方向、三尺度、十质量矩差分全部通过，最大误差5.414e-11；保存后独立GMP积分最大差5.551e-17。所有实际几何为2877顶点/5870面，但索引连接有47种身份；数值通过不能直接授予材料分支或完整优化门。原失败source批次保持原结论，不能把新代数profile通过写回旧false profile。

当前CGAL FullRank native对象已经保存同次Boolean的exact material坐标、birth faces和source faces，但Python入口对非rational-source运行丢弃附件，导致已完成81点没有这些必要取证。新增显式`capture_material_proof=false`参数，默认返回内容/数学/哈希完全保留；true仅导出已有同次字段，不能重做Boolean、积分或新导数。binding-only候选`cgal-proof-capture-build-job-01`约20.02秒独立构建/链接，全部核心对象从projective-interval-build-job-02按摘要保留，无结构布局变化，旧native未覆盖。

真实回归先在旧模块观察opt-in入口缺失失败；候选首次测试采用legacy material配置，exact tape为空而失败。定位现有精确tape只在exact source/material协议收集，测试改为显式exact material配置（正式Noodles本来已开启），不改数学门槛或将空tape当作来源。随后真实附件回归通过，默认和捕获调用的全部原科学字段逐值相等，exact坐标/出生面覆盖实际输出。核心目录51项局部回归实际通过17.162秒；未执行全项目测试或声明完整核心/案例通过。

当前唯一大型作业`noodles-same-call-proof-capture-job-01`：只为已有81点的基点及首个不同连接点补两次必要同次CGAL取证，绑定原清单/请求/材料和隔离binding候选；必须与原保存的全部科学字段相等，否则拒绝使用。不是重跑整套数值差分、不启动完整Newton。随后利用精确出生面、坐标与来源分析重编号/共面重三角化，不能用近邻匹配替代证明。资源仍24GiB/swap0/60分钟、串行，不修改运行依赖。Heart既有数值证据不重跑，来源分支与正式Cherchi接口待办保留；两个模型仍无完整优化交付。

### 2026-10-13 最新用户批准：案例先行、后端收尾并行

用户批准《三维Faster核心最短收尾计划》，明确CGAL路径满足必要门后立即用于案例，不等待Cherchi全部工程收尾。核心正确性、核心性能、案例优化/OBJ/动力学分别验收；不恢复旧长跑、旧81点/Nef队列，不自动提交、合并或删除。执行采用现有工作区与本账本；隔离native候选不覆盖历史二进制。

`noodles-same-call-proof-capture-job-01`实际成功终态73.122秒，两点全部原科学字段逐值一致；先前“仍在运行”的记录已过时。新增精确birth-face来源约束重建，4个手算/反向测试先RED后GREEN，覆盖真实边面交点运动、原输入顶点、1e-60错误出生面及秩不足拒绝。`noodles-offline-proof-pair-job-01`仅使用保存源面/精确材料，0次Boolean/源面调用，约6.305秒证明两状态各2877顶点（1799显式顶点、1078边面交点）精确重建，并在各自固定几何上证明有向表面链一致；仅两点，不授予整批/优化完成。

`noodles-saved-lift-triage-job-01`继续只读原81点，约168.878秒完成所有原源面系数/科学身份、精确来源重建、唯一binary64视图匹配及重建表面链比较。全部通过，但79点原记录缺少实际native精确坐标，报告明确`original_exact_material_certification=false`。Ruling: 浮点视图与精确重建匹配不能替代未保存的实际Boolean有理输出；只补缺失原生附件，不重新计算FD权重或独立积分整批，不以两点推广全部。

公共`evaluate_spin_it_faster_fullrank_quadric`新增显式`capture_material_proof=False`，默认返回/asdict保持；true传递已有同次Boolean附件，以不可变canonical JSON保存并提供fresh-view属性。检查精确策略、geometry身份、坐标视图、birth-face覆盖；附件不是材料/分支准入。真实公共测试先观察缺入口失败后通过。`DurableFullRankMap/run_durable_fullrank`新增同名显式选项，逐次材料附件立即保存；终态按实际系数和最后接受状态几何身份匹配已存求值及摘要，不盲取最后试探、不重新求值。终态引用依旧标记材料/导数未独立验收。真实映射及错误终态、文件修改拒绝测试先RED后GREEN。

本目录58项回归实际通过19.252秒。首次运行3个精确并集测试因未提供`EXACT_BOOLEAN_PROBE`而报配置错误，保留该事实；指向已有exact-input-union-build-job-01的探针后重测全部通过。未运行全项目基线，不宣称全项目或核心收尾完成。

**当前正在运行的串行队列：**`noodles-native-tape-certification-job-01`复用已有两份native附件，仅对其余79点补抓同次精确坐标/birth faces。每点必须与原科学字段逐值一致、精确来源重建为实际native材料的双射、有向表面链一致、source构造分支一致；失败即停止，不缩步、不升级旧结果。24GiB/swap0/60分钟。其依赖（spin.py、birth_face_recipes.py、exact_planar_patches.py及隔离native等）已经冻结，运行期间不得改动。

`short-closeout-sequence-01`只等待该作业成功终态，随后运行一个`small-cgal-proof-closure-job-01`：使用现有小型top的seed1请求、真实CGAL Newton、同次终态附件；通过已有MeshSnapshot、partition和spin_bundle的Writer/Reader完成设计域减同一保留材料的空腔构造、五角色材料写入与独立重载、质量矩及分区检查。新worker仅完成语法检查，实际闭环结果尚未产生，不能预先宣称通过。等待队列绑定下一个worker摘要，修改即拒绝；没有自动启动完整Heart/Noodles作业。后续先读两个终态报告，满足准入后再按原种子0/1和8×60协议开展就绪案例。正式Cherchi接入与性能门继续待办；心形与面条均无完整优化、独立重复或真实动力学交付。

### 2026-10-13 12:33后：有限81来源门及小型正式闭环完成，启动面条完整求解

`noodles-native-tape-certification-job-01`已实际完成全部0–80点，约1193.261秒；复用两点、补79次缺失原生附件，不重新跑差分积分整批。全部原科学字段逐值一致，实际同次有理材料与来源重建双射、有向表面链及源构造分支检查通过，`complete_frozen_sample_certification=true`。只覆盖登记的有限81点，不宣称任意参数邻域可微，不改旧报告中尚无证据时的false结论。

`small-cgal-proof-closure-job-01`正式闭环实际通过：原seed1真实Newton约82.546秒，worker总85.379秒；原始残差2.5667963906969476e-15，独立总质量矩最大差5.329070518200751e-15。终态从同次已存Boolean附件选择，0次终态quadric重求值；保存后fresh Reader重建，设计域/外域完整精确分区均通过且无正体积重叠，有实际接受更新。wrapper约86.006秒、MemoryPeak91373568 bytes、swap0。它是小型接口回归，不是复杂案例成功、重复或性能目标通过。

原串行验证队列已经终态且无残留重作业。新增`noodles_core_run.py`将冻结请求/输入、原81点数值报告、补充精确分支证书、源门及小型闭环组合成一次明确的登记运行准入；旧输入和局部报告不改标签。4项准入测试先观察缺实现失败后通过，覆盖完整8×60协议、两个原种子、缺失样本、旧请求/清单、错误构建及数值失败拒绝。无预设quadric，不继承CGAL证据到Cherchi，不宣称MPFR1024支持。

**当前大型作业：**`noodles-registered-core-seed-0-job-01`由资源wrapper启动，原完整seed0、8次尝试/每次60求值、raw1e-12、level2、显式CGAL、paper projective+algebraic finite+exact binary64源/材料边界。使用隔离`cgal-proof-capture-build-job-01`，同次材料proof与耐久日志开启；24GiB、swap0、60分钟、数值单线程。输入、Python求值器、journal、worker、原native和所有准入报告均绑定；运行期间不改冻结依赖、不并行新重作业。后续先读真实终态/失败状态，再安排seed1或定点诊断；不无修复重跑相同一小时作业。

心形材料分支取证和正式Cherchi入口仍未关闭；当前面条刚进入完整优化，不是已有优化成功。终态独立验收、OBJ、重复及MuJoCo+Blender分别待办，初次收敛不得写为完整案例交付。

本轮退出前的实时快照：seed0已完成4次实际求值，第4次以步长0.25接受真实更新，raw残差由初始523.6990110761269降至342.4630640311879；正在第5次求值，固定设计边界proof缓存实际命中。仍远未达到1e-12，不能写成收敛或案例完成。旧journal的topology标记仍混入连续source身份，因此`topology_transition_relinearized`仅按现有求解事件原样留存，不解释为已证实真实拓扑切换；本轮不修改运行中的journal。

### 2026-10-13 13:03心跳：首轮预算终止，保留原失败并登记定点取证

面条seed0第一次尝试已到60次求值上限，原始残差最后0.25340604743808426，15次实际接受更新，未收敛。按原8次尝试协议自然进入restart1，初值残差385.11411031693655；本次快照约24.3分钟，尚未触及整作业60分钟上限，程序仍运行。没有手动重启、额外seed或并行几何作业。

`noodles-first-core-failures-registration-01/registration.json`只读原耐久日志，固定第3次源面内外符号不一致及第7次首个材料边界失败的完整系数、原异常、已完成阶段、映射/native身份和日志前缀摘要。0次几何/求解调用，未改冻结依赖。后续日志追加不破坏前缀身份；这不是失败修复或材料准入。

源码审查确认`nondegenerate=0`在exact material路径也可能是整体region无效所致，而日志`degenerate_faces=0`；因此不能仅凭最小浮点叉积为0断言真正零面积。该试探还明确出现closed=false、vertex_manifold=false。待原作业终态后，优先对第7次做同生产源面/封口/内外侧的原有理Boolean输出及binary64对应取证，比较canonicalization前后顶点link，区分原实体非流形、舍入/合并改变邻接和检查器问题；在捕获前不猜测修复、不删薄面或关闭检查。第3次源面符号问题单独保留，不与材料问题混为一因。

本轮没有启动心形完整求解或重复检查，也没有新优化终态、OBJ或动力学成果。后续首先读取当前seed0终态，再决定已登记诊断或seed1的调度，不无修改重启同一配置。

### 2026-10-13 用户新指令：暂停两案例，聚焦三维论文核心

用户明确要求先暂停心形/面条处理，先核对“完美复现论文算法”的剩余缺口。已通过应用工具将自动跟踪`automation`置为PAUSED，保留原提示及规则。心形无运行进程；向确认属于当前面条seed0作业的worker PID1658843发送SIGTERM。wrapper随后完整退出：returncode143、PROCESS_FAILED（仅wrapper技术状态）、耗时2779.510491679秒、MemoryPeak617512960 bytes、swap0。新增`user-interruption.json`明确USER_REQUESTED_STOP_CASE_WORK，保留原result/complete不改写，不能误称资源超时、研究性不可行或完整优化失败。

中止前已有24次接受更新，最后evaluation93/restart1，raw残差35.45564856153745，尚不达标；不续拼成完整成功。已经核实无面条worker/launcher或心形调度器残留。两案例不自动恢复，也不启动seed1、终态加工或仿真。保留全部输入、失败、材料及证据；三维核心的来源/几何/正式后端/性能缺口继续单独评估。此前heartbeat指令与本节冲突时，以本次用户暂停为准。

### 2026-10-13：实施10种子核心收尾——首个可复现材料缺陷已修复

最新计划存于 `docs/superpowers/plans/2026-10-13-three-dimensional-faster-core-ten-seed-closeout.md`。只推进核心；两案例和 `automation` 仍 PAUSED，不恢复任何案例队列。

新增原生 `ScopedQuadricGeometryObserver` 与 `capture_spin_it_faster_fullrank_geometry`，在实际生产调用中记录 source-before-side、boolean-before-normalization、boolean-after-normalization；异常仍拒绝求值，诊断输出始终 `production_admission=false`。线程局部 RAII 恢复 observer，不改变 settings/result 的既有布局，不替代失败现场进行另一套 Boolean。`core-failure-tape-build-job-01` 独立编译约52.4秒，未覆盖旧模块。

`core-registered-failure-capture-job-01` 已完成两个注册状态，原异常逐字匹配：第3次约31.99秒，第7次约32.27秒。所有原始有理坐标、binary64坐标、面、来源、源面法向/理想射线及正规化前后索引均保存。零次完整优化。分析记录 `tape-analysis.json`、补充精确binary64叉积的 `tape-analysis-v2.json`，旧文件保留。

**第7次根因：**实际原始 Boolean 有44856面、22254个使用顶点，每边两个反向使用，全部顶点link闭合；无坐标别名或重复精确顶点。`canonicalizeBooleanMesh` 用浮点边差叉积为零作为删除条件，误删32502号面。该面在原有理坐标、精确解释的binary64坐标中均非零；正规化后44855面、3条缺口边、3个开放顶点link。因此这是正规化判定缺陷，不能归咎于原Boolean的“真正零面积面”，也不应放宽闭合检查。

最小修复：对携带同次精确坐标附件的路径，在浮点叉积为零时用原binary64坐标的精确叉积判定；非零面原样保留及保留birth facet。真正精确退化或重复索引明确拒绝，无扰动、补面或重三角化；旧无附件路径保持旧行为。候选 `core-small-face-build-job-01`，native SHA256 `77c5bb07fe5203307a34fea8db721d59c8bbb900db6bb6a534d442459614c8ae`。

真实回归先 RED：`core-material-normalization-red-job-01`，44855 != 44856；再 GREEN：`core-material-normalization-green-job-01`，84.208秒，44856个原精确有向面及birth来源全部保持。该回归只证明正规化不再丢面。随后 `core-small-face-verification-job-01` 通过实际正式求值：44856面、`derivative_defined=true`、`VALID_FIXED_BOOLEAN_TOPOLOGY`，总84.166秒，实际结果及同次材料附件已保存；这不是完整优化、终态封存或适用81点导数准入。

**第3次定位：**精确源边界审计通过closed/edge+vertex manifold/orientation/self-intersection，0自交。内外选择失败的正梯度面固定为8、9、96、97，其余104面负。保存了各面的有限点/理想射线/区间锚点及精确梯度点积。尚未证明选择器在这些粗鞍面弦片上的适用性，不能用多数投票或局部翻面通过；下一步核对构造侧身份与选择器，不声称源面方向问题已解决。

**实测阶段开销：**同第7次新状态，输入准备至source begin约26.80秒，source约3.46秒，Boolean约1.48秒，材料检查约13.59秒，首次transport设计检查约21.06秒。出现两次source motion/transport，feature assembly分别约8.04和8.09秒；源码FullRankMomentEvaluator的canonical profile确实重做导数传递。此为可定位的后续优化项，不把不同范围或冷启动耗时当作论文性能。

验证：诊断版本63项核心目录测试21.731秒；修复版本64项快速回归21.527秒（明确排除已经单独完成的约84秒注册几何回归）。最初局部通配测试曾缺 `_cherchi_core_boolean` 路径而失败，配置正确独立模块后复测通过；不隐瞒、不计为算法失败。独立只读审查 `/root/review_core_failure_tape` 未发现本批生产阻断，提出零点积分号分类和注册输入摘要核验两项小问题；均修正，零点积分类测试先RED再纳入快速GREEN。审查范围不涵盖全项目、真正退化所有类型、Cherchi正式闭环或10种子门。

**当前串行作业：**`core-small-face-independent-moments-job-01` 正在从保存材料做独立精确四面体积分（binary64及原有理两种表示），不重做Boolean，不运行Newton。保留24GiB/swap0/60分钟限制；结果未出前不授予积分通过。正式Cherchi、受影响导数/分支、三个规模各10种子、独立重复及性能门仍未完成。不得据本节恢复两案例或其自动跟踪。

上述独立积分随后已终态通过，实际40.299秒；binary64表示和原有理表示分别独立积分，与生产变量质量矩的最大误差均为 `9.298117831235686e-16`。仅授予该固定材料的独立变量质量矩检查，不授予生产导数、完整优化、重复或封存准入。该作业已结束，无新案例队列。下一批优先关闭第3次源面选择证据，再处理正式Cherchi与同次数据复用；原10种子和性能门均保持待办。

### 2026-10-13 最新顺序：正式Cherchi优先，本批实际结果

用户已将顺序改为正式Cherchi→几何/复用→三规模10种子→性能审查。最新详细账本为 `docs/three-dimensional-faster-cherchi-first-progress-20261013.md`，覆盖前文旧调度。

正式Python入口及独立bindings/cherchi构建已接通；小型实际种子求解、同次终态封存、独立重载及完整进程重复通过，raw2.04e-14，18次求值，两次约54.4/50.8秒。全部科学事件和材料一致，未达性能门。

新物理Boolean映射81点独立积分差5.55e-17、有限连接证书通过；登记三尺度五点导数门仍失败（最大5.63e-5）。保存节点九点诊断误差约5.63e-12、原误差呈四阶缩放，仅支持截断误差解释，不升级原门。保留原协议，不启动30种子。源面侧选择/复用和必要反向测试仍待办。所有本批作业已终态，无案例队列；心形、面条与automation继续暂停。
