# 五目标与两个复杂案例收尾执行账本

## 两项阻塞修复接续

- 已实际读取仅惯量job终态：RESOURCE_LIMIT_WALL，elapsed3600.117543秒，MemoryPeak171331584 bytes，swap0；第65次Boolean未返回。随后systemd状态inactive/dead、MainPID0、ControlGroup空，确认旧worker退出；未启动质量加权高度后续（原控制器遇失败停止）。
- 面条job03已实际准入，worker3314075；前缀验证0.5134秒、后端/冻结输入预检0.0974秒通过，进入locked_problem_from_frozen_input。不是只登记启动意图。下一步读此job真实输出，失败定点诊断，成功继续seal/OBJ。
- 诊断停止机制修正：Python BaseException跨pybind后仍是std::exception，会被Newton吞成普通失败；实际重放worker改为取证和报告durable publish、journal.close后os._exit(0)，终止隔离诊断进程。不能把此PROCESS_EXITED0当优化成功，报告明确optimization_complete=false。
- 接受64二进制系数的独立精确有理数临界值=-1.1660117297668108797977182501034386e-13（非零），A determinant约-0.01513791055。NumPy独立下一步预测arc5.67814e-7且临界值浮点舍入0，只是下一步候选诊断线索，不冒充实际65输入或根因。

- 新局部线索（仅接受64，不是试探65）：按literal monomial约定A对角q4/q5/q6、非对角q9/q8/q7的一半，A特征值约[-0.68064822434,0.14912717887,0.14913734821]，中心约[0,0,0.5]，临界值c-b^T A^-1 b/4约-1.17e-13。提示接近锥面退化边界，应对实际65捕获后核对；不得因小阈值直接拒绝或据此宣称根因。
- 重放取证包补充从实际q65与保存Newton基点计算dual-sphere弧长，注明独立几何测量而非native直接步长事件。运行依赖仍未覆盖，尚无实际取证结果。

- 本目标回合归类progress：新增真实保存操作数重放worker与受限原生栈采集launcher。先校验packet canonical SHA及实际bridge二进制身份，再原样调用Boolean，不运行Newton/不切后端；返回精确输出存入新诊断job目录，不写回原capture目录。
- gdb从启动时跟踪子进程可用；在/bin/sleep 30上20秒定时SIGINT实际捕获clock_nanosleep/nanosleep栈并终止诊断子进程，退出0。此结果只证明采集机制，不是Cherchi热点证据。perf_event_paranoid=4，未修改安全配置。
- 会话55033等待capture-job-01成功result后启动inertia-boolean65-stack-job-01；捕获失败则明确停止，不运行替代输入。仅资源占用可以重试准入。调试样本20秒后终止是定向取证，不改变生产60分钟上限/回溯规则。
- 新真实参数转发与packet篡改拒绝测试通过；联合11 passed/0.28秒。前次仅惯量worker再次确认38分37秒CPU99.6%。面条job03和capture尚无result，不能授予完成；目标保持active。

- 后续新增replay_inertia_boolean65.py及独立资源launcher：用相同native13/bridge02/arithmetic03、同seed1和原协议重放；RunJournal逐项核对旧NEWTON_EVENT科学投影与所有返回材料。第65次进入实际bridge前保存literal操作数、系数、请求、来源运动与后端身份。诊断停止使用BaseException，避免被Newton普通失败回溯吞掉；不修改正常求解器。当前尚未实际运行或授予前缀一致。
- capture目录必须先创建（RunJournal要求现存真实目录）；静态语法检查和两项桥代理测试通过。真实取证还需验证完整执行路径，调用前回溯步长尚没有原生直接记录，不冒称字段齐全。
- 面条控制器会话33626；取证控制器会话43689，后者先等待noodles-terminal-recovery-job-03/result.json，再通过项目锁及unit/cgroup准入，仅重试资源占用，不重试失败作业。实际solver2915886复查34分56秒仍在计算。目标保持active，尚无新的材料包、热点或修复完成结论。

- 修正离线恢复启动器：新目录noodles-terminal-recovery-job-03，移除classified-boundary，保留reuse-input-proof/owned-terminal；加载native12匹配compiler-runtime。旧失败目录不改写。
- 恢复worker在大型输入构造前验证实际Cherchi/arithmetic模块身份及冻结输入文件摘要；仍从已验证接受前缀恢复，不调用Newton。原生生产检查器要求保持。
- 新增独立诊断桥CapturingBridge：在实际Boolean调用前通过调用者耐久publisher保存原样操作数及上下文摘要，支持登记调用序号前停止。当前仅工具和轻量测试完成，尚未接入实际65次重放，也没有热点证据。
- 配置、输入/后端篡改拒绝、接受前缀及取证工具共10项测试通过0.26秒；git diff --check通过。不能当作实际面条封存通过。
- 仅惯量worker2915886检查时仍在第65次Boolean，运行31分20秒，CPU99.5%；没有修改绑定依赖。新恢复控制器仅在PROJECT_JOB_BUSY/PROJECT_WORKLOAD_ACTIVE时等待，实际准入仍由项目锁和unit/cgroup检查决定；作业本身失败不重跑。
- 当前面条材料包、65次输入实际捕获与Boolean热点/修复均未完成。准备时间单列，不阻挡正确性；未拆分时间不倒推论文对应耗时。

## 批准范围

执行聊天批准的《高效收尾规划：统一 Faster 核心、完成五类功能与两个复杂案例》。正式生产采用当前核心精度与精确几何，高精度独立验证；旧 MPFR1024 证据保留。本轮不得自动提交、推送或恢复历史队列。

复用当前 spin-it-strict-computational-v2 工作树。既有未提交报告保持原样。大型作业继续项目互斥、24 GiB、swap 0、60 分钟；冻结运行依赖不得覆盖。

## 本轮实际实现

- 新增不可变 Spin2014ObjectiveContract，绑定五目标选择、权重、单位、模式、质量矩顺序和新精度配置。复用已有独立精确有理数目标/梯度/约束公式，不重复实现原公式。该对象暂未接入生产 Newton/终态封存，不能算五目标正式闭环通过。
- 七项测试通过：手算惯量比25/8、平方高度项4、质量加权项16，验证五类原目标值与约束数量、长度单位变换及模式失配拒绝。首次运行缺少新模块失败，随后通过。
- 新增 MeshSnapshot.interior_difference：使用实际 owned outer/design 原生对象构造固定材料，并从同一个结果导出快照，避免固定材料导出后再次重建 Nef。仍要求完整输入分区检查，旧默认构造不变。新身份不冒用旧输入摘要。
- 新的真实原生测试先因缺少方法失败；实现后，在 native12 与匹配 libstdc++ 下联合测试12 passed。实际精确差集、独立重建集合相等、外部/接触输入拒绝、真实完整分区均通过。尚未测量面条实际输入成本，不授予速度通过。
- 扩展相关验证：目标契约、共享差集、owned快照、终态封存配置、保存操作及重载元数据六组共27 passed/5.51秒；git diff --check通过。既有pytest-asyncio配置弃用提示保留。未运行整个项目测试，不能声称全项目回归通过。

## 接续顺序

### 实际面条输入作业

native13已实际PROCESS_EXITED/0，编译/链接/加载收据完整。新原目标入口模式负例及契约9项通过。spin2014-shared-small-lowering-seed-1-job-01真实CONVERGED，19求值18更新，raw2.5098699634155215e-14，正常5.416秒，总14.033秒，独立原目标终态、材料封存及保存重载通过。独立repeat正常5.461秒，总14.068秒；spin2014-shared-lowering-repeat-comparison-01.json八项严格比较全部true。仅此小型原目标平方高度策略闭环，不转授五目标冻结输入或所有性能门。

随后实际运行TOP_INERTIA_ONLY同seed1：已进行一次重启、64次返回求值，最近接受raw3.9683818973233975e-5；第65次求值source已返回，Boolean长时间未返回，worker持续单核计算。不是首个试探（早先commentary定位已纠正）。当前不能判为不可能收敛，也不能以已有残差下降授予成功；原资源上限保留，运行依赖冻结。会话15747按顺序在本项成功后执行TOP_MASS_AND_LOWERING，失败则停止后续；须读取实际result后继续处理，不重复启动。

当前明确速度缺口：面条实际输入主要阶段536.18秒，余下正式准备/终态尚未计入；仅惯量小型路径某个Cherchi Boolean已超正常单次速度门。应固定该试探源面与实际输入后定位，不能以小型平方高度策略5.42秒代替其他目标性能。

接续已取得实际终态：noodles-owned-input-probe-job-01 PROCESS_EXITED/0，独立partition268.564秒，四个主要阶段合计536.180秒，完整精确分区通过、冻结共享固定边界一致。冷输入总体仍须包含加载和正式准备，不能称完整600秒门通过。新表示fixed_identity已单独保存在report，尚未接入历史终态恢复。

正式原目标继续接入：Python external moment wrapper增加显式原目标配置与模式拒绝；Cherchi正式配置/结果绑定目标契约，实际Newton传入原适配器，终态用既有独立2014公式重算梯度、约束及raw驻点残差。静态保持独立平面入口。33项相关目标契约/终态检查点/记录/owned封存回归通过5.31秒；native13新路径仍须实际构建及运行验证。

新增spin2014共享小型真实seed1目标TOP_LOWERING_ONLY启动器，沿用原8次/60求值协议、实际Cherchi与同次材料封存/独立重载，不使用保存解初始化。须先等native13构建与新目标模式拒绝测试通过，再串行运行。当前构建正在正常推进，不授予构建成功。

接续实测：outer32.225秒、design89.597秒、owned_fixed145.794秒，当前进入independent_partition。固定材料构造旧321.714秒；该项单次收益成立，完整准备与600秒门尚待最终分区和报告。运行进程持续单核计算，未重复启动。

五目标原生链路发现并修改：native12正式构建为SPIN_IT_2014=OFF，原目标适配器未加载。新增外部材料Newton入口显式objective_settings，旋转2014配置调用既有makeSpin2014QuadricMomentProblem和solveMassMomentSystem；模式失配、静态和受限子空间拒绝，旧Faster默认路径保留。结果绑定目标profile/hash。该源码尚未构建验证，不授予通过。新测试在native12因缺少2014目标类型按预期失败。

已准备独立native13构建配置sharedvertex-slabs-spin2014，显式启用2014。控制器会等待面条输入作业实际result后，通过现有互斥运行器启动；当前未启动大型构建。控制器会话52819；面条控制器会话69488。下一回合必须核查真实result/进程，不凭会话号判定存活。构建后须先执行新目标回归，修正编译/运行问题，再接Python正式请求与终态审计；不能把源码接入计为正式五目标求解成功。

已通过项目级串行资源运行器启动 `experiments/faster_core_repair_20261012/noodles-owned-input-probe-job-01`，worker为`profile_noodles_shared_input.py`，绑定native12和冻结输入摘要。实际日志已进入outer构造。此作业比较真实共享固定边界并执行独立完整分区，生成分段耗时；尚未返回结论。运行依赖冻结，下一步必须读取该作业result及output/report，不能重复启动或凭启动状态宣称成功。该作业没有Newton调用，也不授予案例或性能完成。

1. 将原目标适配器接入原生 external moment evaluator，绑定正式请求、运行记录和终态独立检查；静态保留平面路径。现有 Faster-only 默认行为保持兼容。
2. 在面条实际冻结输入上使用上述 owned 构造，对照旧固定材料边界/实体/质量矩，登记新输入身份和分段成本；随后接入准备与终态恢复。
3. 根据保存81点的有向面片与来源继续判定心形29个不同连接样本，检查真实分支后的新Jacobian消费与实际残差下降。
4. 相关小型闭环、重复通过后冻结版本，再按计划运行两个案例及五目标；最后统计、OBJ、动力学。

当前五目标统一正式链路、面条完整材料包、心形完整优化及登记速度门均未完成。本轮未启动复杂案例长跑，未修改旧作业状态。
