# Spin-It 2014：恢复后的材料导出与独立封存

范围：现有阶段恢复接口与材料导出链路衔接；不是新的论文算法、原论文案例完成、
可行性证明或物理实验。本轮未修改目标函数、优化规则、舍入规则或验收容差。

## 接口变化

- `solveSourceAdaptiveSpinMaterialStages` 使用现有求解器的私有 terminal endpoint。
  只有 `BINARY_MOMENT_VALIDATED` 才尝试材料导出；暂停不生成几何产品。
  原来的全程求解入口委托同一实现，科学摘要规则不变。
- 原生 Python 入口新增 `maximum_new_passes`、`resume_checkpoint`；返回只读 bytes
  或 None。每次恢复新建 source/request，由现有 codec 重新核验历史状态和 trace。
- `solve_and_seal_source_spin2014` 增加相同执行参数及显式磁盘 trace 配置。
  恢复终态仍由原来的独立 `SealedSpinArtifact.create` 检查材料，不从 JSON 授予通过。
- 材料 checkpoint 为 `SCMAT001`：绑定几何导出设置（含 minimum_outer_wall），
  包含原生 solver checkpoint，并对整个 envelope 加 SHA-256。
  改变壳厚、边界复杂度预算、数值请求、trace 预算，或丢失历史 trace，均应拒绝。
  导出策略失配在新的 native trace 目录创建前拒绝。

这仍是完整性检查，不是认证签名：能修改所有记录并重算摘要的同用户恶意写入者
不在威胁模型内。外层 `StageCheckpointStore`、`StageRunSupervisor` 和 `run_source_job`
必须绑定不可变输入、配置、实际扩展、源码与检查器文件。

## 验证范围

小型解析体回归实际检查：先求完整对照，再暂停、编码、解码并继续，最终外表面、
材料、空腔数组字节、完整数值诊断、科学摘要及独立材料封存摘要均相同。
暂停时 artifact/material 必须为空；非法控制和状态不得被当作成功输出。
这是 `SYNTHETIC_REGRESSION`，不是以立方体冒充论文模型。

重建茶壶的跨进程实验使用此前运行的实际配置：

- 输入 `build/teapot-input-20260913-run2`。
- 配置 `build/spin2014-primal-box-20260915/low.json`。
- 配置 SHA-256：`f8c2871b6f2bc85f2fb93c32b7c1aef6400ecb3afe2e284eca4d5c6716d05abb`，
  与此前 `build/spin2014-partition-restore-20260916/g6-run-1/provenance.json` 一致。
- 原输入、姿态、设计域、G6 分辨率和全部数值参数保持不变。
  包括 `FINITE_BOX_PRIMAL_V2` 和 Laplacian 权重 0.0001。
- 两次重复，每次有全程 baseline；封存首阶段后，在下一阶段 trace 尚未封存时
  SIGSTOP/SIGKILL 求解子进程，由仍存活的监督器计费，再用新进程恢复。
- 对照判据预先为 `BINARY_VALIDATION_FAILED`、无几何材料，与旧协议相同。
  恢复一致不等于优化成功。没有更改舍入或放宽约束来生成一个“成功陀螺”。
- 两轮所有尝试处于同一个外层 cgroup：360 秒、24 GiB，启动空间 100 GiB、
  运行最低可用空间 20 GiB、采样磁盘上限 2 GiB。单次最多 120 秒。
  这些是本轮工程验证预算，不冒称论文设置。

`tests/run_source_material_recovery.py` 是自动检查程序，不是任意 OBJ 的用户级调度器。
checkpoint 的 trace-byte 计数只涵盖已封存历史；失败部分文件仍由监督器的全目录
磁盘计账保留和限制。没有通过删掉失败尝试重置同一作业预算。

## 保留的失败记录与修正

证据目录：`build/spin2014-resumed-material-20260916`。

1. `red.log`：新接口不存在时，7 项测试实际失败。
2. 初次测试误写暂停状态名，改为现有 `STAGE_BOUNDARY_PAUSED_NOT_VALIDATED`，
   没有改变生产状态。
3. `material-identity-red.log`：审查指出壳厚/边界预算未绑定，两项实测未拒绝。
   随后增加材料 envelope；`envelope-splice-red.log` 进一步复现头部拼接漏洞，
   再增加整体校验和。
4. `green-second.log`：旧 runner 回归因增量构建目录缺少真实 CMakeCache 而失败。
   在同目录实际配置 CMake 后重测，没有伪造或复制一个缓存来骗过证据收集。
5. `teapot-first`：误选旧 `online_trace_g6.json`，其方向策略与权重不同，
   baseline 在 256 MiB trace 限额停止。全部日志保留。这不是上述 low.json 协议的失败。
   更正后使用新目录 `teapot-correct-protocol`，未增加 trace 配额或修改 low.json。

最终测试以 `verification-final-v2.json`、`final-v2-*.stdout/stderr` 为准；
茶壶结果以 `teapot-correct-protocol/payload/result.json` 和外层 `result.json` 为准。
缺少这些完成记录不能认为对应验证通过。

最终实测 140 项 Python 测试全部通过，包含新增 12 项恢复/封存测试；目标 CTest
1/1 实际执行通过，无 skip/Not Run。茶壶 G6 的二值化归一化约束残差为
`4.647020231279599e-6`，超过不变的 `1e-6` 门槛；独立质量矩误差仅
`2.439749362756517e-16`。这说明当前失败项是二值对象的约束，而不是恢复后质量矩
积分不一致。对应 trace/output hash 与旧 low.json 运行一致，但旧源码证据不升级为当前证据。

真实输入两轮均完成：每轮 2 个已封存阶段边界、5 次尝试（含 baseline 与一次强杀），
最终完整 diagnostics 与不中断对照逐字节相同，两轮摘要均为
`38fdf7bc1962a901c217fe45c2ff7f93a8c9ef402bb3bd833c1871054ccdfdb6`。
整作业实测 290.704 秒，MemoryPeak 404475904 字节（约 385.7 MiB）。
两轮均无 sealed artifact；过程成功不等于优化成功。`evidence-check.json` 是收据、
运行中固定身份、科学结果和源码摘要重新核验后的记录。

## 构建与剩余限制

采用可审计增量构建：复制上一轮 codec core archive，仅重新编译材料 wrapper 和
Python binding，链接上一轮的主 pybind 对象及固定依赖。保留编译命令和冻结 archive
摘要；没有修改旧 archive/扩展。CMake 在当前证据目录真实配置，仅执行选定的
Python/native 测试，不声称已全量重建或全项目测试全绿。

Build Graph 本轮增量解析 16 个文件，无解析错误；图包含 463 文件、6008 节点、
64321 边，语言 bash/C/C++/Python。增量图仍不能保证包含所有未跟踪新文件，
实际接口与测试已直接读源码核对；图节点数不表示论文复现进度。

仍需解决的科学问题是重建茶壶的二值化后约束失配，以及更高分辨率/原始输入的
完整实验。恢复基础设施不修正这些问题。高层 Python 材料桥仍将诊断展开到内存，
不声称整条链路总内存恒定；整个作业的硬内存门仍必须启用。
controller 死亡后的自动接管、任意模型统一调度和新的可视化视频不属于本轮完成声明。
原论文模型与关键参数缺失继续阻塞严格完整复现。主工作区、历史日志不覆盖；无 commit/tag。
