# Spin-It 2014：完整阶段历史的跨进程编码与恢复

## 本轮结论与范围

在 `make-it-balance-paper-sim` 独立工作树实现了生产 C++ checkpoint codec，
把已有 native 阶段边界令牌与 Python `StageCheckpointStore` 持久化载体接通。
新进程可以重新加载全部已完成阶段，核验后从下一阶段继续，不重跑先前 SLQP。

**这是合成回归验证通过的恢复基础，不是新的原论文案例成功结果。**
本轮没有运行 G7、改变容差或论文终止规则，没有新增材料对象、仿真视频或
Spin-It Faster 结果，也不升级论文复现状态。

## 接口与一致性

- `encodeSourceAdaptiveCheckpoint` / `decodeSourceAdaptiveCheckpoint` 使用
  `SCSPV001`、小端字段、明确长度及 SHA-256；元数据上限 64 MiB。
- 完整迭代保存在外部封存 trace 中。`SlqpIterationStore::openSealed` 只读重载，
  验证逐记录校验和、结束标记、记录数及剩余字节预算，不驻留全部迭代。
- 解码从请求的源几何重建分区和质量矩，逐阶段重算历史状态的目标、梯度、约束、
  Jacobian，以及最终 KKT 见证；重放 split/merge，检查前后分区摘要及完整历史摘要。
- 保留 beta、Givens phi、原始全局阶段编号与成功阶段累计 trace 字节。
  新输出目录只在验证完成后创建。检查失败不能作为重新分配全额预算的理由。
- 完成的 SLQP 必须具有有效的原生 LP KKT 见证状态，残差满足原请求容差；
  末条 trace 的完整 KKT 与最终解一致。收敛历史最多包含 `maximumIterations` 条记录。
- 继续令牌依旧不能由调用者手填状态直接构造。编码/解码不负责授予论文完成状态。

## 实际验证

四组配置包含非对称静态体、非零 phi 的 yo-yo，以及两套已有终止规则。
它们是 `SYNTHETIC_REGRESSION`，不是论文模型的替代品。

| 配置 | 每次保存边界数 | 重复次数 |
|---|---:|---:|
| 静态体 / unchanged-partition | 3 | 2 |
| 静态体 / paper-post-transition | 2 | 2 |
| yo-yo / unchanged-partition | 2 | 2 |
| yo-yo / paper-post-transition | 1 | 2 |

每个暂停、往返编码检查及继续动作使用独立 native 进程；每个边界通过真实
`StageCheckpointStore` 发布与加载。全部 16 个边界恢复后的最终科学输出和
逐条完整迭代字节与不中断运行完全一致。

每套矩阵执行 66 次 native 调用，包含 17 项反向检查和 5 项只读重载检查：
截断、附加字节、版本、校验和、计数、长度、预算、目标、非有限 phi、dual witness、
分区、请求身份、重算校验和的历史目标，以及完整重算所有摘要的三项语义伪造。
反向检查通过表示伪造被拒绝，而不是伪造输入获准继续。

独立审查发现最终 KKT 状态、末条 KKT 一致性以及迭代数上界三处漏洞。
三个定向测试先在修复前二进制中实际失败（`forged semantic checkpoint accepted`），
修复后全部通过；最终只读复核没有剩余 Critical / Important 问题。

同一矩阵另以 UBSan 执行通过。**仅 codec 和记录解码相关翻译单元被插桩，
并非整个冻结几何依赖库均经过 sanitizer 检查。**
49 项 Python 存储测试、原 native 阶段恢复回归、源自适应回归也通过；
后者 stdout 与上一轮冻结基线逐字节一致。

## 构建与证据

证据目录：`build/spin2014-checkpoint-codec-20260916`。
保留变更前快照、初始 RED、三项语义 RED、失败尝试、构建命令、正常与 UBSan
矩阵结果、原回归日志以及源码/输出摘要清单。`verification-final.json` 保存最后一次
检查的真实子进程退出码，不接受手填 passed 参数。

本轮仍是可审计增量构建：复制已冻结 core archive，在复制品中更新
`SourceAdaptiveSpinSolver.cpp`、受内部接口影响的 `SourceMaterialBoundary.cpp`、
`SlqpIterationStore.cpp` 和新增 codec，再链接测试。冻结 archive 未改动。
CMake 注册和配置已完成，但本轮未全量重建所有依赖或执行完整发布 CTest。

没有创建 commit/tag，没有改动受保护主工作区。以前绑定旧源码的实验不能作为
当前源码的新鲜证据；历史日志保留不覆写。

## 边界与下一步

1. 当前跨进程集成入口是合成测试驱动，不是通用 OBJ/Python 长实验启动器。
2. 没有跨进程资源监督器；测试只观测正常子进程耗时。尚未累计或强制限制失败尝试、
   整个运行的墙钟、重复分支和全部磁盘占用。
3. 进程在完整阶段边界正常退出；尚未验证阶段中途 SIGKILL、超时、掉电恢复。
4. 解码复核源几何评价、最终 KKT 和转移，不重新演算每个历史优化方向与 merit 决策。
5. SHA-256 是完整性绑定，不是身份认证；要求历史文件私有且不可变，不防御并发恶意
   文件替换。trace 使用绝对路径，迁移整个历史目录需另行设计显式重绑定协议。
6. 当前保证限同环境逐位重启；不把不同编译器/依赖下的数值差异当作可无条件接受。

下一重点：通用生产续跑入口与全局资源监督器，计入失败尝试并执行真实杀死/恢复
反向测试。只有这条链路通过后，才重新安排原长实验；不能用分段重启绕开预算。
