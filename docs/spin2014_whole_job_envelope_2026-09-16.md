# Spin-It 2014：完整作业资源门（2026-09-16）

本轮推进长实验的运行保障，不新增论文原始案例成功声明，不修改科学目标、容差或分辨率。
开发目录仍为 `/data/xvxinyao/moment-guided-spin-v1/code/make-it-balance-paper-sim`。
受保护主工作区未修改；未创建 commit/tag。

## 新增能力

`python/make_it_stand/verification/source_job.py` 提供一次性 `run_source_job()`：

- 整个 payload（求解、重试、checkpoint 发布、材料导出、封存）处于同一个临时
  user-systemd service/cgroup v2。子进程另建 session 也仍在该 cgroup。
- `memory.max` 是整个组的内核内存计账上限，不是把单进程 RLIMIT_AS 冒称组 RSS。
  `memory.swap.max=0`；OOMPolicy=kill。默认 24 GiB、60 分钟。
- 启动前检查至少 100 GiB 可用空间；运行中低于 20 GiB 停止。
  输出目录默认采样上限 80 GiB；这是工程默认值，不是论文参数。
- 完整输入身份检查后、真正启动 payload 前，辅助程序核验自身 cgroup、内核限额、
  boot ID 及原始文件路径对应的哈希。
- Python 3.13 的 BOOTTIME 绝对 timerfd 配合组内守护进程持有的 `cgroup.kill`，
  将截止时间绑定到持久化配置，而不是启动服务后的相对时间。
  systemd RuntimeMaxSec 作为第二道超时限制。控制器被 SIGKILL 后，作业仍到期终止。
- UUID service 只在本次作业内创建与清理，不修改系统配置或其他用户服务。
- 目录独占创建并 fsync 父目录；所有失败日志保留，不允许通过复用目录重置预算。
  `PROCESS_EXITED` 仅表示进程成功退出；科学验收始终为 `NOT_GRANTED`。

环境要求：Linux cgroup v2、可用 user-systemd、受委派且可写的 cgroup.kill、
支持 timerfd 的 Python 3.13、本地可用的 link/flock/fsync 文件系统。
不满足时失败，不自动退化成单进程内存限制。

## 茶壶入口接入

`experiments/teapot_spin2014/supervised.py` 封装原有 `run.py`，保留三种 strategy
及所有既有科学参数。输入 OBJ、manifest、声明来源、Python/C++ 源码、指定扩展
和显式 preload 库绑定哈希。启动脚本一律使用绝对路径。

`--native-module` 必须是当前解释器可识别的扩展文件名；启动器传入
`--expected-native-module`。现有 runner 在创建输出或求解前核对真正加载模块的
`__file__`，避免指定 A、实际加载 B。旧直接调用不强制新增参数，保持 legacy 兼容。

现有 `build/teapot-input-20260913-run2` 已完成输入/命令预检：292 项身份绑定。
本轮没有重新运行茶壶求解；预检不证明扩展的编译来源完整，也不是材料封存通过。
该输入永久属于 `RECONSTRUCTED_EXISTING_TEAPOT_PROTOCOL_NOT_PAPER_ORIGINAL`。

## 验证与审查记录

证据目录：`build/spin2014-job-envelope-20260916`。
以 `verification-final-v2.json` 及 `final-v2-*.stdout/stderr` 为本轮最终检查记录。

最终实测：102 项 Python 测试全部通过（整作业 17、茶壶入口 8、单次监督 28、
存储 49）；指定 CTest 1/1 实际执行通过，未跳过。源码及结果清单重新校验通过。
完整 native 中断恢复作业耗时 16.514 秒，内核 MemoryPeak 为 105132032 字节；
该次配置为 60 秒、1 GiB、小型合成输入，不代表 24 GiB 原始模型实验性能。

覆盖的检查包括真实内核内存限额、两个子进程合计 OOM、求解后封存仍计时、
控制器死亡、服务延迟启动、独立 session 子进程清理、磁盘采样超限、低可用空间、
不可重用作业目录、执行文件软链接改指向、payload PYTHONPATH 污染辅助程序导入、
原生扩展加载身份错配和相对入口路径。

独立只读审查促成修复：辅助程序使用 Python `-I` 并固定本地验证模块搜索路径；
校验原始路径而非只校验第一次解析得到的目标；绝对截止看门狗；父目录持久化；
真正加载的扩展与指定文件一致。未据此宣称所有潜在缺陷均已排除。

首次最终检查发现 CTest 在未注册测试时仍返回 0，**该项证据作废**。
保留 `verification-final.json` 和原日志；新验证必须使用 `--no-tests=error`。
CMake 现在拒绝“启用 whole-job 测试但缺少其注册依赖”的配置。
`registration-guard.log` 保存该配置被拒绝的实测证据。

整个既有 native 强杀/恢复驱动也在新作业门内实际执行：两个合成配置各重复两次，
共 4 次阶段中途 SIGKILL、18 次含失败的恢复链尝试及 4 次不中断对照。
科学输出及完整最终轨迹与各自不中断对照逐字节一致。
这些仍是 `SYNTHETIC_REGRESSION`，不能增加论文案例完成数。

## 限制与下一重点

1. 磁盘仍是采样监测与发布预留，不是硬配额；快速写入可以瞬时超量。
2. 控制器的预检/清理/最终审计不在 payload cgroup 的内存计账内。
   其耗时纳入截止准入和最终状态检查，但最后收据的 fsync 开销不保证在截止前完成。
   独立看门狗约束科学 payload 的生命周期，不承诺操作系统调度零延迟。
3. 控制器死亡后没有完成收据；记录保持未完成，不自动接管或重开一个新预算。
   本轮恢复成功指 payload 内监督器仍存活时恢复被强杀的求解阶段。
4. 要求同用户合作、不恶意修改配置/源码或刻意迁移 cgroup；不是恶意代码安全沙箱。
5. C++ 未修改；复用上一轮 codec 二进制，核验 src/include 与旧源码清单一致。
   本轮只配置并执行目标 Python CTest，没有全量重编译或全项目测试通过声明。
6. 茶壶入口已接入整个作业门，但任意原始网格的 native 阶段 checkpoint、材料封存、
   Python 通用恢复入口仍需贯通。下一步应先以一个固定重建输入验证这一完整桥接，
   再启动受预算约束的长实验。作者资产/关键参数缺失的原论文案例继续保持阻塞。

源码与结果 SHA-256 清单不包含互相引用的报告循环；旧输出仍绑定旧源码，不做升级。
本轮完成状态仅限工程模块，不是 `VERIFIED_COMPUTATIONAL_REPRODUCTION`。
