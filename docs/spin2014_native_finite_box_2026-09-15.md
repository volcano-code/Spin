# Spin-It 2014：原生有限盒 QP 接入与两次完整入口验证

## 结论

已接入独立可选的 `FINITE_BOX_QP_V1`，完成小型解析、非线性集成、
默认兼容、存储及独立证据复核测试，并以原低 γ 配置从输入完整启动两次。
**两个实验均为 `RELAXATION_UNKNOWN`，不是优化成功，更不是论文完整复现。**

两次均在第一个自适应阶段停止，44 个填充变量、1 个 Givens 角。
不能将此处的残差与旧策略第三阶段的 94 个填充变量直接作改进对比。
没有进入二值化，没有新封存材料、OpenGL 成功结果或物理仿真视频。

项目目录：`/data/xvxinyao/moment-guided-spin-v1/code/make-it-balance-paper-sim`。
证据目录：`build/spin2014-finite-box-native-20260915`。

## 接入内容与边界

局部问题为 `min g'd + .5 d'Bd`，约束 `Ad+c=0`、
`-beta <= d_beta <= 1-beta`，角变量不受盒约束。
保留实际 L-BFGS 历史和原有活动身份变化重置规则。
使用矩阵无显式展开的 BFGS 作用和有资源上限的 Schur 工作集求解。
这是独立数值实现策略，不宣称为作者公开的原始数值默认实现。

- 旧默认及 `METRIC_TANGENT_RESCUE_V1` 保持独立。
- 新策略每次用完整有限盒 QP，而不只在零步长时救援。
- 失败的 QP 停止为 `NUMERICAL_UNKNOWN`，不使用该方向，也不证明不可行。
- 原有非线性 merit 回溯及显式一阶 fallback 仍可使用并记录；本次实际为零次。
- 共享历史二进制记录载体，以 `BOX_QP_*` 状态与 `finite_box` JSON 区分；
  上下界 face ID 为 `2*i` / `2*i+1`。全部方向、乘子、工作集、预算、标志入 hash。
- Python 从实际已应用的历史独立构造直接 BFGS Hessian，核查等式、驻点、
  原始有限盒、乘子符号、互补性、实际接受试探和下一状态。不只相信 checked 字段。

本次 add/release 工作集过程有明确局限：中间工作集可能病态或无法产生盒内候选。
有界预算及循环/秩检查会拒绝该次求解，不保证所有可行凸 QP 都被此实现求出。

## 新鲜验证

- `python-tests-final-v2.log`：96 项专项 Python 测试通过。
- `ctest-final-v2.log`：7 个定向 CTest 测试组通过。
- 上述测试存在重叠，不相加为 103 项独立测试，也不代表全仓库发布门。
- 六个四维凸非线性问题达到 KKT 门，目标与独立枚举全部 `3^4` 个最终盒面一致。
- 完整内存/磁盘诊断、迭代读回、科学 hash 一致。
- 旧冻结二进制与新二进制的小型默认/tangent 兼容摘要相同：
  `frozen-old-compat.log` 与 `new-compat.log`。
- 独立只读审查修复了精确边界零步长、上界 sub-ULP 漂移、清空接受标志
  绕过检查、大乘子与极小 slack 的互补性，以及浮点消去/非法容差问题。
  均先观察回归测试失败，再修复。详见 `REVIEW.md`。

旧二进制仍未覆盖，SHA-256 为
`7f84620395ec7a6ce1c2cba348a37a58ad71c774e24a06e3d823c8c51efe2146`。
新二进制 SHA-256 为
`4f14dd1c8b5229c7c25845cb2b946c65186f31f8e26c537404590ceb7a75c135`。

## 两次实验

预注册协议见 `PROTOCOL.md`。输入、γ=0.0001、G6/O6 上限、初始化、
960 次迭代/阶段、12 阶段、12000 叶节点和全部容差不变。
配置与旧 low.json 的机械比较确认：仅 direction policy 及对应预算键名不同。
顺序运行，各限制 24 GiB 地址空间与 300 秒，OMP/BLAS 单线程。
两次均正常退出并完整保存失败状态，不是超时或资源限额。

| 指标 | run-1 | run-2 |
|---|---:|---:|
| 运行器总耗时 | 28.1415 s | 28.0997 s |
| 峰值 RSS | 336220 KiB | 336024 KiB |
| 自适应阶段 | 0 | 0 |
| 状态记录 | 109 | 109 |
| QP 调用 / 通过 / 接受 | 109 / 108 / 108 | 109 / 108 / 108 |
| 接受步长 | 全部 1 | 全部 1 |
| fallback / restoration | 0 / 0 | 0 / 0 |
| 终态驻点残差 | 4.2619516e-4 | 4.2619516e-4 |
| 终态等式残差 | 6.08185e-14 | 6.08185e-14 |
| 状态 | RELAXATION_UNKNOWN | RELAXATION_UNKNOWN |

驻点门仍是 `1e-7`。第 109 次 QP（iteration=108）工作集检查报
`RANK_DEFICIENT`，没有接受该方向。二值化字段为 null/未执行。

每份完整归档为 111 个事件（header、pass、109 个状态）。
两次独立 UUID，逐条科学事件及完整归档 SHA-256 完全一致：
`59142750f29aa7b7b745ad49798dc76e1d6909dc2a1951c5789d1b7fb20edbf4`。
每次全部 108 个接受更新经独立复核，下一状态误差为 0，
完整盒 QP 最大独立 KKT 误差为 `3.598358590262141e-15`。
这不验证 QP 的每个中间工作集选择、非线性全局最优或最终几何。

`verification.json` 自动绑定并检查当前源码观察集合、未提交差异、输入、
配置、新二进制、构建设置、完整 trace、输出和检查器。
绑定沿用明确的 `CURRENT_SOURCE_AND_BUILD_FILES_OBSERVED_NOT_COMPILATION_ATTESTATION`，
不冒称历史二进制的编译来源已经得到外部认证。任一有效源码变化后应重新标记
STALE_RESULT，不能改写旧 provenance。报告/输出不进入循环源码 hash。

## 失败点的只读诊断与下一重点

`rank_diagnostic.py` 重建失败时真实的一对 BFGS 历史以及 39 个工作面。
45 个变量、5 条等式、39 个工作面，共 44 行。
以独立 SVD 的相对 `1e-12` 门检查：

- 约束矩阵 C 的秩为 44，最小奇异值约 `3.73e-7`。
- 未缩放 Schur 系统的数值秩为 43，最大/最小奇异值约
  `7917.20 / 1.58e-10`。
- 对角缩放后 Schur 数值秩为 44。

因此不能把原生秩检查失败解释为原设计问题不可行。
但也不能仅凭缩放后满秩就提交方向：一次诊断性缩放求解的等式残差为
`5.98e-8`，超过 `1e-9`；完整盒违反量约 `3958.77`。这个候选没有被应用。
工作集线性方程可解，不等于完整受界 QP 已解。

下一重点是将失败点固定成回归夹具，改进**保持原始可行性的工作集更新**，
结合适当的行尺度处理或数值更稳健的分解，并继续独立核查完整 KKT。
不能只降低 rank threshold、扩大容差或把无界候选截断后声称成功。
新过程通过该失败夹具与既有默认回归后，再冻结新二进制并从同一输入重复两次。
此前局部 94 变量 probe 的八个成功方向不等于从原始输入可完整收敛。

## 保全与复查入口

`pre-change.tar.gz` 保留本轮修改前的目标文件，
`implemented-source.tar.gz` 保留本轮改动文件及新增测试。
没有更改主工作区，没有创建 commit/tag，也没有声称工作树干净。

环境需按既有设置使用 Conda Python 与系统 libstdc++，例如：

```bash
cd /data/xvxinyao/moment-guided-spin-v1/code/make-it-balance-paper-sim
env LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6 OPENBLAS_NUM_THREADS=1 PYTHONPATH=python \
  /home/lin/anaconda3/bin/python build/spin2014-finite-box-native-20260915/verify_runs.py
```

算法接入已完成，专项检查已通过，两次工程实验均诚实记录未收敛；
原论文实验、二值几何与物理验证仍未因此升级。
