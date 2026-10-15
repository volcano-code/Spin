# Spin-It 2014：细分阈值与二值舍入对照（2026-09-16）

## 结论与范围

完成预先固定的三组阈值、每组两次独立运行。**六次均为
`BINARY_VALIDATION_FAILED`，没有新的封存几何、仿真通过或论文原案例通过。**
减小细分阈值确实消除了此前的粗层舍入贡献，但没有降低总残差。
因此保留原默认配置；本轮没有修改求解器、舍入规则或 `1e-6` 验收门槛。

这是重建茶壶的 `RECONSTRUCTED / IMPLEMENTATION_DEFAULT` 参数敏感性实验，
不是作者原始资产实验，也不是不可行性证书。更细分区上的一次局部解不能证明
精度、目标函数或二值残差随阈值单调改善。

## 固定协议与来源

- 证据目录：`build/spin2014-fractional-sensitivity-20260916/`。
- 预注册文件：`PROTOCOL.md`、`planned.json`；执行器：`campaign.py`。
- 基准：`build/spin2014-primal-box-20260915/low.json`，SHA-256
  `f8c2871b6f2bc85f2fb93c32b7c1aef6400ecb3afe2e284eca4d5c6716d05abb`。
- 复用 `build/spin2014-resumed-material-20260916/python/` 内冻结的实际 native 扩展。
- 固定 G6O6、`FINITE_BOX_PRIMAL_V2`、`gamma_L=0.0001`、论文 split/merge 后终止策略；
  只改变 `fractional_epsilon`，不改变模型、姿态、设计域、壳厚、目标或初始化。
- 每次独立 whole-job cgroup：300 秒、24 GiB、swap=0、4 GiB 采样磁盘上限、
  512 MiB native trace；顺序执行，BLAS/OpenMP 单线程。不是沿用或延长旧超时任务。
- 原始论文 `/home/lin/下载/SpinItPaper.pdf` 第 5.3 节说明了 epsilon 筛选、
  split/merge、最大层级终止与二值舍入；该节没有提供数值 epsilon。
  原文提取保存在 `paper-section5.txt`，不得把本次参数称为作者原值。

## 完整结果

| epsilon | 重复次数 | 自适应 passes | 最终叶节点 | 舍入改变单元 | 舍入前归一化残差 | 舍入后归一化残差 | 固定门槛 |
|---|---:|---:|---:|---:|---:|---:|---|
| 0.05 | 2 | 3 | 252 | 186 | 2.372930e-17 | 4.647020e-6 | 失败 |
| 0.025 | 2 | 3 | 252 | 186 | 2.372930e-17 | 4.647020e-6 | 失败 |
| 0.01 | 2 | 11 | 1498 | 1342 | 8.637923e-18 | 8.644816e-6 | 失败 |

每组两次的完整归档 SHA-256 相同，事件与科学输出比较通过；不同 epsilon 的
输入、trace/output 哈希并不相同，即便前两组的终态数值相同。所有六次均实际运行，
没有超限、跳过、丢弃或替换失败记录。每次约 30.3–51.6 秒，峰值内存最高
413,257,728 bytes（cgroup MemoryPeak）。时间与内存不进入确定性摘要。

独立质量矩检查的归一化误差分别为约 `2.44e-16`、`2.44e-16`、`4.88e-16`。
这支持积分一致性，但不代替最终二值材料满足约束或封闭几何审计。

## 为什么降低阈值没有成功

`stage-audit.json` 对真实终态逐单元重算舍入增量，而不是从期望答案填写结果。
以下是控制最大残差的横向一阶矩约束的**带符号**归一化贡献，不是各层误差范数：

| epsilon | level 4 舍入贡献 | level 6 舍入贡献 |
|---|---:|---:|
| 0.05 / 0.025 | +1.001354e-6 | +3.645666e-6 |
| 0.01 | 0 | -8.644816e-6 |

两个 beta 约为 0.012/0.015 的粗层单元，在 0.05 与 0.025 下未触发细分。
0.01 改变了分区及优化路径，最终所有粗层舍入贡献为零，但 level 6 的 1342 个
改变单元产生更大的净偏移。不能将消除粗层贡献解释为总残差改善。

最后一次优化、split/merge **之前**的记录如下；不是舍入后或最终分区重新优化的目标值：

| epsilon | paper objective | regularization energy | total objective |
|---|---:|---:|---:|
| 0.05 / 0.025 | 1.0780369060526591 | 0.0033767858980411356 | 1.0814136919507003 |
| 0.01 | 1.07523856441194 | 0.010695236112349156 | 1.0859338005242891 |

分区/邻接改变，不能仅比较总目标值就判断同一离散问题优化得更好。
预转换一阶 KKT 证据也不能升级为转换后分区或二值解的最优性证明。

## 验证和审查

`verify_study.py` 重新核对全部六个预注册身份、单变量配置、源码/输入/扩展摘要、
完整 provenance、输出清单、归档 footer、独立舍入分解和三组重复；生成
`verified-analysis.json`。脚本明确不授予几何、物理或论文验收状态。

本轮实际执行检查器测试：

```bash
PYTHONPATH=python:experiments/teapot_spin2014 PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 \
  /home/lin/anaconda3/bin/python -m pytest -q \
  tests/python/test_stage_audit.py tests/python/test_source_trace_archive.py \
  tests/python/test_resolution_rounding_audit.py \
  experiments/teapot_spin2014/test_online_compare.py
```

结果：**56 passed**（`checkers-final.log`）。这不是全项目测试通过声明。
只读审查确认单变量控制与归档检查成立，并发现 `campaign.py` 的潜在异常路径：
任意控制器/证据异常都会令下次重复跳过，却将原因标成资源规则。
本轮六次均正常退出，未触发该路径；原始脚本及其绑定保持不变，不能修改后仍声称
是同一执行器。**该脚本不得未经修正和异常路径测试直接复用于后续批次。**

## 下一重点

不继续事后搜索 epsilon，也不通过改舍入阈值争取通过。
优先完成此前因 300 秒预算未完成的真实 G7O7 对照，使用已验证的分阶段恢复机制，
另行预注册预算、固定配置和两次重复；旧超时记录保留。
该实验尚未在本轮启动，也不保证更细分辨率一定通过。
只有二值约束通过后，才进入独立几何封存、质量矩核验和仿真。

本轮边界：求解器生产代码未改、无提交/tag、无主工作区覆盖，Spin-It Faster 未变更。
