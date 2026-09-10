# Task5i：数值边界修正与挖空链路复验

日期：2026-09-10。实现提交：`17f113ac84fcb5fa058f101b186f2accd3664de5`。
开发仍位于 `/data/xvxinyao/make-it-balance-paper-sim`，未覆盖主项目未提交修改。

## 本轮实现

修正原始表面八叉树的精度准入：不再仅因第一个最细单元正常，就允许后续可能重合的网格面进入细分。真实原生输入复现了旧代码先通过配置检查、随后报坐标重合的行为。现在提前明确拒绝不受支持的精度，不改变模型、不降分辨率，也不判定模型物理不可行。

此边界是保守的：某些实际上可用的等ULP对齐输入也会被拒绝。普通输入的来源、配置、分区及下游细网格输出保持不变。它不是新的挖空目标函数或完整BBW求解器。

## 已执行验证

- 修正前原生回归真实失败，修正后通过。
- 主配置与关闭Faster配置的相关集成测试分别 **8/8通过**。
- 最终全部目标构建成功，CTest **201/201通过，474.14秒**；测试编号无重复、无遗漏，无NotRun/Skipped记录。
- 独立代码审查通过，没有发现本轮阻塞性问题。
- 独立有理数舍入检查覆盖1,400组数值配置；获准配置没有发现重合区间。有限扫描不代替全域数学证明。
- 两次实际细网格导出逐字节一致，且与上轮普通输出一致：1,728个单元、2,197个节点、6,084条边、完整218个控制点。

## 新生成的小型挖空结果

使用原有固定解析方体回归输入，没有重新选择姿态、设计域或初值。两次当前源码下的实际优化、材料生成、独立封存、重新加载与积分均通过：

| 项目 | 两次结果 |
|---|---|
| 分类 | `SYNTHETIC_REGRESSION`，不是论文原始模型 |
| 数值状态 | `BINARY_MOMENT_VALIDATED` |
| 记录新鲜度 | `CURRENT` |
| 求解轨迹 | 每次39条记录 |
| 挖空单元 | 每次72个cuboid |
| 实际合并空腔 | 每次1个 |
| 最终材料质量 | `503/64 = 7.859375`，原质量8 |
| 最终质心 | `(0,0,2003/2012)` |
| 外表面 | 顶点与面字节不变 |
| 两次输入、轨迹、科学结果 | 分别逐字节一致 |

结果目录（均在开发项目内）：

- `build/task-5i-controller-W6WCqr/material-run2/`
- `build/task-5i-controller-W6WCqr/material-run3/`

每份含 `summary.json`、`inputs.json`、`native-run.json`、`scientific.json`、`evidence.json` 和 `sealed/` 数字材料包。封存几何可用于后续查看器加载验证；**本轮没有重新打开或验证真实桌面OpenGL，也没有进行陀螺旋转仿真**。

实现提交后再次独立加载两份结果，均仍为 `CURRENT`；输入、轨迹、几何及质量矩检查保持通过。

首次手动启动漏用了现有CTest配置的系统C++运行库，产生导入失败，保留在 `material-run1/failure.json`。随后与CTest统一 `LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6`；未安装依赖或修改算法。失败尝试没有计入上述两次成功结果。

## 仍未完成

整体仍为 **`PARTIAL_BLOCKED`**。几何封存通过不等于目标稳定性通过；本例仍明确记录 `solver_target_feasibility=UNKNOWN`、`build_source_attestation=NOT_VERIFIED`、`physical_validation=NOT_IN_SCOPE`。

下一核心工作仍是BBW离散能量、质量矩阵、完整面约束和求解器验证，再接联合变形。本轮核查还明确了四面体质量矩阵规则，以及三角形边界权重经三线性插值后的误差问题，详见 [BBW契约核查](spin2014_bbw_contract_audit.md)。Faster剩余算法、CAD和两篇论文原始案例/关键参数缺口没有因此关闭。

完整源码、二进制、日志摘要与独立检查记录见 `.superpowers/sdd/2026-09-05-spin-it-strict-computational-v2/task-5i-controller-verification.md`；任务独立代码审查单列记录，不能由测试数量替代。
