# Spin-It 联合形变与双材料定向审计（2026-09-23）

本记录是严格计算复现中的诊断证据，不是原论文案例的完成声明。所有输入都仍是
`SYNTHETIC_REGRESSION` 或重建资产；项目状态保持 `PARTIAL_BLOCKED`。

## 联合形变：primitive 0 的包含关系失败

在冻结的 98-handle source-cage 配置上，新增了 primitive-subset 回放见证。它保存：

- 发生拒绝的 callback state 和 state SHA-256；
- 同一状态重新调用 source evaluator 后的外表面和 primitive 0 OBJ；
- primitive 原始点索引和全部 98 列 BBW 插值权重；
- guard 的输入／输出摘要，以及重新求值的输出摘要。

两次独立运行的可比较科学摘要完全相同：

| 运行 | 科学 SHA-256 | 状态 | callback / geometry stage |
| --- | --- | --- | --- |
| `source-cage-joint-primitive-witness-20260923-d` | `99613001d49a22f4dc792336d113b401c26641993ba9cf702f26dea2d177afd7` | `RESTORATION_STALLED`, `NOT_ADMITTED` | `1 / 14` |
| `source-cage-joint-primitive-witness-20260923-e` | `99613001d49a22f4dc792336d113b401c26641993ba9cf702f26dea2d177afd7` | `RESTORATION_STALLED`, `NOT_ADMITTED` | `1 / 14` |

两份见证的外表面和 primitive OBJ 也逐字节一致。每份见证中
`replayed_evaluation_output_sha256 == geometry_output_sha256`，因此 guard、导出器和
诊断回放使用的是同一 source evaluation，而不是不同的点索引、坐标或插值路径。

`primitive_subset_witnesses.jsonl` 中的五个 BBW 诊断系数可出现最大
`9.56e-10` 的重建级差异；它们不改变 callback 轨迹、外表面或 primitive 网格。该附件
仍封存在原始证据包中，但已从用于跨进程确定性比较的**科学输出摘要**中隔离，避免将
诊断附件的浮点噪声误报为求解器不确定性。

结论：当前配置遭遇的是经 exact-material guard 重现的实际包含关系边界，尚无“求解器、
检查器和导出器路径不一致”的实现缺陷证据。它也不是连续设计域不可行证明；状态为
`RESTORATION_STALLED`，不应据此修改权重、关闭包含检查或宣称不可行。

## 双材料分区与质量矩

双材料导出 manifest 现同时封存输入 material/carveable mask、固定壳／移除／两种材料的
label partition、target frame、binary constraints 和约束无穷范数。独立检查器执行：

1. 从输入 mask 验证固定壳没有被删除、dense insert 不越出 carveable domain，且 labels
   是互斥的 voxel partition；
2. 对每一标签的体素 cuboid 解析积分；
3. 对实际导出 OBJ 独立积分并检查封闭二流形、有限顶点和域边界；
4. 以真实密度组合十项质量矩，并独立重算当前 target-frame 下的 Spin-It 约束。

运行 `dual-material-partition-audit-20260923-b` 的结果：

| 案例 | 低密度 / 高密度 cells | 固定壳 / 已移除 carveable | 最大几何—标签矩误差 | 约束无穷范数 |
| --- | ---: | ---: | ---: | ---: |
| `dense_regression` | 56 / 4 | 56 / 4 | `8.88e-15` | `0` |
| `sumone_reconstructed` | 192 / 0 | 192 / 196 | `1.61e-15` | `5.74e-18` |

`dense_regression` 因而仅获得合成双材料方法闭环证据。重建 Sumone 仍是
`NO_RETAINED_DENSE_INSERT`：高密度内件为空，不能被标为双材料优化成功，更不能计入作者
Sumone 论文案例。

## 后续

## Faster：全秩严格映射

全秩严格路径现将实际生效的导数方法写入 `FullRankMomentEvaluation`，并由回归门验证为
`discrete-fullrank-quadric-boolean-transport-v2`。这意味着该路径的 Boolean 值、十项质量矩、
Newton Jacobian、残差和线搜索来自同一个离散材料映射；旧 generic solver 的连续导数默认值
保留为兼容接口，严格 runner 不调用它。

两次完整合成方法运行的摘要均为
`01f01fee08be871e14a22e1f6164379b2ba8f799546ad4f0a528828294cab680`。它们覆盖 adapter、
多步长链式差分、受保护失败路径、Newton budget、独立体积分、五种 mode dispatch、off-target
yo-yo 的 Boolean topology relinearization 和进程内回放。off-target yo-yo 收敛时的 raw
residual 是 `3.21e-15`，并通过了独立 final-material/KKT 审计。

该结果仍只说明合成全秩 Boolean→moment→Newton 方法链可重复，不能替代秩亏曲面、Cherchi
backend 对照、二维／浮力／CAD 链路或任一作者原始 Faster 案例。

本批次还发现平面族已有 4 参数（不是伪造十参数）离散值／导数实现和多步长差分夹具，但它
此前未进入 CTest。已将 `PlaneMomentEvaluator.cpp` 和
`spin-faster-plane-family` 注册到 Faster 构建；实际运行结果为 **1/1 通过**。这只关闭了
平面家族的测试注册缺口，柱面、抛物面及其各自低维 Newton adapter 仍未完成。

联合形变该配置不再盲目重跑；下一主线继续完成 Faster 的秩亏曲面与独立 Boolean 后端，
联合形变只在有新的、可测试的几何缺陷假设时重新进入修复批次。
