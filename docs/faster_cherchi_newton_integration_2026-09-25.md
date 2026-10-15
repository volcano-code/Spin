# Faster: Cherchi Boolean 到 Eq.14 Newton 集成诊断（2026-09-25）

## 结论

严格离散求解入口现可调用外部材料几何评估器。一次限预算合成诊断在 Newton 每次函数评估中实际调用固定 Cherchi Boolean，再从该次 Cherchi 输出恢复来源运动、计算十项原始质量矩及其 Jacobian，并将几何／拓扑身份写入 Newton 事件。

本次**未收敛**：初期 3 次评估 smoke run 验证了拓扑切换接受/拒绝；随后当前二进制绑定的完整单次 60 次预算运行推进到 Cherchi 输出边界不满足闭合绕向的 Jacobian 前置条件，最终到达评估上限。结果是算法集成与阻塞证据，不是成功优化、不可行证明或论文案例。

## 输入和协议

- 类别：`SYNTHETIC_REGRESSION`。
- 输入：既有固定设计域盒体、固定 ellipsoid 初始二次曲面、`YO_YO` 目标；未读取或注入答案。
- 初期 smoke run：`DUAL_MONOMIAL_RAW_V1`，1 个重启、3 次函数评估。主诊断随后使用同一初值与目标、1 个重启、最多 60 次函数评估。
- 每个回调：重建 quadric source → 调 Cherchi `intersection/subtraction` → 读取精确来源 sidecar → 计算该输出网格质量矩及固定拓扑导数。
- 终止状态：`RESOURCE_LIMIT_REACHED`；诊断为最后一次拓扑转换试探不降低残差并达到评估预算。该有限预算是本方法探针配置，不表示主机资源耗尽。
- 最终对象：未通过收敛门，未运行最终独立材料约束验收；求解器明确返回 `final_geometry_verified=false`。

## 实际结果

有效的限预算拓扑烟测在 [r3 证据目录](/data/xvxinyao/moment-guided-spin-v1/code/make-it-balance-paper-sim/build/strict-faster-cherchi-newton-integration-20260925-r3/)；当前二进制的完整预算诊断在 [r11 证据目录](/data/xvxinyao/moment-guided-spin-v1/code/make-it-balance-paper-sim/build/strict-faster-cherchi-newton-integration-20260925-r11/)。

| 评估 | Cherchi 顶点数 | 面数 | 几何摘要前缀 | 拓扑摘要前缀 | 结果 |
|---:|---:|---:|---|---|---|
| 0 | 170 | 332 | `e80888a2278e` | `f931b3d1260e` | 初始几何评估通过 |
| 1 | 160 | 320 | `6d06e80b53a7` | `1fd7c411373a` | 拓扑改变；独立残差下降，接受并在下一轮重新线性化 |
| 2 | 194 | 388 | `7b4f193c4cd3` | `453763c26952` | 拓扑改变；残差不降，拒绝 |

求解事件终结于 `maximum_evaluations_reached`。拓扑变化后没有沿用旧 Jacobian：新一轮使用重新计算的 map；最后一次试探仍因残差条件被拒绝。

## 实现与验证

- 新增严格 Python/native evaluator bridge：要求回调同次几何给出 10 项目标坐标 raw moments、10×10 系数 Jacobian、非空几何与拓扑 hash、导数有效状态。
- Bridge 返回数值轨迹、事件、重启计数和 hash，并固定声明 `numerical_newton_only_external_geometry_not_sealed`；它不能产生封存对象或授予最终几何通过。
- 反向测试验证 `derivative_defined=false` 时无法通过，且事件保留几何／拓扑摘要。
- native extension 构建通过；本轮最终相关 CTest 通过 4/4（quadric triangulator、Cherchi 来源、Cherchi motion、完整 Faster Python API；总计 70.83 秒）。新增 fail-closed 测试包含在 API CTest 内。它们不等同项目全量测试。

## 保留的失败尝试

- `...-r1`：Newton 调用结束后报告序列化拒绝 `inf`；目录保留。后来修正为显式字符串表示非有限诊断值。
- `...-r2`：Cherchi 子进程收到相对路径，无法打开设计域；目录及错误结果保留。后来将输出根路径规范为绝对路径。
- `...-r3`：上述修复后的有效实际 Boolean→Newton 运行。
- `...-r4`：扩展至 60 次求解器评估，但回调仍调用 CGAL 来逐次判定符号；第 60 次遇到 CGAL 截断封口无法证明位于设计域外。该运行保留，且不作为“同一几何映射”证据。
- `...-r5`：移除迭代中的 CGAL 求值；在固定正定椭球支上以解析中心符号选择 intersection/subtraction，仅保留一次初始 CGAL 对照。60 次 Newton 评估中 Cherchi Boolean 实际产出 32 次；其余回调在求解器重建四次曲面截断边界时失败。末次诊断为截断轮廓某顶点出现两个出边。Newton 按容差拒绝拓扑不降残差试探，随后遇到该几何失败并达到评估上限。
- `...-r6`：同 r5 的数值调用在序列化新增失败试探系数时遇到 Eigen view 不能直接返回 pybind 的绑定错误；源代码修正为物化 Eigen 向量后通过反向测试。目录保留。
- `...-r7`：同配置取证成功归档。r5 的六个 `two outgoing edges` 失败试探系数已随求解事件完整保存。调用层级 2 下，固定曲面经 4 次半径扩展仍不能形成合法截断环。
- `...-r8`：层级 3（项目默认）完整预算运行。失败由 runner 过严的“跨试探源三角索引完全相同”门引起；检查器确认该门不必要，因为每个 q 都重建了当前源网格、Boolean 与匹配 Jacobian，Newton 自身按材料拓扑身份处理重线性化。该失败促成移除 runner 级跨迭代索引相等限制。
- `...-r9`：撤销索引限制后的重跑暴露了临时文件编号复用（失败回调不进入成功列表）；目录保留。编号已改用独立 evaluator-call 计数。
- `...-r10`：层级 3、固定 60 次预算、无 CGAL 迭代混入的有效 Newton 集成。60 次 callback 中有 29 次生成并积分 Cherchi 材料网格；同一拓扑转换附近，Cherchi 输出在某次试探上未能形成闭合一致绕向边界，因此该回调的 Reynolds Jacobian 被 fail-closed 拒绝。该 q 及完整失败事件保存在 result JSON；这不是不可行证明。
- `...-r11`：当前二进制绑定的最终重复；与 r10 同为 60 次求解器评估、29 次 Cherchi 材料网格输出。第 22 次评估（迭代 5）的 Cherchi 边界不闭合且绕向不一致，导数被拒绝，最终预算耗尽。二进制 SHA 已与输入清单核对一致。

不得覆盖这些诊断目录。有效输入/输出是受控合成夹具，不具有作者原始案例身份。

## 剩余缺口

1. 该桥接尚未证明存在非平凡收敛解；当前轨迹暴露了 Newton 附近的拓扑切换和预算敏感性。
2. Cherchi 对不同二次曲面与 CGAL 的精确实体差异仍独立存在；不能将该 Boolean 后端宣称为论文实现等价。
3. `LNC` 来源运动仍不可用；不能跨该构造使用失效 Jacobian。
4. 成功数值解还需对最终 Boolean 材料重新积分并独立检查约束；本次未到此阶段。
5. Spin-It 的联合形变完整材料闭环和本计划其他 Faster 应用算法仍待验收；原论文资产、参数、正式案例矩阵、大规模统计、CAD/浮力/二维集成及最终制造与仿真门仍未完成。

## 2026-09-26 status addendum

The historical statement above that LNC source motion was unavailable is
superseded for a fixed-parameter, fixed-arrangement construction: the pinned
Cherchi source stores the interpolation parameter as a constant and evaluates
the point as `(1-t)P+tQ`. The current provenance parser retains that parameter,
the local motion differential applies `(1-t)dP+t dQ`, and its tests pass at
three finite-difference steps. No LNC point occurs in the inspected retained
Newton fixture outputs, so this update is method-unit coverage, not proof of a
Cherchi Boolean-to-Newton LNC trajectory. The earlier non-closed boundary,
full-run nonconvergence, and Cherchi/CGAL solid-discrepancy findings remain
unchanged.

对 r11 的解释：`RESOURCE_LIMIT_REACHED` 是本次 60 次 evaluator 调用上限触发的求解器状态，不是主机内存/时间上限，也不是数学不可行证书。最终未收敛的直接障碍为 Cherchi 输出边界不满足闭合、一致绕向的导数前置条件；更早的 level-2 native truncation-boundary 失败是独立低分辨率构造缺口。r10 与当前 extension hash 不符，保留但视为旧二进制证据，不用于当前结论。不要调权重或静默切换 Boolean 后端。

总体复现状态仍为 `PARTIAL_BLOCKED`。
