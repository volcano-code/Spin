# Faster 无答案初始化扩展批次（2026-09-23）

## 范围与证据等级

这是固定重建立方体问题上的合成方法验证，不是 Spin-It Faster 作者案例，也不计入论文案例完成数。原始网格、设计域和控制参数在运行前固定；求解器未接收用户提供的二次曲面、旧终态或已知解。全部失败、重启和事件均保留。

预注册：
[`faster-default-start-preregistration.json`](../experiments/spin_it_reproduction/multi_case_2026-09-23/faster-default-start-preregistration.json)，SHA-256 `de3869206f1c5253b45d841a860eaad03d352e505d6ae490dc008dbf9870cd53`。

runner SHA-256：`bb89451a05e3d19157282926fda9521698f4112156c1793641628cb574418e60`。实际加载扩展来自
`build/strict-algorithms-exact-current-20260922/python/_make_it_stand.cpython-313-x86_64-linux-gnu.so`，SHA-256 `2ec0efd0c508b2cc15af17b04547185e11dee20bc2fc0408f26ee51eead21949`。运行时预加载与该构建绑定的 GNU `libstdc++`，以避免 Conda 较旧运行库先行加载。

新种子 20–24 每个在独立 Python 进程中运行两次。已观察过的 0–2 exploratory 和此前预注册的 10–12 结果不被混入本批次，也不用于替换种子。

## 结果

| Seed | 两次结果 | 重启内函数评估 | 归一化残差 | 独立原始矩残差 | 输出 SHA-256 |
|---:|---|---|---:|---:|---|
| 20 | 两次 `RESOURCE_LIMIT_REACHED`；封口认证失败 | `[60, 1, 60]` | `5.5835486e-1` | 不适用（无最终几何） | `b30b45e1a9dcb98fde2fa58925b16613529c763754ea77099526d31b5318d85b` |
| 21 | 两次 `CONVERGED`，最终几何独立验证通过 | `[34]` | `5.4745851e-16` | `2.7870254e-15` | `cf4f6934faa3145db61d1e9e22869f60fd8c6090eaa91f114460bab3158a0786` |
| 22 | 两次 `CONVERGED`，最终几何独立验证通过 | `[1, 17]` | `8.1153746e-16` | `4.1547747e-15` | `55fc8df147365d231572602e85ad976b97a8b3c1bbe2744744de8d4f3d80571f` |
| 23 | 两次 `CONVERGED`，最终几何独立验证通过 | `[18]` | `4.3362132e-14` | `1.6985349e-13` | `6009fc25294d0dbba72d85a6c77a073e3e0be067e978de0586defad179e38c8e` |
| 24 | 两次 `CONVERGED`，最终几何独立验证通过 | `[21]` | `4.9800492e-16` | `8.5958073e-15` | `5a0a219cdb94c032492fa75d3d8c29d824d9de5fa626d9b83c9d88e4b555c9a0` |

所有 5 组重复对在 `status`、最终几何验证、输出／事件／轨迹摘要、完整事件、各次重启评估数上完全一致。Seed 20 两次均报告 `truncation cap is not certified outside design domain`；这是当前求解尝试失败，不是不可行证书，也没有被剔除。

单次结果文件位于 `build/faster-default-start-seeds-20-24-20260923/seed-{20..24}-run{1,2}`。独立读取脚本已核验 10 份 JSON 的 seed、扩展哈希、失败诊断及重复字段；核验通过。运行只证明当前同一离散 Boolean—质量矩—Jacobian—Newton 路径在该合成立方体例子上，对部分预登记默认初值可重复求解，并暴露了另一初值下的封口认证阻塞。

证据限制：这些单次 JSON 是既有诊断 runner 直接生成的记录，绑定了 Git commit、runner 与 native 扩展摘要，以及运行时观察到的 worktree-status 摘要；它们没有在每次运行开始前封存 `PaperRunEvidenceV2` 的完整有效源码差异和依赖快照。因此，本批次不进入严格发布验收器的 `CURRENT` 证据集合；不能用报告中的重复一致性替代更强的源码证据链。

## 结论边界

这不是广泛初始化鲁棒性证明，不是完整曲面族求解器验证，也不是作者模型或 Table 1 实验。Seed 20 的失败不表示问题不可行；它指向截断封口策略在该重启序列上的认证限制。整体复现状态仍为 `PARTIAL_BLOCKED`。下一项算法缺口仍是秩亏曲面对应的真实低维优化路径，以及 Cherchi／CGAL 的离散材料几何对照。
