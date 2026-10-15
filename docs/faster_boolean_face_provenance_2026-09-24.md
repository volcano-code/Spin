# Cherchi Boolean 输出来源审计（2026-09-24）

## 本次新增

固定版本的 Cherchi 参考探针在 `--dump-exact-rational` 模式下额外写出
`<output.obj>.faces.tsv`。每一行按规范化 OBJ 面顺序给出生成该面的输入**网格对象**位掩码。原 OBJ 字节保持不变；常规探针模式不产生 sidecar。

Python 诊断器 `verification.cherchi_face_provenance.recover_face_sources` 使用精确有理数，对每个输出三角形枚举其完整落在的输入三角面，并将候选输入对象与 Cherchi 对象标签核对。输出状态为：

- `UNIQUE`：恰有一个输入三角面候选，且对象标签一致；
- `AMBIGUOUS`：存在多个候选（例如重合输入面），不能据此选择唯一导数来源；
- `UNMATCHED`：没有任何输入三角面包含输出三角形；
- `LABEL_MISMATCH`：几何候选与 Cherchi 对象标签不相符。

坐标匹配通过 sidecar 中的 OBJ 坐标键完成，因为 OBJ 顶点会被规范化重排；sidecar 内部顶点序号不可直接当作 OBJ 顶点序号。

## 当前证据

- 精确几何盒体交集：28/28 输出三角面均唯一映射到一个原始输入三角面，且 Cherchi 对象标签一致；两个输入对象都有输出边界面。
- 新增两项单元测试覆盖重合面歧义和无来源面的 fail-closed 行为。
- Cherchi 合成交集／差集质量矩与重复输出回归通过；本次 sidecar 不改变 OBJ SHA-256。
- Faster 五类曲面合成 Boolean runner 已接入诊断：先保留两次普通输出作为正式交叉比较，再单独运行 exact-sidecar 路径；仅当诊断 OBJ 哈希与普通输出一致时才生成候选源面记录，否则报告 `DIAGNOSTIC_OUTPUT_MISMATCH`。
- 该 provenance runner 输出已计入运行证据，并由注册的 CTest `spin-faster-cherchi-face-provenance` 检查 Python 来源恢复单元测试。

## 明确限制

这只是**合成夹具上的几何来源恢复诊断**，不是通用的布尔导数映射，也没有接入 Faster Newton 求解器。对象标签只表示输入网格，不包含上游原始三角面 ID。精确点被映射到源三角面仍不足以推导所有离散 Jacobian：必须进一步验证三角形切分、交点来源、源顶点／交线特征、拓扑与封口分支变化，并在实际 Cherchi 输出映射上做多步长有限差分。任何非唯一、未匹配或标签不一致都不能继续使用为可信导数来源。

本证据属于 `SYNTHETIC_REGRESSION / BOOLEAN_COMPONENT_ONLY`，不代表作者定制 Boolean 版本，也不提高论文案例复现状态。
