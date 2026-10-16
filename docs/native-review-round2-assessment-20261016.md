# 外部 round2 更新包评估（2026-10-16）

来源：`/home/lin/下载/Spin_stage1_native_review_round2`。
当前基线：`b0a13e3f7050630e84e3fbe19e7cf4004540a06e`。
本次仅选择性合入已本地复现并验证的面索引修复，没有整体覆盖项目，没有提交或推送。
17 个交付文件均与包内 CHANGE_MANIFEST.json 的 after_sha256 一致。该校验仅证明包内完整性，不证明算法正确性。

## 逐项结论

| 更新 | 判断 | 实际处理 |
|---|---|---|
| spin.py 的 NumPy 整数面索引、负数/int32 溢出检查 | 明确改进；当前代码可本地复现问题 | 合入该函数修改及原包 17 项测试 |
| GMP 精确坐标导出 ExactPointExport.h、Boolean 适配及 CMake | 有价值的候选；本机独立 C++ 测试通过，但未完成当前生产链集成验证 | 正式后端暂不替换，原包完整保留 |
| 五目标统一入口、双密度配置及导出 | 提供易用入口，但不足以替代当前性能路线 | 暂不安装为正式入口 |
| 配套报告/示例/安装器 | 有参考价值，验收未完成 | 保留在原包，不运行安装器或累计补丁 |

## 本地验证

1. 当前 spin.py 运行原包面索引测试：10 failed / 7 passed。
   八种 NumPy 整数类型被误拒绝；负索引和超过 signed32 的 Python 整数未在该入口拒绝。
2. 合入后，面索引测试和现有 test_spin2014_objective_contract.py：25 passed，0.07 秒。
   使用项目 Python 源码、Anaconda Python 3.13 和 build17 原生模块搜索路径。
3. 在本机以 C++20、O2、frounding-math、项目固定 Indirect_Predicates 头文件和本机 GMP 编译原包 C++ 回归：
   12 项精确舍入和 7 项构造/退化检查通过。二进制为 `/tmp/spin-round2-export-review-20261016`。
   这没有编译/验证完整 Boolean 适配器，不是生产 Jacobian 通过证据。
4. 早先额外尝试 test_spin2014_quadric_interface.py 时未设置原生扩展路径，18 项以 `_make_it_stand` 不可用失败；
   该次共 30 passed / 18 failed，属于测试环境配置失败，不能记为接口回归通过。本轮没有运行完整原生接口套件或整仓测试。

## 为什么不能整体替换

新 `_unified_spin2014_runtime.py::_rotation` 硬编码了
`PYTHON_FRACTION_V1`、`EXACT_NEF_REGION_V1`、`INDEPENDENT_EACH_EVALUATION_V1`、
`RECHECK_INPUT_V1`、`INDEXED_BOUNDARY_REBUILD_V1` 和 `EXACT_GMP_V1`。
这些是现有可用的审计路线，但会重新引入正常迭代的独立检查与重复准备成本；
不能将这套配置认定为优于当前分层检查、证明复用和认证积分配置的快速生产版本。
原包报告也明确记载：旋转验收没有完整终态，第一个 YOYO 用例被 900 秒限制中断，
其余旋转未开始；静态仅数值验收，完整材料未封存。O0 主扩展计时不能用于速度达标。

GMP 导出实现绕过 legacy bignatural 坐标除法，有明确的修复动机，公式审阅和局部本机测试支持继续评估。
它同时将有理数字符串从二进制形式改为十进制形式；当前解析器支持二者，
但摘要/来源身份以及真实交点导数仍需同操作数对照。完整崩溃输入和 sanitizer 日志
不在本文件夹（README 提到的独立 evidence ZIP 未随该文件夹提供），
因此本轮没有本地重现报告中的原始崩溃或独立 CGAL 集合对照。
报告另列 Shewchuk 谓词越界尚未解决；该坐标导出改动不能据此宣称全部内存问题解决。

后续将此原生候选正式替换前，需独立构建并对实际保存操作数比较精确实体、
质量矩、来源、Jacobian、终态重载和耗时；保留当前已冻结模块及历史运行身份。
目前没有证据表明这个包解决了 Heart 收敛、面条完整封存或五目标整体速度门。

## 回退与保护

本次生产代码修改仅限 `spin.py::_faces`，原内容保留在上述 Git 基线。
新增 `tests/python/test_spin_face_integer_inputs.py`。
未修改现有 C++ 源面、Newton、Cherchi 正式模块、运行依赖或用户原有未提交报告。
未删除原包或历史结果。
