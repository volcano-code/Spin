# Spin-It 2014 后续 BBW 验收契约核查

核查日期2026-09-10，基线59ed96c；本轮只实现细网格及一般场插值，不宣称已实现下述求解器。原始Spin-It第6节/Figure4要求较细网格权重及表面/新增单元顶点的三线性传递。

## 需区分的论文路径

[BBW 2011 原文](https://igl.ethz.ch/projects/bbw/bounded-biharmonic-weights-siggraph-2011-compressed-jacobson-et-al.pdf)第3.1、3.4节同时给出完整耦合问题和逐handle独立求解后归一化的加速变体。归一化本身不是“论文禁止的简化”，但两种结果不能声称具有同一最优性证明。Eq.8是包含质量矩阵的FEM能量；Eq.4还要求cage面上的线性边界条件，不仅仅固定角点。

后续应分别记录原始列解、归一化变化、最终权重与各自适用的残差。完整耦合问题检查其等式/盒约束及KKT；独立列变体检查原始列KKT，并测量后处理造成的变化，不能把原始解的最优性转授给归一化结果。尚无证据确定Spin-It作者采用哪个具体变体及全部网格参数，不能冒称作者原值。

## 当前源码的具体接口风险

本地libigl提交`f06fb13014aadff7eb1d7423ac018620441a3657`：

- `libigl/include/igl/bbw.cpp` 默认`partition_unity=false`，明确断言未实现true路径；release构建可能禁用assert，不能只把该开关设true就宣称完成耦合约束。
- 同一函数遇到`SOLVER_STATUS_MAX_ITER`只输出警告，其最终bool不是“每列已收敛”的充分证据。严格适配器必须保留每列实际终态并独立检查。
- `libigl/include/igl/harmonic.cpp` 的几何路径使用cotmatrix与massmatrix，二阶能量含逆质量矩阵。项目旧`BoundedBiharmonicWeights.cpp`是无权邻接图LᵀL，不能未经核验称为同一个离散问题。
- 旧项目函数的`converged`来自归一化前列状态，之后截断、归一化、重写约束；旧diagnostics验证bounds/行和/固定值，不含最终最优性检查。保留兼容但不赋予新的严格身份。

核查源码SHA-256：

```text
57b5c884bcb821cacdd743345d9bb26aa8d92d46913caf5b73e7b8e03c5a6ad1  libigl/include/igl/bbw.cpp
1d5ac07831f92c3a42b98c14b89c5759fd77e93646e8e19ba01edbf4f3c7f242  libigl/include/igl/harmonic.cpp
```

## 必须完成的后续工作

1. 从本轮细网格建立一致体积离散、质量矩阵、刚度矩阵和独立能量检查；固定单元细分/边界条件，未知作者设置注明来源。
2. 让全部cage面上的细节点获得明确线性边界条件，不将完整角点集合误当作完整边界约束。
3. 实现并分别验收适用的耦合/独立列求解协议，保存失败与后处理差异。归一化不是隐形修复步骤。
4. 将实际已验收权重传给原始表面和新增叶节点，再接beta/phi/全部handles联合优化与几何回滚。

本轮一般插值器对多重线性场的精确再现不证明BBW权重具有相同的线性精度；两项测试的数学对象不同。原论文案例缺失及完整复现发布门没有改变。
