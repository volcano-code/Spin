# 封存材料的原生 OpenGL 查看器（Task7v1）

本次接通了「实际 Spin-It 2014 源表面求解 → 独立材料封存 → 原生 OpenGL」流程。开发目录仍为 `/data/xvxinyao/make-it-balance-paper-sim`，没有整合或覆盖主目录的未提交修改。

当前展示对象是 **SYNTHETIC_REGRESSION 解析方体**，用于检查真实求解、挖空、质量矩和显示是否连接正确，不是论文中的原始大象、心形或其他实验模型。本次没有增加 Faster 原论文案例完成数；两篇论文整体仍为 `PARTIAL_BLOCKED`。

## 打开这次结果

在本机真实 X11 桌面运行以下命令；窗口已打开时直接切换到标题为 `Sealed spin | HOLLOW_ONLY | PARTIAL_BLOCKED` 的窗口即可。

```bash
cd /data/xvxinyao/make-it-balance-paper-sim
env -u LD_PRELOAD -u LD_LIBRARY_PATH DISPLAY=:1 \
  build/task-7v1-viewer/make-it-stand --sealed-spin-bundle \
  build/task-7v1-source-run-1/sealed \
  2f45a44a9c74fd3d02bc68a7c36e11b8d7fdee0a0fe22a9857236338d6ed3929
```

最后一个参数是从原输入独立锁定的 problem digest，不应从未知包的 manifest 中读取后自动信任。此启动命令只适用于当前对应构建和封存结果；更改参与计算或检查的源码／扩展后，旧包可能被拒绝，应在新目录重新运行，不能改旧摘要绕过检查。

将 `--sealed-spin-bundle` 换成 `--check-sealed-spin-bundle` 可无窗口重新积分并检查。这个原生入口不需要为 Conda 预加载 libstdc++；当前 Python 测试驱动仍使用单独记录的运行时设置。

## 如何看 carving

- 左侧是原始模型，右侧是封存后的材料；两侧共享观察相机，模型坐标和姿态不改写。
- `material`：金色，全部保留材料。`cavity`：红色，被挖掉区域的边界表示；红色不是新增实体。
- 关闭 `Display cutaway`：可以比较完整外表面，内部空腔会被不透明外壁遮挡，这是正常现象。
- 打开 `Display cutaway`：仅隐藏显示三角面，便于看内部。当前不是精确截面运算：跨截面的三角面也会被隐藏，没有补截面盖。画面上缺少的外壁不代表算法删除了外壁。
- `Isolate cavity`：只看挖空区域。`Explode cavity (display only)`：将红色显示副本移开；不会移动实际空腔或保存修改。
- `Fixed shell material` 与 `Outer surface (not shell)` 分开：前者是有厚度的固定材料，后者只是外边界。`retained` 是设计域中剩余的材料，不能与全体 material 重复计算质量。
- 红点为整体质心，绿线只表示唯一确定的惯量主轴，紫线为目标轴，黄点为输入支撑点。本例横向两个特征值重复，所以只显示唯一的 z 主轴；一个支撑点不代表通过了有限面积支撑多边形检查。

状态栏始终保留几何验证范围、`UNKNOWN`／`NOT_VERIFIED` 和 `PARTIAL_BLOCKED`。窗口不进行优化、轨迹认证或物理仿真。视觉旋转是相机交互，不是陀螺动力学。

## 本次实际验证结果

- 主构建全部已注册 CTest：199/199 通过，518.63 秒；最终原生加载用例独立重跑 3/3 通过。
- 两次从源表面实际求解均完成：各 39 条求解记录、72 个材料挖空 cuboid，几何最终合并为 1 个实际空腔。
- 两次 `inputs.json`、`native-run.json`、`scientific.json` 分别逐字节一致。
- 最终材料质量 `503/64 = 7.859375`，质心 `(0,0,2003/2012)`；独立重新积分的归一化质量矩误差最大约 `1.39e-17`。此误差只针对这个解析回归，不代表所有论文模型精度。
- 真实桌面 `DISPLAY=:1`、物理 DP-0 显示器、NVIDIA OpenGL 上验证了完整材料、剖切、空腔隔离和显示副本移位。交互前后封存包的 18 个文件摘要完全相同。

本机结果目录：`build/task-7v1-source-run-1`、`build/task-7v1-source-run-2`；真实桌面截图：`build/task-7v1-desktop/full-material.png`、`isolated-cavity.png`、`exploded-cavity.png`。这些是实际窗口观测，可能含宿主桌面浮层，不是离线论文渲染。

## 尚未完成

这个增量只完成封存几何的原生显示，不关闭 Task7/M6，更不授予 `VERIFIED_COMPUTATIONAL_REPRODUCTION`。完整 cage／BBW 联合变形、双材料严格链路、Faster 剩余几何／二维／浮力／CAD、原始资产和关键参数、原论文数值矩阵、认证迭代回放及独立仿真仍需各自通过验收。后续应继续补算法和验证缺口，不能用本例展示替代原论文实验。
