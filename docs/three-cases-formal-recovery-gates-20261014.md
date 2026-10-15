# 三案例正式恢复前置证据

对应2026-10-13核心与三案例分层验证计划；本文件为清单，不能授予准入。

## 公共入口

当前正式CGAL对照与Cherchi三固定状态证据已存在，cached gen04的146回归通过；小型seed0/1与seed0同版本独立重复均完成，全部科学事件和同次材料一致。旧gen02证据保留，不冒充新版本重复。gen05面积复用候选正在验证，不能提前继承gen04同版本闭环标志。

## 带柄心形

- 冻结请求：rounded-clearance-v3-faster-frozen-v2；输入：bounded-admission-job-03/derived_problem_digest 6b405b3547989c52134f16d04cf11ba7e19bd723c8493afc9409c12758e5d43c。
- 两个失败状态分别在repaired-heart-production-moments-job-01与repaired-heart-exact-production-moments-job-01完成81点数值检查；独立矩误差2.55e-15/1.20e-14；full_material_branch_certification=false。
- 原点记录只含material_vertices/faces和源摘要，缺同次完整精确来源；不能仅靠相近矩值完成分支认证。优先补现有点确实缺失的附件，绑定新正式映射并区分真实分支变化。
- 旧请求外域最大跨度2.57、目标原点0、设计密度1、权重1使用INPUT_MODEL_UNITS_TARGET_ORIGIN_ZERO_V1。正式SpinArtifactProblem默认归一化会改变长度尺度；必须先独立核验坐标、质量矩、Jacobian及高度项转换，不能直接沿用数值权重造成目标变化。
- 源码核对：SpinMomentProblem的Faster目标为惯量比加weight*(r3/r0)。相对于同一个目标原点，若x_model=L*x_normalized，则质量、一次矩、二次矩分别缩放R*L^3、R*L^4、R*L^5；惯量比不变，高度需weight_normalized=weight_model*L。这是待实际材料和导数验证的等价转换，不是改权重寻求成功。轴、目标原点、二次曲面系数及单位同时对应，不能只改一个标量后授予准入。
- 保留握柄/轴杆/底尖和.375mm间距，未完成上述门前不得启动案例完整求解。

## 密集面条

- 81固定部件平衡构造已完成；仅full-input报告不授予准入。后续balanced-shared-export-job-02明确derived_input_partition_passed=true，已有输入OBJ可用证据。
- 既有CGAL准入由noodles_core_run.py组合：请求/清单摘要、原81点数值/source报告、native-tape-certification全部81点及对应构建、小型封存闭环。源码显示其绑定CGAL cgal-proof-capture-build-job-01，不能转授Cherchi。
- 原请求外域最大跨度1、中心有小偏移、目标原点(3.2106657655911306e-6,0,-.48925228022445055)、设计密度1。仍需核验正式入口的精确实体与归一化、目标轴/原点和高度权重对应。
- 接入新正式Cherchi须取得自身同输入材料/来源/导数证据，保留旧真实分支分类和用户中止记录，不恢复旧进程或续拼。

## 花盘

- 分件delivery-job-03：14文件/2归档摘要匹配；table-scene-v1的38交付文件摘要匹配；已有64 MuJoCo轨迹摘要匹配。
- flower-complete.obj/flower-rotor.obj存在，其assembly-check标为同坐标多部件装配，不是轴承焊接实体；需将其与最终物理角色/论文式三层显示统一核对。
- 以上只证明文件身份和存在，优化/重复/几何准入及轨迹惯量的完整审查仍独立。

三案例只在各自正式准入成立后串行运行/交付。新的原始输入、作者CAD、制造资格、全参数可行性均未由该清单建立。
