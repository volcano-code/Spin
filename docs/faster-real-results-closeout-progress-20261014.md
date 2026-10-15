# Faster真实结果收尾执行账本

Spec / plan: /home/lin/.codex/attachments/243f29c0-ec04-4952-8caf-a9981051803a/goal-objective.md

- 上一中断回合无可确认实现变更，按no progress处理；本轮重新检查实际PID343736存活，面条已数值收敛而终态尚未完成，不重启。
- 已有linked worktree spin-it-strict-computational-v2；保留全部既有修改。无自动提交、合并或删除。
- 执行方式：inline，按executing-plans与TDD。前置依赖：调度保护须在后续大型作业入口使用；求解模块复用改动须等待冻结作业结束。
- Ruling: 本轮先新增独立调度模块，不修改当前运行绑定的source_job.py或求解Python；后续入口统一调用包装器。当前历史入口不能靠新锁自动获得串行保证，必须同时扫描真实systemd单位与控制组。
- Ruling: 使用项目锁加持久pending记录，缺少终态或观察失败均拒绝新作业；不靠PID消失判定完成。资源终止可释放调度但不授予科学成功。
- 测试首先覆盖跨进程互斥、控制器退出后工作进程存活、结果未落盘、资源终止、控制组仍非空及观察失败。尚未实现或通过。
- 用户批准仅重验受影响内容；不运行无关全论文/2014大型套件。完整核心与案例门保持未完成。

## 本轮执行证据

- 调度新模块`serial_source_job.py`已实现跨进程flock、持久pending注册、真实systemd单位扫描及cgroup.events populated复核；缺终态、单位身份不符、观察失败均拒绝新作业。原运行资源模块未修改。
- 红测：新模块缺失导致测试收集失败；实现后13项测试通过，包括锁持续至结果发布、控制器异常后保留待核对记录。pytest已有asyncio配置warning保留。
- 实际非模拟现场检查调用observe_units/assert_project_idle，正确返回PROJECT_WORKLOAD_ACTIVE: spin2014-job-8d13aea11ed547d4919b11ef6616af8a.service；PID343736仍真实运行。未启动新大型作业。后续launcher尚待迁移，不能声明所有旧入口已受保护。
- 面条evaluation15已核对：accepted、raw=6.698726400984546e-14、mu=249.53355510240303>0、objective=1.437814460042704，物理相对残差2.285199499561518e-17；attempt_terminal reason=converged。native_newton_solver为569.557492111秒，后续互补Boolean已返回；终态完整材料报告仍缺失。已有接受附件保留，不用终态初始化新优化。
- Heart正式准备PROCESS_EXITED/0，request等价，两个困难状态有理材料一致、Jacobian差0、矩差3.885780586188048e-16；optimization_invocations=0，case_admission=false。不把准备通过标成优化成功。
- 后续动作：等待当前作业真实终态；若超资源上限保留原状态，仅为实际保存终态安排独立后处理诊断。待冻结依赖可修改后迁移后续入口至新调度器，并实施一次准备复用。当前作业重叠计时不纳入严格性能验收。
- 本轮联合回归：test_serial_source_job、test_source_checkpoint_store及v2共82 passed/4.76秒；未运行全项目或新大型实验，不扩大通过范围。最后真实面条PID343736为56:17、CPU99.7%、stderr为空，仍非终态。
- PreparedBoundary此前串行链微基准已结束：两个保存材料输出逐位一致；2856面材料legacy约0.0743秒、prepared约0.0681秒，仅首128面算术选择范围。收益有限且与面条运行重叠，不作为正式性能验收；优先一次输入准备，而非据此宣称数量级加速。
