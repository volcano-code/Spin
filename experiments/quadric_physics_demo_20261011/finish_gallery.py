"""Render saved results only; no optimization or altered physics."""
from pathlib import Path
import json, hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.path import Path as MPath
from matplotlib.patches import PathPatch, Circle

P=Path(__file__).resolve().parent/'multiple-cases-v1'
def read(p): return json.loads((P/p).read_text())
a,b=[read(f'flower-run-{i}/application.json') for i in (1,2)]
assert a['scientific_hash']==b['scientific_hash']
checks=[read(f'flower-check-{i}/output/scientific.json') for i in (1,2)]
for check in checks:
    assert all(t['moments_match'] and t['balance_constraints_pass'] for t in check['cases'][0]['attempts'])
rows=a['cases'][0]['attempts'];s=rows[0]['solver_result'];g=read('flower-geometry.json')
def compound(ax,loops,color):
    vertices=[];codes=[]
    for loop in loops:
        q=np.asarray(loop);vertices.extend([*q,q[0]]);codes.extend([MPath.MOVETO]+[MPath.LINETO]*(len(q)-1)+[MPath.CLOSEPOLY])
    ax.add_patch(PathPatch(MPath(vertices,codes),facecolor=color,edgecolor='#28364d',lw=.7))
def sample(cycle):
    pts=[]
    for seg in cycle:
        if seg['type']=='LINE_SEGMENT':pts.append(seg['start'])
        else:
            t=seg['start_angle']+np.linspace(0,seg['sweep'],max(3,int(abs(seg['sweep'])*50)),endpoint=False)
            pts.extend(np.asarray(seg['center'])+seg['radius']*np.column_stack((np.cos(t),np.sin(t))))
    return pts
fig,axs=plt.subplots(1,3,figsize=(15,5),facecolor='#edf2f8')
for ax in axs[:2]:
    ax.set_aspect('equal');ax.set_xlim(-1.15,1.5);ax.set_ylim(-1.35,1.35);ax.axis('off')
compound(axs[0],[g['outer'],g['bore'][::-1]],'#e2a953')
compound(axs[1],[g['outer'],g['inner'][::-1]],'#a8b7ca')
compound(axs[1],[sample(c) for c in s['boundary_cycles']],'#28a79b')
compound(axs[1],[g['bearing'],g['bore'][::-1]],'#a8b7ca')
axs[1].add_patch(Circle(s['center'],s['radius'],fill=False,ls='--',ec='#d74470',lw=1.8))
for ax in axs[:2]:ax.plot(0,0,'+',color='#182840',ms=12,mew=2)
axs[0].set_title('Original eccentric flower + bearing bore')
axs[1].set_title('Actual dogleg result / circular cut boundary')
axs[2].axis('off')
err=max(t['maximum_normalized_moment_error'] for t in checks[0]['cases'][0]['attempts'])
axs[2].text(0,.95,f'Faster planar application\nNot the author\'s original Flower\n\nSeeds 0 / 1: CONVERGED\nEvaluations: 31 / 36\nAccepted steps (seed 0): {s["accepted_steps"]}\nTwo processes: scientific hashes match\nIndependent moment error: {err:.2e}\n\nPolygon + disk, not voxel variables\nBearing ring and fixed rim preserved\nSaved analytic extrusion: thickness 0.15\n\nFull topology / manufacturing:\nNOT APPROVED',va='top',fontsize=12,linespacing=1.5)
fig.tight_layout();fig.savefig(P/'flower.png',dpi=150);plt.close(fig)
report={'classification':'RECONSTRUCTED_METHOD_VALIDATION','scientific_hash':a['scientific_hash'],'repeat_match':True,'independent_moment_checks':True,'full_topology_verified':False,'manufacturing_approved':False,'attempts':[{'seed':t['seed'],'status':t['status'],'accepted_steps':t['solver_result']['accepted_steps'],'evaluations':t['function_evaluations']} for t in rows]}
(P/'flower-comparison.json').write_text(json.dumps(report,indent=2))
html='''<!doctype html><meta charset="utf-8"><title>Faster 方法：多案例结果</title><style>body{max-width:1150px;margin:40px auto;background:#101b2b;color:#e9eff7;font:17px/1.7 system-ui;padding:20px}section{background:#1c2b40;padding:24px;border-radius:16px;margin:25px 0}img,video{width:100%;border-radius:8px}a{color:#71ded0}.warning{color:#ffce83}h1{font-size:36px}</style>
<h1>二次曲面方法 · 多案例结果</h1><p>全部为重建方法验证，不是作者原始模型的逐例复现。旧结果重新读取；五瓣指尖陀螺为本次真实新求解。</p>
<section><h2>01 五瓣指尖陀螺｜Faster 二维算法</h2><img src="flower.png"><p>实际带孔 polygon–disk + dogleg；固定外围及轴承环。种子 0、1 均收敛，两进程重复一致，保存后独立积分通过。虚线圆是实际求解边界，不是手工画出的目标空腔。当前未通过完整制造拓扑验收。</p><a href="flower-comparison.json">重复与检查摘要</a> · <a href="flower-run-1/application.json">实际求解及解析挤出实体</a></section>
<section><h2>02 Cuboid 悠悠球｜二次曲面材料优化</h2><img src="cuboid.png"><p>Spin-It 2014 悠悠球目标 + Faster 风格二次曲面 Newton。两个设计域与最终材料来自保存结果，不是体素结果平滑。</p><details><summary>PyBullet 无外力矩旋转诊断：未通过守恒检查</summary><video controls preload="metadata" src="yoyo-physics/output/yoyo-free-rotation.mp4"></video><img src="yoyo-physics/output/axis-curves.png"><p class="warning">该视频是实际计算，但原始物体出现严重能量及角动量漂移，不能用它证明优化优越性。无重力、无绳索、无持续驱动；不是悠悠球沿绳上下运动。完整失败数据保留。</p><a href="yoyo-physics/output/summary.json">四场景守恒误差</a></details></section>
<section><h2>03 静态平衡｜退化二次曲面（平面）</h2><img src="static.png"><p>真实保存材料重新积分。2014 平方质心高度目标，平面专用求解路径；重建多面体，不是作者 tilted vase。这里展示静态几何与质量属性，不伪装为动态仿真。</p></section>
<section><h2>04 Heart｜原有旋转对比及更正</h2><img src="../presentation-v1/cutaway.png"><p>质心高度 25.715 → 25.849 mm，确实上升 0.134 mm。R 从 2.333 降至 2.091，R+h² 从 3.391 降至 3.160；联合目标改善不等于高度单独下降。此旧输入没有握柄，不能在优化后加柄并沿用旧惯量。</p><a href="../presentation-v1/index.html">原有四面板接触仿真（时间步未收敛）</a></section>
<section><h2>与 Spin-It 体素方案的区别</h2><p>本页的设计变量是二次曲面系数、平面参数或二维圆参数，而不是为大量体素分配填充率。可见的三角网格只是几何表示，不代表使用了体素优化。Heart、Cuboid、静态采用2014目标的二次曲面重实现；新 Flower 则运行 Faster 的二维应用算法。没有在相同输入上执行两种优化，因此不宣称速度或效果优于2014。</p></section>'''
(P/'index.html').write_text(html)
files=['index.html','flower.png','flower-comparison.json','cuboid.png','static.png','yoyo-physics/output/summary.json','yoyo-physics/output/yoyo-free-rotation.mp4']
(P/'gallery-manifest.json').write_text(json.dumps({f:hashlib.sha256((P/f).read_bytes()).hexdigest() for f in files},indent=2))
print(P/'index.html')
