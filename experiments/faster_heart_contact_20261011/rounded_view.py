"""Dependency-free interactive input preview; orbiting is not physical motion."""
import argparse
import hashlib
import json
from pathlib import Path
import trimesh


def build(case):
    case = Path(case)
    spec = json.loads((case/'spec.json').read_text())
    meshes = []
    for name, color, kind in (('heart-body.obj', [197, 86, 148], 'body'),
                             ('pin.obj', [219, 199, 154], 'pin'),
                             ('inner-source.obj', [233, 184, 67], 'domain')):
        mesh = trimesh.load(case/name, force='mesh', process=False)
        meshes.append(dict(vertices=(mesh.vertices*25).tolist(), faces=mesh.faces.tolist(),
            color=color, kind=kind, sha256=hashlib.sha256((case/name).read_bytes()).hexdigest()))
    data = json.dumps(meshes, allow_nan=False).replace('</', '<\\/')
    page = '''<!doctype html><html lang="zh-CN"><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>圆润心形 v2 · Faster 输入重建</title><style>
*{box-sizing:border-box}body{margin:0;background:#f5f7fb;color:#24324b;font:16px system-ui,sans-serif}
main{max-width:1160px;margin:32px auto;padding:0 24px}h1{font-size:30px}.sub{color:#59667d;line-height:1.7}
.tag{color:#a03175;font-size:13px;font-weight:700;letter-spacing:1px}.card{background:white;border:1px solid #dce3ec;border-radius:16px;padding:20px;margin:20px 0}
canvas{width:100%;height:580px;display:block;touch-action:none;border-radius:10px;background:linear-gradient(#f9fafc,#e9edf4)}
button{border:1px solid #ced5df;border-radius:7px;background:white;padding:9px 14px;margin:5px;cursor:pointer}
button.active{background:#ad367e;color:white}.metrics{display:flex;gap:32px;flex-wrap:wrap;margin:18px 0}.metrics b{display:block;font-size:21px}
.pending{padding:14px;background:#fff3dd;border-radius:9px;line-height:1.7}img{max-width:100%}
@media(max-width:700px){canvas{height:400px}h1{font-size:25px}}</style>
<main><div class="tag">SPIN-IT FASTER · RECONSTRUCTED METHOD VALIDATION</div>
<h1>圆润心形陀螺：冻结输入 v2</h1>
<p class="sub">拖动观察真实重建网格；滚轮缩放。固定竖直轴、圆润握柄与倾斜心形。此页面展示输入设计，视角转动不是动力学仿真。</p>
<div class="metrics"><div>心形姿态<b>−35°</b></div><div>内域比例<b>0.90</b></div><div>前后厚度倍率<b>1.60</b></div><div>PLA 配置<b>1200 kg/m³</b></div></div>
<div class="card"><div><button id="front">正面</button><button id="side">侧面</button><button id="persp" class="active">透视</button><button id="domain">显示可优化域</button></div>
<canvas id="view"></canvas><p class="sub">粉色为外形；米色为固定轴与握柄。黄色开关显示初始可优化域，未生成优化空腔。</p></div>
<div class="pending"><b>当前状态：输入构造与准入检查待完成。</b><br>
没有优化终态，也没有本例 PyBullet 结果。实际连接、精确分区、隔离距离和质心补偿容量由独立作业检查；现有导数队列结束后顺序执行。</div>
<div class="card"><b>后续求解协议</b><p class="sub">使用 Faster 惯量—质心高度目标，CGAL 二次曲面 Boolean、level 2、预登记种子 0 和 1。优化后的内部材料由实际求解生成，随后独立积分及重复。该输入是重建模型，尚未取得作者原始心形。</p>
<a href="spec.json">冻结配置</a> · <a href="PROGRESS.md">实施记录</a> · <a href="input-preview-shaded.png">四视图预览</a></div>
</main><script>const meshes=__MESHES__;
const canvas=document.getElementById('view'),ctx=canvas.getContext('2d');let yaw=-.65,pitch=.15,zoom=1,showDomain=false,drag=null;
function projection(v){const x=v[0],y=v[1],z=v[2]-27;
const a=Math.cos(yaw)*x-Math.sin(yaw)*y,b=Math.sin(yaw)*x+Math.cos(yaw)*y;
return [a,Math.cos(pitch)*z-Math.sin(pitch)*b,Math.cos(pitch)*b+Math.sin(pitch)*z]}
function draw(){const w=canvas.clientWidth,h=canvas.clientHeight,dpr=devicePixelRatio||1;canvas.width=w*dpr;canvas.height=h*dpr;ctx.setTransform(dpr,0,0,dpr,0,0);ctx.clearRect(0,0,w,h);
const scale=Math.min(w/80,h/80)*zoom;let tris=[];
for(const mesh of meshes){if(mesh.kind==='domain'&&!showDomain)continue;const points=mesh.vertices.map(projection);
for(const face of mesh.faces){const t=face.map(i=>points[i]);const ab=t[1].map((v,i)=>v-t[0][i]),ac=t[2].map((v,i)=>v-t[0][i]);
const n=[ab[1]*ac[2]-ab[2]*ac[1],ab[2]*ac[0]-ab[0]*ac[2],ab[0]*ac[1]-ab[1]*ac[0]],len=Math.hypot(...n)||1;
const shade=.55+.45*Math.max(0,(n[0]*-.25+n[1]*.45+n[2]*-.85)/len);tris.push({t,depth:t.reduce((s,p)=>s+p[2],0)/3,color:mesh.color.map(c=>Math.round(c*shade)),alpha:showDomain&&mesh.kind==='body'?.18:1});}}
tris.sort((a,b)=>b.depth-a.depth);for(const {t,color,alpha} of tris){ctx.beginPath();t.forEach((p,i)=>{const x=w/2+p[0]*scale,y=h/2-p[1]*scale;i?ctx.lineTo(x,y):ctx.moveTo(x,y)});ctx.closePath();ctx.globalAlpha=alpha;ctx.fillStyle='rgb('+color.join(',')+')';ctx.fill();}ctx.globalAlpha=1;
ctx.fillStyle='#54617a';ctx.font='13px system-ui';ctx.fillText('输入预览 · 未优化 · 未进行物理仿真',16,h-18);}
canvas.onpointerdown=e=>{drag=[e.clientX,e.clientY];canvas.setPointerCapture(e.pointerId)};canvas.onpointermove=e=>{if(!drag)return;yaw+=(e.clientX-drag[0])*.008;pitch=Math.max(-1.4,Math.min(1.4,pitch+(e.clientY-drag[1])*.008));drag=[e.clientX,e.clientY];draw()};canvas.onpointerup=()=>drag=null;
canvas.onwheel=e=>{e.preventDefault();zoom=Math.max(.55,Math.min(2,zoom*Math.exp(-e.deltaY*.001)));draw()};
function setCamera(y,p,id){yaw=y;pitch=p;for(const key of ['front','side','persp'])document.getElementById(key).classList.toggle('active',key===id);draw()}
document.getElementById('front').onclick=()=>setCamera(0,0,'front');document.getElementById('side').onclick=()=>setCamera(Math.PI/2,0,'side');document.getElementById('persp').onclick=()=>setCamera(-.65,.15,'persp');
document.getElementById('domain').onclick=()=>{showDomain=!showDomain;document.getElementById('domain').classList.toggle('active',showDomain);draw()};window.onresize=draw;draw();</script></html>'''
    with (case/'index.html').open('x') as out:
        out.write(page.replace('__MESHES__', data))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('case', type=Path)
    build(parser.parse_args().case)
