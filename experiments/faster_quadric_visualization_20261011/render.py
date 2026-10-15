"""Read-only visualization of a saved Faster/Cherchi method run.

No solver, repair, smoothing, or new scientific acceptance is performed.
OBJ coordinates are display-only derivatives of the archived exact materials.
"""
from pathlib import Path
import base64
import hashlib
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.collections import LineCollection
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
SOURCE = ROOT / "experiments/spin_it_reproduction/faster/results/cherchi-centered-axis-index-binding-20261001-seed1-01/solver-output"
REPEAT = SOURCE.parent.parent / "cherchi-centered-axis-index-binding-20261001-seed1-02/solver-output"
FONT = "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc"
font_manager.fontManager.addfont(FONT)
plt.rcParams.update({"font.family": font_manager.FontProperties(fname=FONT).get_name(), "axes.unicode_minus": False,
                     "font.size": 10, "axes.titleweight": "bold", "savefig.facecolor": "#f6f8fc"})
TEAL, ORANGE, NAVY, GRAY = "#078d92", "#e67838", "#142944", "#8e9bac"
tracked = {}


def read_bytes(path):
    data = path.read_bytes()
    tracked[str(path)] = hashlib.sha256(data).hexdigest()
    return data


def read_json(path):
    return json.loads(read_bytes(path))


def mesh(path):
    vertices, faces = [], []
    for line in read_bytes(path).decode().splitlines():
        fields = line.split()
        if fields and fields[0] == "v":
            vertices.append([float(x) for x in fields[1:4]])
        elif fields and fields[0] == "f":
            indices = [int(x.split("/")[0])-1 for x in fields[1:]]
            assert len(indices) == 3, "This display does not retriangulate saved meshes"
            faces.append(indices)
    v, f = np.asarray(vertices), np.asarray(faces)
    assert np.isfinite(v).all() and f.min() >= 0 and f.max() < len(v)
    return v, f


def box(ax, radius, color=GRAY, alpha=1):
    from itertools import product
    for p in product((-radius, radius), repeat=3):
        for k in range(3):
            if p[k] < 0:
                q = list(p); q[k] = radius
                ax.plot(*zip(p, q), color=color, alpha=alpha, lw=1)


def clip_triangles(v, f, axis=1, value=.0):
    """Display cut only. Never persisted or presented as a closed output mesh."""
    out = []
    for tri in v[f]:
        polygon = []
        for a, b in zip(tri, np.roll(tri, -1, axis=0)):
            ia, ib = a[axis] <= value, b[axis] <= value
            if ia:
                polygon.append(a)
            if ia != ib:
                polygon.append(a + (b-a)*(value-a[axis])/(b[axis]-a[axis]))
        if len(polygon) >= 3:
            for i in range(1, len(polygon)-1):
                out.append([polygon[0], polygon[i], polygon[i+1]])
    return out


def section(v, f, plane=.017):
    segments = []
    for tri in v[f]:
        hits = []
        for a, b in zip(tri, np.roll(tri, -1, axis=0)):
            if (a[1] < plane) != (b[1] < plane):
                p = a+(b-a)*(plane-a[1])/(b[1]-a[1])
                hits.append(p[[0, 2]].tolist())
        if len(hits) == 2:
            segments.append(hits)
    assert segments
    return np.asarray(segments)


def polynomial(q, u, v, w):
    # QuadricCoefficients10::quadraticMatrix: off-diagonal entries are q[7:10].
    return (q[0]+q[1]*u+q[2]*v+q[3]*w+q[4]*u*u+q[5]*v*v+q[6]*w*w
            +2*q[7]*v*w+2*q[8]*u*w+2*q[9]*u*v)


def setup3d(ax, lim=.29):
    ax.set(xlim=(-lim, lim), ylim=(-lim, lim), zlim=(-lim, lim),
           xlabel="u", ylabel="v", zlabel="w / 目标轴")
    ax.set_box_aspect((1, 1, 1)); ax.view_init(elev=24, azim=-58)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])
    ax.set_facecolor("#f6f8fc"); ax.grid(False)


def show_material(ax, data, title):
    v, f = data
    tris = clip_triangles(v, f)
    a = np.asarray(tris)
    normals = np.cross(a[:, 1]-a[:, 0], a[:, 2]-a[:, 0])
    lengths = np.linalg.norm(normals, axis=1)
    normals /= np.maximum(lengths[:, None], 1e-30)
    shade = .5+.5*np.abs(normals @ np.array([.2, -.6, .77]))
    colors = np.column_stack((.02+.08*shade, .37+.24*shade, .4+.22*shade, np.ones(len(a))))
    ax.add_collection3d(Poly3DCollection(tris, facecolors=colors, edgecolors="none", rasterized=True))
    box(ax, .25, GRAY, .6)
    ax.plot([0, 0], [0, 0], [-.3, .3], "--", color=NAVY, lw=1.4)
    setup3d(ax); ax.set_title(title, pad=14)


def show_section(ax, q, segments, title):
    x = np.linspace(-.25, .25, 301)
    u, w = np.meshgrid(x, x)
    z = polynomial(q, u, .017, w)
    ax.contourf(u, w, z, levels=[-1e3, 0, 1e3], colors=["#ffe9d9", "#def4f1"])
    ax.contour(u, w, z, levels=[0], colors=ORANGE, linewidths=1.5)
    ax.add_collection(LineCollection(segments, colors=TEAL, linewidths=1.7))
    ax.axvline(0, color=NAVY, ls="--", lw=.8)
    ax.set(xlim=(-.27, .27), ylim=(-.27, .27), xlabel="u", ylabel="w / 目标轴", title=title)
    ax.set_aspect("equal"); ax.grid(alpha=.12)


def main():
    r = read_json(SOURCE / "result.json")
    rr = read_json(REPEAT / "result.json")
    inputs = read_json(SOURCE / "inputs.json")
    trace = r["trace"]
    assert r["status"] == rr["status"] == "CONVERGED"
    assert trace == rr["trace"], "Archived repeated scientific trace changed"
    assert r["sealed_artifact"]["scientific_digest"] == rr["sealed_artifact"]["scientific_digest"]
    assert inputs["mode"] == "YO_YO"
    assert r["external_terminal_geometry_audit"]["float_validation_status"] == "REJECTED"
    assert trace[-1]["raw_equation_residual"] <= 1e-12
    # Failed evaluations remain in the archive; only accepted states have a mesh.
    evals = {e["quadric_hash"]: e for e in r["evaluations"] if "quadric_hash" in e and "oriented_output_path" in e}
    frames = []
    for t in trace:
        e = evals[t["quadric_hash"]]
        index = e["evaluation_index"]
        data = mesh(SOURCE / f"material-{index:03d}-oriented.obj")
        frames.append((t, e, data, section(*data)))
    design = mesh(SOURCE / "design.obj")
    assert np.allclose(design[0].min(0), -.25) and np.allclose(design[0].max(0), .25)
    # Bind exact sources without interpreting a display conversion as exact geometry.
    for suffix in [".exact.tsv", ".faces.tsv", ".construction.tsv"]:
        read_bytes(SOURCE / ("material-010.obj"+suffix))
    read_bytes(SOURCE / "sealed-material-artifact/manifest.json")
    improvement = 100*(trace[0]["objective"]-trace[-1]["objective"])/trace[0]["objective"]
    fig = plt.figure(figsize=(17, 10), facecolor="#f6f8fc")
    fig.suptitle("Spin-It Faster  /  二次曲面决定挖空边界", x=.055, y=.96, ha="left", fontsize=24, color=NAVY)
    fig.text(.055, .905, "真实 Cherchi Boolean 悠悠球方法案例  ·  固定外壳不变  ·  10 个归一化系数 / 9 个独立方向", color=NAVY, fontsize=12)
    gs = fig.add_gridspec(2, 3, left=.055, right=.96, top=.855, bottom=.16, hspace=.40, wspace=.30)
    ax = fig.add_subplot(gs[0, 0], projection="3d")
    box(ax, .5, GRAY); box(ax, .25, TEAL)
    ax.plot([0, 0], [0, 0], [-.6, .6], "--", color=ORANGE, lw=2)
    setup3d(ax, .58); ax.set_title("① 输入：外壳 + 明确的可优化域")
    ax.text2D(.03, -.13, "灰框：原外表面；内框：设计域\n两者之间的材料固定不动", transform=ax.transAxes, color=NAVY, fontsize=10)
    show_material(fig.add_subplot(gs[0, 1], projection="3d"), frames[0][2], "② 随机初始候选：实际保留材料")
    show_material(fig.add_subplot(gs[0, 2], projection="3d"), frames[-1][2], "③ 优化后：实际保留材料")
    show_section(fig.add_subplot(gs[1, 0]), trace[0]["quadric_coefficients"], frames[0][3], "初始截面  /  v = 0.017")
    show_section(fig.add_subplot(gs[1, 1]), trace[-1]["quadric_coefficients"], frames[-1][3], "最终截面  /  v = 0.017")
    ax = fig.add_subplot(gs[1, 2])
    ax.semilogy([t["iteration"] for t in trace], [t["raw_equation_residual"] for t in trace], "o-", color=TEAL, lw=2)
    ax.axhline(1e-12, color=ORANGE, ls="--", label="严格阈值 1e-12")
    ax.set(xlabel="已接受状态编号（不是物理时间）", ylabel="原始方程残差", title="球面 Newton：保存的真实轨迹")
    ax.grid(alpha=.18); ax.legend(fontsize=9)
    fig.text(.055, .067, f"目标值  {trace[0]['objective']:.6f} → {trace[-1]['objective']:.6f}（下降 {improvement:.2f}%）    |    9 次接受更新    |    11 次求值    |    残差 {trace[-1]['raw_equation_residual']:.2e}", fontsize=12, color=NAVY)
    fig.text(.055, .033, "青色：保留材料 / 实际网格截线；橙色：解析二次曲面边界；浅色：解析场示意。3D 为显示剖切，外壳在②③中隐藏。\n历史方法验证，非作者原模型。精确材料检查通过；浮点派生网格被拒绝，仅供显示，非制造或物理仿真结果。", fontsize=9, color="#52627a")
    fig.savefig(OUT / "overview.png", dpi=150); plt.close(fig)
    html_frames = []
    for t, e, data, seg in frames:
        fig = plt.figure(figsize=(11, 5), facecolor="#f6f8fc")
        ax = fig.add_subplot(121, projection="3d")
        show_material(ax, data, f"状态 {t['iteration']}：实际 Boolean 保留材料")
        show_section(fig.add_subplot(122), t["quadric_coefficients"], seg, "二次曲面与实际材料截线")
        fig.tight_layout(pad=2)
        path = OUT / f"state-{t['iteration']:02d}.png"
        fig.savefig(path, dpi=110); plt.close(fig)
        html_frames.append({"iteration": t["iteration"], "objective": t["objective"],
                            "residual": t["raw_equation_residual"], "coefficients": t["quadric_coefficients"],
                            "evaluation": e["evaluation_index"], "quadric_hash": t["quadric_hash"],
                            "image": "data:image/png;base64,"+base64.b64encode(path.read_bytes()).decode()})
    template = '''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Spin-It Faster · 二次曲面材料优化</title><style>
*{box-sizing:border-box}body{margin:0;background:#f6f8fc;color:#142944;font:16px system-ui,"Noto Sans CJK SC",sans-serif}main{max-width:1180px;margin:38px auto;padding:0 24px}h1{font-size:34px;margin:12px 0}h2{margin-top:0}.tag{color:#078d92;font-weight:700;letter-spacing:2px}.sub{color:#52627a;line-height:1.8}.card{background:white;border:1px solid #dce4ed;border-radius:16px;padding:24px;margin:24px 0}.grid{display:grid;grid-template-columns:1.5fr 1fr;gap:28px}img{width:100%;height:auto}button{border:0;border-radius:8px;padding:12px 20px;background:#078d92;color:white;font:inherit;cursor:pointer}input{width:70%;vertical-align:middle}table{width:100%;border-collapse:collapse;font-size:14px}td,th{padding:7px;text-align:left;border-bottom:1px solid #e2e9ee}.metrics{display:flex;gap:30px;flex-wrap:wrap;margin:20px 0}.metrics b{display:block;font-size:24px;color:#078d92}code{overflow-wrap:anywhere}.note{border-left:4px solid #e67838;padding-left:16px;font-size:14px;line-height:1.8}@media(max-width:800px){.grid{grid-template-columns:1fr}h1{font-size:26px}}</style>
<main><div class="tag">SPIN-IT FASTER / ACTUAL SAVED OPTIMIZATION</div><h1>改变二次曲面，重新分配内部材料</h1>
<p class="sub">盒形悠悠球方法夹具 · Cherchi Boolean · 非零固定外壳 · 历史独立重复结果<br>这是优化过程的回放，不是旋转动力学。模型使用方法验证单位，不是实物打印规格。</p>
<div class="card"><h2>方法总览</h2><img src="overview.png" alt="Faster 方法总览"><p>优化变量是 <b>p(u,v,w)=a+b·x+xᵀAx</b> 的 10 个系数；归一化后有 9 个独立方向。连续二次曲面先三角化，再与设计域做真实 Boolean。固定外壳始终不参与削减。</p></div>
<div class="card"><h2>拖动查看真实接受状态</h2><button id="play">播放优化回放</button> <input id="slider" type="range" min="0" max="9" value="0" step="1"><span id="step"></span>
<div class="metrics"><div>目标值<b id="objective"></b></div><div>原始方程残差<b id="residual"></b></div><div>实际求值索引<b id="evaluation"></b></div></div>
<img id="frame" alt="实际材料与二次曲面截面"><div class="grid"><div><h3>如何看图</h3><p class="sub">青色实体及截线来自该状态保存的实际 Boolean 网格。橙线是相同系数的解析 p=0 截线，淡橙表示挖除侧，淡青表示保留侧。解析场与离散材料分别显示，不声称两者完全重合。</p><p class="sub">3D 仅显示设计域内材料的一半，隐藏固定外壳便于观察。显示剖切不产生新的合格网格。状态 0 是随机初始候选，已经包含挖空，并非原始实心输入。</p><p class="note">最终形态接近双锥极限。系数表保留完整数值，不能把四舍五入后的锥面公式当作实际终态。当前回放不新增曲面退化或通用收敛保证。</p></div><div><h3>当前二次曲面系数</h3><table><thead><tr><th>项</th><th>系数</th></tr></thead><tbody id="coeff"></tbody></table></div></div></div>
<div class="card"><h2>为什么这是 Faster 方法</h2><p>本例调用保存的 Faster YO_YO 问题、球面 Newton 和 Cherchi 材料求值链。每次更新改变全局二次曲面，再用其与设计域的交集/差集得到材料边界。图中的三角形用于几何计算与显示，不是逐体素填充变量；没有调用 Spin-It 2014 的体素 SLQP、自适应八叉树或 cage 形变。</p><p>这是 Faster 的方法验证例子，不是作者论文原始模型，也不等同于当前“2014 五目标适配”新版本已完成验收。</p><p class="note">证据范围：历史精确材料验收及独立重复已通过；浮点派生网格验收为 REJECTED。本页面只是显示派生物，不可用于制造或接触仿真。完整轨迹含 11 次求值；滑块仅展示初始状态与 9 个接受状态，失败求值仍保留在源日志。</p><a href="manifest.json">查看来源摘要与证据清单</a> · <a href="overview.png">打开高清总览</a></div></main>
<script>const frames=__DATA__;const names=['1','u','v','w','u²','v²','w²','2vw','2uw','2uv'];let timer=null;const slider=document.getElementById('slider');function draw(){const d=frames[Number(slider.value)];document.getElementById('step').textContent=d.iteration+' / 9';document.getElementById('objective').textContent=d.objective.toFixed(9);document.getElementById('residual').textContent=d.residual.toExponential(3);document.getElementById('evaluation').textContent=d.evaluation;document.getElementById('frame').src=d.image;document.getElementById('coeff').innerHTML=d.coefficients.map((c,i)=>'<tr><td>'+names[i]+'</td><td>'+c.toExponential(8)+'</td></tr>').join('');}slider.oninput=draw;document.getElementById('play').onclick=()=>{if(timer){clearInterval(timer);timer=null;document.getElementById('play').textContent='播放优化回放';return}document.getElementById('play').textContent='暂停';timer=setInterval(()=>{slider.value=(Number(slider.value)+1)%frames.length;draw()},1500)};draw();</script></html>'''
    (OUT / "index.html").write_text(template.replace("__DATA__", json.dumps(html_frames)), encoding="utf-8")
    manifest = {"schema": "saved-faster-method-visualization-v1", "source_run": str(SOURCE),
                "repeat_run": str(REPEAT), "classification": inputs["classification"],
                "mode": inputs["mode"], "is_physics_simulation": False, "is_author_model": False,
                "display_only_float_mesh": True, "float_validation_status": "REJECTED",
                "coordinate_frame": "MANIFEST_TARGET_NORMALIZED: axis is w; original world axis is x",
                "coefficient_convention": "[a,bu,bv,bw,Auu,Avv,Aww,Avw,Auw,Auv]; off-diagonal terms doubled",
                "source_normal_subdivision": inputs["quadric_icosphere_subdivisions"],
                "accepted_updates": len(trace)-1, "function_evaluations": r["function_evaluations"],
                "initial_objective": trace[0]["objective"], "final_objective": trace[-1]["objective"],
                "objective_decrease_percent": improvement, "raw_residual": trace[-1]["raw_equation_residual"],
                "repeated_trace_equal": True, "repeated_material_digest_equal": True,
                "original_scope_field": r["scope"], "final_geometry_verified_field": r["final_geometry_verified"],
                "material_digest": r["sealed_artifact"]["scientific_digest"], "source_sha256": tracked,
                "renderer_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                "outputs": {f.name: hashlib.sha256(f.read_bytes()).hexdigest() for f in OUT.iterdir() if f.suffix in [".png", ".html"]}}
    (OUT / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2)+"\n")
    print(json.dumps({"output": str(OUT), "accepted_states": len(frames), "objective_decrease_percent": improvement,
                      "raw_residual": trace[-1]["raw_equation_residual"], "trace_repeat_equal": True}))


if __name__ == "__main__":
    main()
