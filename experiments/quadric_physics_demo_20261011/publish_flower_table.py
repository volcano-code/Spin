"""Publish the single-file assembly with actual table simulation exports."""
import argparse
import hashlib
import json
from pathlib import Path
import zipfile


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def publish(bundle):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    relative_output = 'physics-job-02/output'
    out = bundle/relative_output
    result = json.loads((bundle/'physics-job-02/result.json').read_text())
    if result['status'] != 'PROCESS_EXITED' or result['returncode'] != 0:
        raise ValueError('TABLE_JOB_NOT_COMPLETE')
    summary = json.loads((out/'summary.json').read_text())
    render = json.loads((out/'render.json').read_text())
    if render['video_sha256'] != sha(out/'flower-table-comparison.mp4') or render['simulation_summary_sha256'] != sha(out/'summary.json'):
        raise ValueError('VIDEO_PHYSICS_IDENTITY_MISMATCH')
    if len(summary['scenes']) != 4:
        raise ValueError('MISSING_COMPARISON_RUN')
    for scene in summary['scenes']:
        path = out/scene['file']
        if sha(path) != scene['sha256']:
            raise ValueError('TRAJECTORY_CHANGED')
        with np.load(path) as saved:
            np.savetxt(path.with_suffix('.csv'), saved['records'], delimiter=',',
                       header=','.join(saved['columns']), comments='', fmt='%.17g')
    config = json.loads((bundle/'table-scene.json').read_text())
    if sha(bundle/'table-scene.json') != summary['source_config_sha256']:
        raise ValueError('TABLE_INPUT_CHANGED')
    for name, digest in config['files'].items():
        if sha(bundle/name) != digest:
            raise ValueError('MESH_CHANGED')
    checks = json.loads((bundle/'assembly-check.json').read_text())
    for name, check in checks.items():
        if sha(bundle/name) != check['sha256']:
            raise ValueError('WHOLE_ASSEMBLY_CHANGED')
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.6), constrained_layout=True)
    for name, color in [('original', '#bd547f'), ('optimized', '#b48a16')]:
        scene = next(row for row in summary['scenes'] if row['variant'] == name and row['dt'] == .00025)
        path = out/scene['file']
        if sha(path) != scene['sha256']:
            raise ValueError('TRAJECTORY_CHANGED')
        data = np.load(path)['records']
        mask = data[:, 0] <= .12
        axes[0].plot(data[mask, 0], data[mask, 3], label=name, color=color)
        axes[1].plot(data[::100, 0], data[::100, 2]/(2*np.pi), label=name, color=color)
    axes[0].set(xlabel='Time (s)', ylabel='Joint-frame reaction Fx (N)', title='Imbalance load')
    axes[1].set(xlabel='Time (s)', ylabel='Rotation rate (Hz)', ylim=(34.9, 35.1), title='Unpowered ideal bearing')
    for ax in axes:
        ax.grid(alpha=.2); ax.legend()
    fig.savefig(out/'comparison.png', dpi=160)
    plt.close(fig)
    assembly_names = ['flower-complete.obj', 'flower-rotor.obj', 'flower-cutaway.obj', 'flower.mtl',
                      'physical-properties.json', 'assembly-check.json', 'README.txt']
    with zipfile.ZipFile(bundle/'flower-complete-obj.zip', 'x', compression=zipfile.ZIP_DEFLATED) as archive:
        for name in assembly_names:
            archive.write(bundle/name, name)
    paths = [bundle/name for name in assembly_names+['table-scene.json', 'run_table_scene.py',
        'flower-table.urdf', 'urdf-reload-check.json', 'open_table_gui.py']]
    paths += sorted((bundle/'parts').glob('*.obj'))
    paths += [p for p in sorted(out.iterdir()) if p.suffix in ('.json', '.npz', '.csv', '.bullet', '.mp4', '.png')]
    with zipfile.ZipFile(bundle/'flower-pybullet-table.zip', 'x', compression=zipfile.ZIP_DEFLATED) as archive:
        for path in paths:
            archive.write(path, str(path.relative_to(bundle)))
    original, optimized = [config['properties'][name] for name in ('original', 'optimized')]
    e0, e1 = [np.linalg.norm(value['center_m'][:2])*1000 for value in (original, optimized)]
    forces = {name: next(row['mean_horizontal_reaction_N'] for row in summary['scenes']
        if row['variant'] == name and row['dt'] == .00025) for name in ('original', 'optimized')}
    html = f'''<!doctype html><html lang="zh"><meta charset="utf-8"><meta name="viewport" content="width=device-width">
<title>花盘完整 OBJ · PyBullet 桌面场景</title><style>
body{{font:16px/1.7 system-ui;color:#25313a;background:#f4f2ec;max-width:1100px;margin:auto;padding:26px}}
h1{{font-size:30px}}section{{background:white;padding:22px;margin:22px 0;border-radius:10px}}
img,video{{max-width:100%;height:auto}}a{{color:#235d84}}nav{{display:flex;gap:20px;flex-wrap:wrap}}table{{border-collapse:collapse;width:100%}}td,th{{padding:10px;text-align:left;border-bottom:1px solid #ddd}}
</style><h1>完整花盘与桌面物理仿真</h1>
<p>同一优化材料的额外交付 · 重建方法验证 · 米制单位</p>
<section><nav><a href="flower-complete-obj.zip" download>完整 OBJ 包（含颜色）</a>
<a href="flower-complete.obj" download>单文件完整装配 OBJ</a><a href="flower.mtl" download>配套 MTL</a>
<a href="flower-rotor.obj" download>单文件转子 OBJ</a><a href="flower-cutaway.obj" download>隐藏上盖的剖视 OBJ</a></nav>
<p>完整装配一文件包含9个部件，导入即显示全部结构。上下盖会遮住内部切口，剖视文件省略上盖供查看。轴承与转子保持各自网格，未焊成不可转动的实体。</p></section>
<section><h2>与论文一致的效果与当前范围</h2>
<p>论文第4.4、6.2节通过圆形切除调整质量分布，使质心与轴承中心重合。我们的二维圆边界／dogleg解、三层挤出和独立材料积分验证了这一效果。重建轮廓、固定区域和支撑条件与作者原始输入有差别，不能据此推断相同实物停转时间。</p>
<table><tr><th>独立保存材料结果</th><th>优化前</th><th>优化后</th></tr>
<tr><td>转子质量</td><td>{original['mass_kg']*1000:.6f} g</td><td>{optimized['mass_kg']*1000:.6f} g</td></tr>
<tr><td>质心偏离轴承中心</td><td>{e0:.6f} mm</td><td>{e1:.3g} mm</td></tr>
<tr><td>35 Hz 平均水平反力</td><td>{forces['original']:.6g} N</td><td>{forces['optimized']:.6g} N</td></tr></table>
<a href="https://pub.ista.ac.at/group_wojtan/projects/2024_Hafner_SIF/sif-final.pdf">原论文：Spin-It Faster</a></section>
<section><h2>桌面场景</h2><img src="{relative_output}/table-overview.png">
<p>两只花盘安装在桌上支座的理想轴承上，各自加载真实质量、质心和完整惯量。初始35 Hz，无持续驱动。桌子和支座固定，转子有净空；该场景对应轴承式指尖陀螺，并未模拟自由落桌碰撞。真实轴承摩擦尚未标定。</p></section>
<section><h2>实际 PyBullet 轨迹</h2><video controls preload="metadata" poster="{relative_output}/preview.png" src="{relative_output}/flower-table-comparison.mp4"></video>
<p>8秒、首相位、两个时间步（0.0005／0.00025秒）、优化前后共4次新增仿真。视频先正常速度，再以0.05倍回放前0.5秒。理想轴承中两者均持续转动，差异体现在偏心载荷。此前64场景的完整工程支撑矩阵仍在原案例页。</p>
<img src="{relative_output}/comparison.png"><nav><a href="{relative_output}/flower-table-comparison.mp4" download>下载视频</a><a href="flower-pybullet-table.zip" download>下载场景、脚本、轨迹与视频</a><a href="{relative_output}/summary.json">本次结果与两时间步对照</a></nav></section>
<section><p>重跑：解压仿真包，在目录执行 <code>python run_table_scene.py --bundle . --output my-table-result</code>。交互窗口执行 <code>python open_table_gui.py</code>；也提供 <a href="flower-table.urdf">URDF 场景</a>。CSV／NPZ 中的反力采用 Bullet 关节坐标系，未直接当作世界坐标分量。环境说明见 README。Bullet快照为诊断附件，可移植重建以脚本为准。</p><a href="../index.html">返回原优化结果及64场景验证</a></section></html>'''
    (bundle/'index.html').write_text(html, encoding='utf-8')
    manifest = dict(schema='flower-table-published-v1',
        files={str(p.relative_to(bundle)): sha(p) for p in paths},
        page_sha256=sha(bundle/'index.html'),
        archives={name: sha(bundle/name) for name in ('flower-complete-obj.zip', 'flower-pybullet-table.zip')})
    (bundle/'delivery-manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--bundle', type=Path, required=True)
    args = parser.parse_args()
    publish(args.bundle.resolve())
