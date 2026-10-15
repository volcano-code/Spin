"""Publish only the authenticated current material, archives and media."""
import argparse
import hashlib
import json
from pathlib import Path


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--case',type=Path,required=True)
    args=parser.parse_args();case=args.case.resolve()
    delivery=case/'delivery-job-03/output'
    archive=json.loads((delivery/'delivery-manifest.json').read_text())
    proof=case/'saved-independent-job-03/output/saved-independent-check.json'
    if archive['saved_check_sha256']!=sha(proof):raise ValueError('ARCHIVE_PROOF_MISMATCH')
    for name,digest in archive['archives'].items():
        if sha(delivery/name)!=digest:raise ValueError('ARCHIVE_CHANGED')
    studio=case/'qualified-studio-02/render-job-01/output'
    images=json.loads((studio/'render-manifest.json').read_text())
    for name,digest in images['images'].items():
        if sha(studio/(name+'.png'))!=digest:raise ValueError('STUDIO_IMAGE_CHANGED')
    display=case/'paper-display-job-03/output'
    figure=json.loads((display/'display-manifest.json').read_text())
    if figure['delivery_manifest_sha256']!=sha(delivery/'delivery-manifest.json') or figure['video_sha256']!=sha(display/'paper-physical-comparison.mp4'):
        raise ValueError('MEDIA_DELIVERY_IDENTITY_MISMATCH')
    physics=json.loads((case/'physics-job-02/output/summary.json').read_text())
    if physics['recorded_scenes']!=64:raise ValueError('PHYSICS_MATRIX_INCOMPLETE')
    snapshots=[('input','输入：透明固定层与粉色设计域'),('cutaway','优化中间层：黄色保留材料'),('exploded','三层展开：右下侧切除清楚可见'),('assembly','完整装配与夹持件')]
    cards=''.join('<figure><img src="qualified-studio-02/render-job-01/output/'+name+'.png"><figcaption>'+title+'</figcaption></figure>' for name,title in snapshots)
    html='''<!doctype html><html lang="zh"><meta charset="utf-8"><meta name="viewport" content="width=device-width">
    <title>Spin-It Faster 花盘 · 已验证优化结果</title><style>
    body{font:16px/1.7 system-ui;color:#252a30;background:#f5f4ef;max-width:1200px;padding:28px;margin:auto}
    h1{font-size:30px}section{background:white;padding:24px;margin:24px 0;border-radius:12px}
    .gallery{display:grid;grid-template-columns:1fr 1fr;gap:20px}figure{margin:0}img,video{width:100%;height:auto}
    figcaption{padding:12px 0}a{color:#235b80}nav{display:flex;gap:22px;flex-wrap:wrap}.note{color:#606872}
    @media(max-width:750px){.gallery{grid-template-columns:1fr}}</style>
    <h1>花盘：真实优化、最终 OBJ 与三维物理对比</h1>
    <p>Spin-It Faster 二维圆边界／dogleg 方法，实际三层结构挤出。重建方法验证（RECONSTRUCTED_METHOD_VALIDATION）。</p>
    <section><nav><a href="delivery-job-03/output/final-assembly-obj.zip" download>下载最终装配 OBJ 包</a>
    <a href="delivery-job-03/output/final-assembly/retained.obj" download>单独下载优化中间层</a>
    <a href="paper-display-job-03/output/index.html">交互三维／二维详解</a>
    <a href="delivery-job-03/output/diagnostic-obj.zip" download>原始／设计域／空腔诊断包</a></nav>
    <p>最终装配包只含9个实际部件与物理属性。查看切除位置时打开 retained.obj，或在完整装配中隐藏上下盖。</p></section>
    <section class="gallery">'''+cards+'''</section>
    <section><h2>中间层固定正面</h2><img src="paper-display-job-03/output/fixed-front.png">
    <p>圆心约在轴承中心右侧2.10 mm、下侧5.39 mm，半径5.78 mm；位置来自实际优化。</p></section>
    <section><h2>真实物理轨迹</h2><video controls preload="metadata" src="paper-display-job-03/output/paper-physical-comparison.mp4"></video>
    <p>35 Hz、8秒、八个初始相位、两个时间步，优化前后共64个场景。理想轴承与柔性支撑共用登记参数；视频包含正常速度与0.1倍慢放。</p>
    <p class="note">柔性支撑为10 Hz、阻尼比0.1的工程模型；具体数值与时间步差异见完整统计。制造资格尚未授予。</p>
    <a href="physics-job-02/output/summary.json">全部物理统计与参数</a> ·
    <a href="delivery-job-03/output/final-assembly/physical-properties.json">质量、质心、完整惯量与材料密度</a> ·
    <a href="saved-independent-job-03/output/saved-independent-check.json">保存后独立材料检查</a></section></html>'''
    target=case/'index.html'
    if target.exists():raise ValueError('CURRENT_RESULT_PAGE_ALREADY_EXISTS')
    target.write_text(html,encoding='utf-8')
    (case/'result-page-manifest.json').write_text(json.dumps(dict(schema='qualified-flower-result-page-v1',
        page_sha256=sha(target),delivery_sha256=sha(delivery/'delivery-manifest.json'),
        studio_sha256=sha(studio/'render-manifest.json'),display_sha256=sha(display/'display-manifest.json'),
        physics_sha256=sha(case/'physics-job-02/output/summary.json')),indent=2)+'\n')


if __name__=='__main__':main()
