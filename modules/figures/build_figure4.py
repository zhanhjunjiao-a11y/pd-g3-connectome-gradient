"""Rebuild Figure 4 from frozen model outputs; no model fitting or correction.

Run: python code/build_figure4.py
Dependencies: numpy, matplotlib, pillow (pypdf optional for export audit).
Input CSV files are bundled in ../source_data; previous figures are not touched.
"""
from pathlib import Path
import os
import sys
import json
import csv
import hashlib

ROOT = Path(os.environ['GRADIENT_FIG_OUTPUT'])
DATA = Path(os.environ['GRADIENT_FIG_INPUT'])
OUT = ROOT / 'figures'
QA = ROOT / 'qa' / 'figure4'
OUT.mkdir(parents=True, exist_ok=True)
QA.mkdir(parents=True, exist_ok=True)
os.environ.setdefault('MPLCONFIGDIR', str(QA / 'mplconfig'))
# Dependencies are supplied by the declared environment.

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from PIL import Image

SOURCE = DATA / 'Figure_4_clinical_associations.csv'
source_hash = hashlib.sha256(SOURCE.read_bytes()).hexdigest()
rows = list(csv.DictReader(SOURCE.open(encoding='utf-8-sig', newline='')))
metrics = ['fixed19_g3', 'compression', 'gradient_distance_3d', 'gradient_distance_10d',
           'profile_angle', 'within_sm_fc', 'within_cb_fc', 'between_fc']
labels = ['Fixed 19-region G3 score', 'Signed G3 difference\n(SMN − cerebellum)',
          '3D centroid distance', '10D centroid distance',
          'External-profile similarity', 'Within-somatomotor FC',
          'Within-cerebellar FC', 'Between-system FC']
lookup = {(r['question'], r['metric']): r for r in rows}
assert len(rows) == len(lookup) == 16
for question, n, df in [('severity', 334, 329), ('improvement', 235, 228)]:
    check_file = DATA / f'Figure_4_audit_{question}.csv'
    audit_rows = {r['metric']: r for r in csv.DictReader(check_file.open(encoding='utf-8-sig', newline=''))}
    for metric in metrics:
        r = lookup[question, metric]
        check = audit_rows[metric]
        assert int(r['n']) == int(check['n']) == n
        assert int(r['family_n']) == int(check['family_n']) == 8
        assert int(r['df']) == int(check['df']) == df
        assert r['standard_errors'] == check['standard_errors'] == 'HC3'
        for field, other in [('standardized_beta','standardized_beta'),
                             ('standardized_ci_low','standardized_ci95_low'),
                             ('standardized_ci_high','standardized_ci95_high'), ('q_BH8','q_BH8')]:
            assert float(r[field]) == float(check[other]), (question, metric, field)
        b, lo, hi, q = [float(r[k]) for k in ['standardized_beta', 'standardized_ci_low', 'standardized_ci_high', 'q_BH8']]
        assert np.isfinite([b, lo, hi, q]).all() and lo <= b <= hi and 0 <= q <= 1
assert all(float(lookup['improvement', metric]['q_BH8']) > .90 for metric in metrics)

matplotlib.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'DejaVu Sans'],
    'font.size': 8.0, 'axes.labelsize': 8.0, 'xtick.labelsize': 8.0,
    'ytick.labelsize': 8.0, 'legend.fontsize': 8.0,
    'axes.linewidth': .6, 'axes.spines.top': False, 'axes.spines.right': False,
    'axes.spines.left': False, 'xtick.major.width': .6, 'xtick.major.size': 2.5,
    'svg.fonttype': 'none', 'pdf.fonttype': 42,
    'figure.facecolor': 'white', 'savefig.facecolor': 'white'
})
BLUE, DARK, GREY = '#3B68A8', '#40464D', '#A8ADB2'
fig = plt.figure(figsize=(185 / 25.4, 142 / 25.4))
axes = [fig.add_axes([.31, .19, .285, .58]), fig.add_axes([.685, .19, .285, .58])]
fig.text(.03, .967, 'Clinical associations of gradient geometry and connectivity',
         fontsize=9.2, fontweight='bold', va='top')

fmt = lambda x: f'{x:+.3f}'.replace('-', '−')
plotted = []
for ax, question, panel, title, predictor, n in zip(
        axes, ['severity', 'improvement'], ['A', 'B'],
        ['OFF motor severity', 'Acute motor improvement'],
        ['OFF-state imaging', 'imaging change (ON − OFF)'], [334, 235]):
    x0, y0, width, height = ax.get_position().bounds
    fig.text(x0 - .035, .91, panel, fontsize=9.5, fontweight='bold', va='center')
    fig.text(x0, .91, title, fontsize=8.5, fontweight='bold', va='center')
    fig.text(x0, .875, 'Predictors:', fontsize=8.0, color=DARK, va='center')
    fig.text(x0, .844, predictor, fontsize=8.0, color=DARK, va='center')
    fig.text(x0, .808, f'n = {n}', fontsize=8.0, color=DARK, va='center')
    ax.vlines(np.zeros(8), np.arange(8)-.02, np.arange(8)+.42,
              color=GREY, linewidth=.7, linestyle=(0, (3, 3)), zorder=1)
    ax.axhline(3.5, color='#DDE0E3', linewidth=.6, zorder=0)
    for i, metric in enumerate(metrics):
        y = 7 - i
        r = lookup[question, metric]
        b, lo, hi, q = [float(r[k]) for k in ['standardized_beta', 'standardized_ci_low', 'standardized_ci_high', 'q_BH8']]
        ax.errorbar(b, y+.13, xerr=np.array([[b-lo], [hi-b]]), fmt='o',
                    markersize=4.4, markerfacecolor=BLUE if q < .05 else 'white',
                    markeredgecolor=BLUE, markeredgewidth=.9, ecolor=BLUE,
                    elinewidth=.85, capsize=2.0, capthick=.8, zorder=3)
        display = f'{fmt(b)}  [{fmt(lo)}, {fmt(hi)}]'
        ax.text(.5, y-.21, display, transform=ax.get_yaxis_transform(),
                ha='center', va='center', fontsize=8.0, color=DARK)
        plotted.append({'question': question, 'metric': metric, 'n': n,
                        'standardized_beta': b, 'standardized_ci_low': lo,
                        'standardized_ci_high': hi, 'q_BH8': q, 'display': display})
    ax.set_xlim(-.36, .36)
    ax.set_ylim(-.55, 7.6)
    ax.set_xticks([-.3, 0, .3], ['−0.3', '0', '0.3'])
    ax.set_yticks([])
    ax.set_xlabel('Standardized β (95% CI)', labelpad=5)
    if question == 'severity':
        for i, label in enumerate(labels):
            yf = fig.transFigure.inverted().transform(ax.transData.transform((0, 7-i)))[1]
            fig.text(.03, yf, label, va='center', fontsize=8.0, color=DARK, linespacing=1.05)

handles = [Line2D([], [], marker='o', linestyle='None', markersize=4.4,
                   markeredgecolor=BLUE, markerfacecolor=BLUE, label='q < 0.05'),
           Line2D([], [], marker='o', linestyle='None', markersize=4.4,
                   markeredgecolor=BLUE, markerfacecolor='white', label='q ≥ 0.05')]
fig.legend(handles=handles, loc='lower left', bbox_to_anchor=(.025, .075),
           ncol=2, frameon=False, handletextpad=.6, columnspacing=1.8, borderaxespad=0)
fig.text(.46, .091, 'BH correction: 8 tests per outcome', fontsize=8.0, va='center', color=DARK)
fig.text(.03, .037, 'Improvement: OFF − ON motor score. Imaging change: ON − OFF.', fontsize=8.0, color=DARK)

fig.canvas.draw()
renderer = fig.canvas.get_renderer()
outside = []
font_sizes = []
for text in fig.findobj(matplotlib.text.Text):
    if text.get_visible() and text.get_text():
        box = text.get_window_extent(renderer)
        font_sizes.append(text.get_fontsize())
        if box.x0 < -.5 or box.y0 < -.5 or box.x1 > fig.bbox.width+.5 or box.y1 > fig.bbox.height+.5:
            outside.append(text.get_text())
assert not outside, outside
assert min(font_sizes) >= 8
for ext in ['pdf','svg','png']:
    fig.savefig(OUT / f'Figure_4_revised.{ext}', dpi=600)
fig.savefig(OUT / 'Figure_4_preview.png', dpi=300)
plt.close(fig)

qa = {
    'source_sha256': source_hash,
    'source_values_crosschecked_against': ['Supplementary Data S09', 'Supplementary Data S10'],
    'source_files': {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in DATA.glob('Figure_4*.csv')},
    'rows_used': len(plotted), 'plotted_values': plotted, 'numeric_mapping_exact': True,
    'minimum_font_size_pt': min(font_sizes), 'text_outside_canvas': outside,
    'size_mm': [185, 142], 'raster_dpi': 600,
    'BH_families': {'severity': 8, 'improvement': 8},
    'model_refit': False, 'statistics_recomputed': False,
    'predictor_labels': ['Predictors: OFF-state imaging', 'Predictors: imaging change (ON − OFF)'],
    'signed_difference_label': 'Signed G3 difference (SMN − cerebellum)',
    'caption_changes_needed': 'No model changes. Explicitly identify A predictors as OFF-state imaging and B predictors as ON−OFF imaging change; signed G3 difference is SMN minus cerebellum.',
    'exports': {}
}
for ext in ['pdf','svg','png']:
    p = OUT / f'Figure_4_revised.{ext}'
    qa['exports'][ext] = {'bytes': p.stat().st_size, 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}
    if ext == 'png':
        with Image.open(p) as image:
            qa['exports'][ext].update({'pixels': list(image.size), 'dpi': list(image.info.get('dpi', ()))})
assert hashlib.sha256(SOURCE.read_bytes()).hexdigest() == source_hash
(QA / 'build_checks.json').write_text(json.dumps(qa, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps({'rows': len(plotted), 'size_mm': qa['size_mm'], 'minimum_font_size_pt': min(font_sizes), 'source_sha256': source_hash, 'outputs': list(qa['exports'])}))
