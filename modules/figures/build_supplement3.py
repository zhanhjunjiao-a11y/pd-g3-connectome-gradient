"""Render verified localization estimates only; no model fitting or ROI selection."""
from pathlib import Path
import os, sys, json, hashlib, platform

HERE = Path(os.environ['GRADIENT_FIG_OUTPUT'])/'supplement3_work'
HERE.mkdir(parents=True,exist_ok=True)
ROOT = Path(os.environ['GRADIENT_FIG_OUTPUT'])
# Use the declared environment; no historical runtime injection.
os.environ['MPLCONFIGDIR'] = str(HERE / 'mpl_config')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from PIL import Image
# PDF geometry is fixed by figsize; structural PDF verification is performed by the package verifier.

NETS = ['Visual', 'DorsAttn', 'VentAttn', 'Frontoparietal', 'Default', 'Limbic', 'Subcortex']
LABELS = ['Visual (29)', 'Dorsal attention (28)', 'Ventral attention (26)',
          'Frontoparietal (25)', 'Default mode (46)', 'Limbic (26)',
          'Subcortex /\nbrainstem (42)']
FILES = {'distance': 'network_distance_results.csv',
         'centered': 'network_centered_distance_results.csv',
         'raw_fc': 'network_raw_FC_results.csv',
         'gap': 'network_raw_gap_results.csv'}
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
DATA = Path(os.environ['GRADIENT_FIG_INPUT'])/'supplement3'
tables = {k: pd.read_csv(DATA / v) for k, v in FILES.items()}
for key, table in tables.items():
    assert set(table.network) == set(NETS)
    assert table.n.eq(475).all() and table.n_HC.eq(141).all() and table.n_PD.eq(334).all()
    assert table.standard_errors.eq('HC3').all()
    assert table.family_n.eq(14 if key == 'raw_fc' else 7).all()
    assert np.isfinite(table[['estimate', 'ci_low', 'ci_high', 'q_BH']]).all().all()
    assert (table.ci_low <= table.estimate).all() and (table.estimate <= table.ci_high).all()

OUT = ROOT / 'figures' / 'Supplementary_Figure_3'
OUT.parent.mkdir(parents=True, exist_ok=True)
BLUE, ORANGE, INK = '#0072B2', '#B34D00', '#243B53'
STYLE = {'font.family': 'Arial', 'font.size': 8.5, 'axes.labelsize': 8.5,
         'xtick.labelsize': 8, 'ytick.labelsize': 8.5, 'legend.fontsize': 8,
         'axes.spines.top': False, 'axes.spines.right': False,
         'axes.spines.left': False, 'axes.linewidth': .7,
         'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none',
         'savefig.facecolor': 'white', 'figure.facecolor': 'white'}
records = []

def forest(ax, rows, values='estimate', lo='ci_low', hi='ci_high',
           offset=0, color=INK, marker='o', label=None, stars=True, fill=True, panel=''):
    rows = rows.set_index('network').loc[NETS]
    x = rows[values].to_numpy()
    err = np.vstack([x-rows[lo].to_numpy(), rows[hi].to_numpy()-x])
    y = np.arange(7)+offset
    ax.errorbar(x, y, xerr=err, fmt=marker, ms=4.3, color=color,
                mfc=color if fill else 'white', mec=color, mew=.85,
                elinewidth=1.0, capsize=2, capthick=.85, label=label, zorder=3)
    for i, (net, row) in enumerate(rows.iterrows()):
        if stars and row.q_BH < .05:
            ax.annotate('*', (row[hi], y[i]), xytext=(3, 0), textcoords='offset points',
                        ha='left', va='center', color=color, fontsize=10)
        records.append({'panel': panel, 'network': net, 'series': label,
                        'estimate_column': values, 'estimate': float(row[values]),
                        'ci_low': float(row[lo]), 'ci_high': float(row[hi]),
                        'q_BH': float(row.q_BH), 'family_n': int(row.family_n),
                        'star_displayed': bool(stars and row.q_BH < .05)})

def base(ax, show_labels):
    ax.set_ylim(6.55, -.55)
    ax.set_yticks(np.arange(7), LABELS if show_labels else ['']*7)
    ax.tick_params(axis='y', length=0, pad=6)
    ax.tick_params(axis='x', length=3, pad=3)
    ax.axvline(0, color='#747474', lw=.7, ls=(0, (3, 3)), zorder=1)
    ax.grid(axis='x', color='#E3E7EB', linewidth=.45, zorder=0)
    for y in np.arange(7):
        if y % 2 == 0:
            ax.axhspan(y-.47, y+.47, color='#F6F7F9', lw=0, zorder=0)

def title(ax, letter, text, y=1.215):
    ax.text(0, y, letter, transform=ax.transAxes, ha='left', va='bottom',
            fontsize=11, fontweight='bold', color='#161616')
    ax.text(.13, y, text, transform=ax.transAxes, ha='left', va='bottom',
            fontsize=9.1, fontweight='bold', color='#161616', linespacing=1.15)

with plt.rc_context(STYLE):
    fig = plt.figure(figsize=(6.5, 7.1))
    axes = [fig.add_axes([.235, .600, .292, .285]),
            fig.add_axes([.653, .600, .310, .285]),
            fig.add_axes([.235, .140, .292, .285]),
            fig.add_axes([.653, .140, .310, .285])]
    a, b, c, d = axes
    for ax, labels in zip(axes, [True, False, True, False]): base(ax, labels)
    title(a, 'A', 'Changes in target\nconnectivity')
    title(b, 'B', 'Signed connectivity\ndifference between systems')
    title(c, 'C', 'Normalized profile\ndistance contributions', 1.15)
    title(d, 'D', 'Centered profile\ndistance contributions', 1.15)
    raw = tables['raw_fc']
    forest(a, raw[raw.seed.eq('SM')], offset=-.15, color=BLUE, marker='o', label='SMN', panel='A')
    forest(a, raw[raw.seed.eq('CB')], offset=.15, color=ORANGE, marker='s', label='Cerebellum', panel='A')
    a.set_xlim(-.087, .053); a.set_xticks([-.08, -.04, 0, .04])
    a.set_xlabel('Adjusted PD − HC FC\n(Fisher z)', labelpad=5)
    a.legend(loc='lower left', bbox_to_anchor=(-.05, 1.03), frameon=False,
             ncol=2, columnspacing=.8, handlelength=1.1, handletextpad=.35, borderaxespad=0)
    gap = tables['gap']
    forest(b, gap, 'HC_adjusted_mean', 'HC_adjusted_ci_low', 'HC_adjusted_ci_high',
           offset=-.15, color=BLUE, marker='o', label='HC', stars=False, fill=False, panel='B')
    forest(b, gap, 'PD_adjusted_mean', 'PD_adjusted_ci_low', 'PD_adjusted_ci_high',
           offset=.15, color=ORANGE, marker='s', label='PD', stars=False, panel='B')
    b.set_xlim(-.16, .275); b.set_xticks([-.1, 0, .1, .2])
    b.set_xlabel('Adjusted SMN − cerebellum FC\n(Fisher z)', labelpad=5)
    b.legend(loc='lower left', bbox_to_anchor=(-.035, 1.03), frameon=False,
             ncol=2, columnspacing=1.5, handlelength=1.1, handletextpad=.35, borderaxespad=0)
    for i, row in gap.set_index('network').loc[NETS].iterrows():
        if row.q_BH < .05:
            b.text(.262, NETS.index(i), '*', ha='center', va='center', fontsize=10, color=INK)
    for ax, key, panel in [(c, 'distance', 'C'), (d, 'centered', 'D')]:
        forest(ax, tables[key], marker='o' if key=='distance' else 'D', color=INK, panel=panel)
        ax.set_xlim(-.071, .014); ax.set_xticks([-.06, -.04, -.02, 0])
        ax.set_xlabel('Adjusted PD − HC distance\nNegative: greater profile similarity', labelpad=5)
    c.text(0, 1.04, 'Uncentered cosine distance', transform=c.transAxes, fontsize=8)
    d.text(0, 1.04, 'Sensitivity to profile mean offsets', transform=d.transAxes, fontsize=8)
    fig.text(.235, .035, '334 PD / 141 HC   ·   Points and bars: estimates and 95% CIs', fontsize=8)
    fig.text(.235, .012, '* BH q < 0.05; correction families are specified in the caption.', fontsize=8)
    fig.canvas.draw()
    # Verify that all text lies within the intended physical canvas.
    renderer = fig.canvas.get_renderer()
    boundary = fig.bbox
    clipped = []
    for txt in fig.findobj(match=matplotlib.text.Text):
        if not txt.get_visible() or not txt.get_text(): continue
        box = txt.get_window_extent(renderer)
        if box.x0 < -1 or box.y0 < -1 or box.x1 > boundary.width+1 or box.y1 > boundary.height+1:
            clipped.append(txt.get_text())
    assert not clipped, clipped
    for ext in ['png', 'pdf', 'svg']:
        fig.savefig(OUT.with_suffix('.'+ext), dpi=600, bbox_inches=None)
    preview = HERE/'Supplementary_Figure_3_preview.png'
    fig.savefig(preview, dpi=180, bbox_inches=None)
    plt.close(fig)

caption = '''Supplementary Figure 3. Localization of SMN–cerebellar connectivity-profile convergence to common target networks.
The atlas-defined somatomotor network (SMN; 30 regions) and cerebellum (31 regions) were represented by their mean functional connectivity (FC) to the same 222 external regions. All seven target categories are shown; numbers in parentheses indicate target-region counts. Direct SMN–cerebellar connections and within-system connections are excluded from these external profiles. (A) Adjusted PD–HC differences in mean FC from each seed system to each target category. Blue circles denote SMN and orange squares denote cerebellum. Positive values indicate higher FC in PD; for negative baseline FC, a positive difference may reflect weaker anticorrelation. (B) Adjusted HC and PD levels of the within-participant signed SMN-minus-cerebellum FC difference, retaining the baseline direction of the difference. Open blue circles denote HC and filled orange squares denote PD. Asterisks at the right edge indicate a significant PD–HC change in this signed difference, not a test of absolute convergence. (C) PD–HC differences in summed network contributions to cosine distance. For external profiles p and q, u=p/||p|| and v=q/||q||; each target contributes 0.5(u_j−v_j)^2, and the sum across all targets equals 1−cos(p,q). Negative group differences indicate smaller normalized profile distances in PD. Contributions are sums and therefore depend on target-category size; per-target averages are supplied in the accompanying source table. (D) The same decomposition after subtracting each profile’s mean across all 222 targets before normalization, a sensitivity analysis corresponding to Pearson correlation distance. Panels C and D use the same horizontal scale. All points and bars represent estimates and pointwise 95% confidence intervals from ordinary least-squares models adjusted for age, sex and mean framewise displacement, with HC3 standard errors (334 PD, 141 HC). Adjusted levels in B are evaluated at pooled sample mean covariate values. Asterisks denote Benjamini–Hochberg q<0.05 across 14 tests in A and across seven tests separately in B, C and D. All analyses are exploratory. Target contributions describe the normalized coordinate differences; they are neither independent causal effects nor an additive decomposition of angular similarity. Centering the profiles is not global signal regression.
'''
(HERE/'figure_caption.txt').write_text(caption, encoding='utf-8')
alt = 'Four forest-plot panels show seven shared target categories. Raw FC differences differ between the SMN and cerebellum; HC/PD signed system differences retain their baseline signs. Normalized and centered profile distances are lower in PD across six categories after correction; the default-mode category does not pass correction. Error bars represent 95% confidence intervals.'
(HERE/'Supplementary_Figure_3_alt_text.txt').write_text(alt, encoding='utf-8')
pd.DataFrame(records).to_csv(HERE/'Supplementary_Figure_3_plotted_values.csv', index=False, encoding='utf-8-sig', float_format='%.14g')
# Matplotlib generated one figure/page at the declared dimensions.
image = Image.open(OUT.with_suffix('.png'))
assert image.size == (3900, 4260)
assert OUT.with_suffix('.pdf').read_bytes().startswith(b'%PDF-')


assert image.getextrema()[3] == (255, 255)
report = {'figure': 'Supplementary Figure 3', 'width_in': 6.5, 'height_in': 7.1,
          'png_pixels': list(image.size), 'png_dpi': image.info.get('dpi'),
          'png_mode': image.mode, 'pdf_pages_expected': 1,
          'minimum_font_pt': 8, 'statistical_fitting': False,
          'all_seven_target_categories_shown': True, 'pointwise_ci_level': .95,
          'sources': {v: sha(DATA/v) for v in FILES.values()},
          'outputs': {ext: {'path': str(OUT.with_suffix('.'+ext)), 'sha256':sha(OUT.with_suffix('.'+ext))} for ext in ['png','pdf','svg']},
          'builder_sha256': sha(__file__), 'manual_visual_review': 'pending',
          'journal_compliance': 'Provisional manuscript supplement layout; not a publisher certification.',
          'python': platform.python_version(), 'matplotlib': matplotlib.__version__}
(HERE/'Supplementary_Figure_3_QA.json').write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding='utf-8')
print(json.dumps(report, ensure_ascii=False, indent=2))
