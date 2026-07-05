# Compute G3 (component index 2) means at the two pole ROI sets.
# Run with:  python compute_g3_pole_means.py
import numpy as np
from pathlib import Path

HERE = Path(__file__).parent
HC_PATH = HERE / "03_analysis" / "a_gradient" / "gradient_results" / "gradients_HC.npy"
PD_PATH = HERE / "03_analysis" / "a_gradient" / "gradient_results" / "gradients_PD.npy"

CORTICAL_POLE = [57, 58, 59, 60, 67, 68, 161, 162, 164, 271]   # 1-based labels
CEREBELLAR_POLE = [225, 257, 258, 259, 260, 261, 264, 272, 275]

def report(name, path):
    if not path.exists():
        print(f"[skip] {name}: missing {path}")
        return
    g = np.load(path)                       # shape (n_sub, 283, n_comp)
    print(f"\n=== {name}  shape={g.shape} ===")
    g3 = g[:, :, 2]                         # (n_sub, 283)
    cort_idx = [r - 1 for r in CORTICAL_POLE]
    cb_idx   = [r - 1 for r in CEREBELLAR_POLE]
    # group-mean topography: average across subjects first
    g3_grp = g3.mean(axis=0)                # (283,)
    cort_mean = g3_grp[cort_idx].mean()
    cb_mean   = g3_grp[cb_idx].mean()
    print(f"G3 mean at cortical pole  (10 ROIs): {cort_mean:+.6f}")
    print(f"G3 mean at cerebellar pole (9 ROIs): {cb_mean:+.6f}")
    print(f"diff (cortical - cerebellar)       : {cort_mean - cb_mean:+.6f}")
    print(f"abs span                            : {abs(cort_mean - cb_mean):.6f}")
    # global extremes for context
    print(f"G3 group-mean min={g3_grp.min():+.4f}  max={g3_grp.max():+.4f}")
    print(f"argmin label (1-based) = {int(g3_grp.argmin())+1}")
    print(f"argmax label (1-based) = {int(g3_grp.argmax())+1}")
    # per-subject pole-difference score (could serve as 'compression score')
    cort_sub = g3[:, cort_idx].mean(axis=1)
    cb_sub   = g3[:, cb_idx].mean(axis=1)
    pole_diff = cort_sub - cb_sub
    print(f"per-subject (cort - cereb) mean = {pole_diff.mean():+.6f}  std = {pole_diff.std(ddof=1):.6f}")

report("HC", HC_PATH)
report("PD", PD_PATH)

# Optional: HC vs PD group test on the per-subject pole-difference score
try:
    g_hc = np.load(HC_PATH)[:, :, 2]
    g_pd = np.load(PD_PATH)[:, :, 2]
    cort_idx = [r - 1 for r in CORTICAL_POLE]
    cb_idx   = [r - 1 for r in CEREBELLAR_POLE]
    pd_diff = g_pd[:, cort_idx].mean(1) - g_pd[:, cb_idx].mean(1)
    hc_diff = g_hc[:, cort_idx].mean(1) - g_hc[:, cb_idx].mean(1)
    from scipy import stats
    t, p = stats.ttest_ind(hc_diff, pd_diff, equal_var=False)
    pooled_sd = np.sqrt(((len(hc_diff)-1)*hc_diff.var(ddof=1) + (len(pd_diff)-1)*pd_diff.var(ddof=1)) / (len(hc_diff)+len(pd_diff)-2))
    d = (hc_diff.mean() - pd_diff.mean()) / pooled_sd
    print("\n=== HC vs PD on (cortical - cerebellar) pole-difference ===")
    print(f"HC mean = {hc_diff.mean():+.6f}  (n = {len(hc_diff)})")
    print(f"PD mean = {pd_diff.mean():+.6f}  (n = {len(pd_diff)})")
    print(f"Welch t = {t:.3f}   p = {p:.3e}   Cohen's d = {d:+.3f}")
    print("If sign of HC mean is positive, cortical pole is G3 positive end.")
    print("If |PD diff| < |HC diff|, 'compression' interpretation supported.")
except Exception as e:
    print("PD/HC group test skipped:", e)
