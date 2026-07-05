# -*- coding: utf-8 -*-
"""
取数 — 供组间特征值/ G3-loading 检验用。
Part A 立即可出(读现成 npy);Part B 需重跑(原 pipeline 丢了每被试 lambdas)。
在能看到实盘 E:\预测模型 的环境运行(本机 python 或 Claude Code)。仅需 numpy/pandas。
"""
import numpy as np, pandas as pd, os
from pathlib import Path

GRAD = Path(r"E:\预测模型\03_analysis\a_gradient\gradient_results")
OUT  = Path(r"E:\预测模型\results\eigentest_export"); OUT.mkdir(parents=True, exist_ok=True)

# ============================================================
# Part A — per-ROI HC-mean / PD-mean G3 loading(283 维各一列)
# gradients_{sess}.npy = (n_subj, 283, 10);列 index 2 = G3。
# 注意:这些 loading 先按特征值缩放(gradients[:,j]*=λ_j),再做 Procrustes(只旋转、
#       scale=False)→ 量级在对齐中被保留,故 HC-vs-PD 的 G3 均值差反映真实几何,不是被重缩放。
# ============================================================
hc      = np.load(GRAD / "gradients_HC.npy")     # (141,283,10) HC
pd_base = np.load(GRAD / "gradients_PD.npy")     # PD 基线
off     = np.load(GRAD / "gradients_ses1.npy")   # PD OFF（手册"主分析用 OFF 态"）
on      = np.load(GRAD / "gradients_ses2.npy")   # PD ON

df = pd.DataFrame({
    "roi_id":               np.arange(1, 284),
    "hc_mean_g3":           hc[:, :, 2].mean(0),
    "pd_mean_g3_baseline":  pd_base[:, :, 2].mean(0),
    "pd_mean_g3_OFF":       off[:, :, 2].mean(0),   # ← 与 OFF 态主分析对齐时用这列
    "pd_mean_g3_ON":        on[:, :, 2].mean(0),
})
df["diff_OFF_minus_HC"] = df["pd_mean_g3_OFF"] - df["hc_mean_g3"]
df.to_csv(OUT / "hcpd_g3_meanloading.csv", index=False)
# 也存成 npy 方便直接用
np.save(OUT / "hc_mean_g3.npy",  df["hc_mean_g3"].values)
np.save(OUT / "pd_mean_g3_OFF.npy", df["pd_mean_g3_OFF"].values)
print("Part A 完成:", OUT / "hcpd_g3_meanloading.csv", df.shape)
print("  ⚠️ 确认 t_map_adjusted 的 GLM 用的是 PD-baseline 还是 ses1-OFF 的 loading,取对应列做组间检验。")

# ============================================================
# Part B(可选,需重跑)— 每被试 λ1..λ10,分 HC / PD-OFF
# 为什么要重跑:原 pipeline `aligned, _ = process_single_subject(...)` 把每被试 lambdas
#             丢了,只存了对齐后梯度。故特征值不在盘上,需从 FC 矩阵重算。
# 嵌入与 grad_utils.py 完全一致(normalized-angle → 90% 稀疏化 → diffusion map alpha=0.5)。
# ============================================================
RUN_PART_B = True   # 想要每被试特征值就改成 True,并确认下面 FC_DIRS 路径
# [verified 2026-06-08] FC_DIRS + glob 'zsub-*.txt' match gradient_full_pipeline.py
# (GretnaSFCMatrixZ, Fisher-z); data lives under 02_data_raw\Session; sanity corr(HC-mean G1)=1.0
if RUN_PART_B:
    import glob
    def nak(fc):
        fc = fc.copy(); np.fill_diagonal(fc, 0)
        n = np.linalg.norm(fc, axis=1, keepdims=True); n[n == 0] = 1e-10
        c = np.clip((fc / n) @ (fc / n).T, -1, 1); return 1 - np.arccos(c) / np.pi
    def spar(a, s=0.9):
        a = a.copy()
        for i in range(a.shape[0]):
            r = a[i]; r[r < np.percentile(r, s * 100)] = 0; a[i] = r
        return (a + a.T) / 2
    def lambdas_of(fc, nc=10, alpha=0.5):
        a = spar(nak(fc))
        D = a.sum(1); Da = np.power(D, -alpha); Da[~np.isfinite(Da)] = 0
        M = a * np.outer(Da, Da); D2 = M.sum(1); D2i = 1 / D2; D2i[~np.isfinite(D2i)] = 0
        P = M * D2i[:, None]; s = np.sqrt(D2); si = 1 / s; si[~np.isfinite(si)] = 0
        Ps = P * s[:, None] * si[None, :]; Ps = (Ps + Ps.T) / 2
        w, _ = np.linalg.eigh(Ps); w = np.sort(w)[::-1]
        return w[1:nc + 1]                       # λ1..λ10(跳过平凡 λ0=1）
    FC_DIRS = {   # ← 确认指向 Fisher-z 的 283×283 FC 矩阵(原 pipeline 用 GretnaSFCMatrixZ)
        "HC":  r"E:\预测模型\02_data_raw\Session\ses-HC-rmCSFWM-Gaussian\GretnaSFCMatrixZ",
        "OFF": r"E:\预测模型\02_data_raw\Session\ses-1-rmCSFWM-Gaussian\GretnaSFCMatrixZ",
    }
    rows = []
    for grp, d in FC_DIRS.items():
        files = sorted(glob.glob(os.path.join(d, "zsub-*.txt")))
        print(f"  {grp}: {len(files)} 被试")
        for f in files:
            lam = lambdas_of(np.loadtxt(f))
            sid = os.path.basename(f).replace(".txt", "").replace("zsub-", "sub-")
            rows.append({"group": grp, "sub_id": sid,
                         **{f"lambda{k+1}": lam[k] for k in range(10)}})
    pd.DataFrame(rows).to_csv(OUT / "per_subject_eigenvalues.csv", index=False)
    print("Part B 完成:", OUT / "per_subject_eigenvalues.csv")
