"""独立复核：用与搜索程序完全无关的纯 Python/NumPy 实现，重新计算候选和剑桥纪录的能量。

- 不调用 thomson.py / energy.c（避免“自己验证自己”）；
- 两两距离用 NumPy float64 计算，求和用 math.fsum（精确舍入求和），总误差约 1e-10 量级；
- 同时检查：点数、是否都在单位球面上（先归一化并报告原始偏差）、有无重合点、最小间距；
- 识别对称性：检测 C2–C6 旋转轴、反演中心、镜面（容差 1e-6）。

用法：
  python3 src/verify_big.py --n 2502 --candidate runs_big/probe_2502_r0/CANDIDATE_xxx/points.xyz
  （剑桥纪录默认读 seeds/cambridge/points_<N>_1.dat；也可用 --record 指定）
输出：终端报告 + 候选目录下 verification.json
"""
import argparse, json, math, os, sys, time
import numpy as np
from scipy.spatial import cKDTree

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DB = {1382: 926576.1030330, 2000: 1950575.5888846, 2192: 2345720.9932072, 2502: 3060842.2403807,
      3182: 4963364.3188600, 3552: 6191355.6866585, 4352: 9311272.9909708}


def read_points(path):
    if path.endswith(".npy"):
        P = np.load(path)
    else:
        rows = []
        for line in open(path):
            parts = line.split()
            if len(parts) >= 3:
                try:
                    rows.append([float(v) for v in parts[:3]])
                except ValueError:
                    pass
        P = np.array(rows)
    return np.asarray(P, dtype=np.float64)


def energy_fsum(P):
    n = len(P)
    terms = []
    for i in range(n - 1):
        d = np.sqrt(((P[i + 1:] - P[i]) ** 2).sum(axis=1))
        terms.append(math.fsum((1.0 / d).tolist()))
    return math.fsum(terms)


def rot(axis, ang):
    a = axis / np.linalg.norm(axis)
    K = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
    return np.eye(3) + np.sin(ang) * K + (1 - np.cos(ang)) * K @ K


def symmetry(P, tol=1e-6):
    tree = cKDTree(P)
    ok = lambda M: tree.query(P @ M.T)[0].max() < tol
    found = {"inversion": bool(ok(-np.eye(3))), "rotations": [], "mirrors": 0}
    cands = [P[0] + P[j] for j in range(1, len(P))] + [np.cross(P[0], P[j]) for j in range(1, 40)]
    _, vecs = np.linalg.eigh(P.T @ P)
    cands += list(vecs.T)
    def refine(a, order):
        # 近似轴 → 用轨道平均精确求轴（特征向量本身只精确到约 1e-5）
        for _ in range(3):
            d, j = tree.query(P @ rot(a, 2 * np.pi / order).T)
            if d.max() > 1e-2:
                return a
            s, idx = P.copy(), np.arange(len(P))
            for _k in range(order - 1):
                idx = j[idx]
                s = s + P[idx]
            s = s[np.linalg.norm(s, axis=1) > 0.3]
            if not len(s):
                return a
            s *= np.sign(s @ a)[:, None]
            a = s.sum(axis=0) / np.linalg.norm(s.sum(axis=0))
        return a

    axes = []
    for a in cands:
        if np.linalg.norm(a) < 1e-8:
            continue
        a = a / np.linalg.norm(a)
        for order in (6, 5, 4, 3, 2):
            b = refine(a, order)
            if any(abs(abs(b @ c) - 1) < 1e-6 for c, _ in axes):
                continue
            if ok(rot(b, 2 * np.pi / order)):
                axes.append((b, order))
                break
    found["rotations"] = sorted([o for _, o in axes], reverse=True)
    normals = []
    for j in range(1, len(P)):
        nrm = P[0] - P[j]
        nrm /= np.linalg.norm(nrm)
        if any(abs(abs(nrm @ m) - 1) < 1e-6 for m in normals):
            continue
        if ok(np.eye(3) - 2 * np.outer(nrm, nrm)):
            normals.append(nrm)
    found["mirrors"] = len(normals)
    return found


def check(P, label):
    raw_dev = float(np.abs(np.linalg.norm(P, axis=1) - 1).max())
    Q = P / np.linalg.norm(P, axis=1, keepdims=True)
    dmin = float(cKDTree(Q).query(Q, k=2)[0][:, 1].min())
    t = time.time()
    e = energy_fsum(Q)
    print(f"{label}: N={len(Q)}  能量 {e:.10f}  （单位球偏差 {raw_dev:.1e}，最小间距 {dmin:.4f}，用时 {time.time() - t:.0f}s）",
          flush=True)
    return Q, e, raw_dev, dmin


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--candidate", required=True)
    ap.add_argument("--record", default="")
    a = ap.parse_args()
    rec_path = a.record or os.path.join(ROOT, "seeds", "cambridge", f"points_{a.n}_1.dat")
    C = read_points(a.candidate)
    R = read_points(rec_path)
    if len(C) != a.n or len(R) != a.n:
        sys.exit(f"点数不对：候选 {len(C)}，纪录 {len(R)}，应为 {a.n}")
    Qc, ec, devc, dminc = check(C, "候选")
    Qr, er, _, _ = check(R, "剑桥纪录")
    db = DB.get(a.n)
    imp = er - ec
    print(f"\n剑桥数据库公布值 {db}；本程序复算纪录 {er:.10f}" + (f"（差 {er - db:+.2e}）" if db else ""))
    print(f"改进量（纪录 − 候选）= {imp:.10f}")
    sym = symmetry(Qc)
    print(f"候选对称性：旋转轴阶数 {sym['rotations'] or '无'}；反演中心 {'有' if sym['inversion'] else '无'}；镜面 {sym['mirrors']} 个")
    verdict = imp > 1e-5 and dminc > 1e-3
    print("\n结论：" + ("✅ 候选能量确实更低（改进量远大于计算误差）" if verdict else "❌ 未确认改进，不要对外宣称"))
    out = dict(N=a.n, candidate_energy=ec, record_energy_recomputed=er, record_energy_database=db,
               improvement=imp, unit_sphere_max_deviation=devc, min_distance=dminc, symmetry=sym,
               method="independent NumPy float64 distances + math.fsum exact summation",
               verified=bool(verdict), time=time.strftime("%Y-%m-%dT%H:%M:%S"))
    target = os.path.dirname(os.path.abspath(a.candidate))
    json.dump(out, open(os.path.join(target, "verification.json"), "w"), ensure_ascii=False, indent=1)
    print(f"已写入 {os.path.join(target, 'verification.json')}")


if __name__ == "__main__":
    main()
