#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""稿件数字回归核对：把正文报告的每一个关键数字，回到结果文件里复算一遍。

用途：投稿前、改稿后、换口径后各跑一次，确认正文数字没有与分析结果脱钩。
设计原则：本脚本不判断"哪个数字对"，只判断"正文写的值是否等于结果文件里的值"。
          任何不一致都会以 FAIL 列出，由人决定改正文还是改分析。

用法：
  python3 check_manuscript_numbers.py [交接包目录]
  默认目录 = /Users/mac/Documents/Codex/2026-08-31/xian/outputs/妇科科研方向/握力与跌倒-投稿交接包

退出码：全部通过 0，存在 FAIL 为 1（可直接挂进投稿前的检查清单）。
"""
import csv
import math
import os
import sys

DEFAULT_PKG = ("/Users/mac/Documents/Codex/2026-08-31/xian/outputs/妇科科研方向/"
               "握力与跌倒-投稿交接包")


def read_csv(path):
    with open(path, encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))


def f(x):
    return float(str(x).strip().replace("%", ""))


def dl_pool(ors, los, his):
    """DerSimonian–Laird 随机效应合并（与正文 Methods 声明的估计量一致）。"""
    yi = [math.log(o) for o in ors]
    sei = [(math.log(h) - math.log(l)) / (2 * 1.959964) for l, h in zip(los, his)]
    wi = [1 / s ** 2 for s in sei]
    ybar = sum(w * y for w, y in zip(wi, yi)) / sum(wi)
    Q = sum(w * (y - ybar) ** 2 for w, y in zip(wi, yi))
    k = len(yi)
    c = sum(wi) - sum(w ** 2 for w in wi) / sum(wi)
    tau2 = max(0.0, (Q - (k - 1)) / c) if c > 0 else 0.0
    ws = [1 / (s ** 2 + tau2) for s in sei]
    pooled = math.exp(sum(w * y for w, y in zip(ws, yi)) / sum(ws))
    i2 = max(0.0, (Q - (k - 1)) / Q) * 100 if Q > 0 else 0.0
    return pooled, tau2, i2


class Audit:
    def __init__(self, pkg):
        self.pkg = pkg
        self.q = os.path.join(pkg, "04-关键结果数据", "队列分析")
        self.mr = os.path.join(pkg, "04-关键结果数据", "孟德尔随机化-FinnGen")
        self.rev = os.path.join(pkg, "04-关键结果数据", "反向孟德尔随机化")
        self.rows = []
        self.notes = []

    def check(self, group, label, claimed, claimed_src, actual, actual_src, tol=0.0005):
        ok = actual is not None and abs(claimed - actual) <= tol
        self.rows.append((group, label, claimed, actual, claimed_src, actual_src, ok))

    def note(self, text):
        self.notes.append(text)

    def run(self):
        # ---- 队列：样本与主结果 ----
        m08 = {r["结局"]: r for r in read_csv(os.path.join(self.q, "08_多队列合并.csv"))}
        m42 = {(r["结局"], r["SE口径"]): r
               for r in read_csv(os.path.join(self.q, "42_合并估计_设计校正.csv"))}
        a = m08["跌倒（任意）"]
        s = "08_多队列合并.csv"
        self.check("队列", "总样本量", 113846, "摘要/Results", f(a["总n"]), s, 0)
        self.check("队列", "跌倒事件数", 15872, "摘要/Results", f(a["总事件"]), s, 0)
        self.check("队列", "合并 OR", 0.862, "摘要/Results", f(a["合并OR"]), s)
        self.check("队列", "合并 OR 下限", 0.829, "摘要/Results", f(a["下限"]), s)
        self.check("队列", "合并 OR 上限", 0.897, "摘要/Results", f(a["上限"]), s)
        self.check("队列", "合并 P（1.1e-13）", 1.1e-13, "Results", f(a["p"]), s, 1e-15)
        self.check("队列", "I²（22.0%）", 22.0, "Results", f(a["I2_百分比"]), s, 0.05)

        b = m42[("跌倒（任意）", "survey 设计校正后")]
        self.check("队列", "τ²（0.00035）", 0.00035, "Results", f(b["tau2"]), "42_…设计校正.csv", 1e-5)
        c = m42[("跌倒（任意）", "校正前（仅加权）")]
        self.check("队列", "仅加权 OR（0.862）", 0.862, "Results", f(c["合并OR"]), "42_…设计校正.csv")
        self.check("队列", "仅加权 I²（29.0%）", 29.0, "Results", f(c["I2"]), "42_…设计校正.csv", 0.05)
        self.check("队列", "仅加权 τ²（0.00045）", 0.00045, "Results", f(c["tau2"]), "42_…设计校正.csv", 1e-5)

        inj = m08["跌倒受伤"]
        self.check("队列", "跌倒受伤 OR（0.905）", 0.905, "Results", f(inj["合并OR"]), s)
        self.check("队列", "跌倒受伤 下限（0.849）", 0.849, "Results", f(inj["下限"]), s)
        self.check("队列", "跌倒受伤 上限（0.965）", 0.965, "Results", f(inj["上限"]), s)
        self.check("队列", "跌倒受伤 I²=0", 0.0, "Results", f(inj["I2_百分比"]), s, 0.05)
        hip = m08["新发髋骨骨折"]
        self.check("队列", "髋骨骨折 OR（0.881）", 0.881, "Results", f(hip["合并OR"]), s)
        self.check("队列", "髋骨骨折 I²（20.4%）", 20.4, "Results", f(hip["I2_百分比"]), s, 0.05)

        # ---- 队列：逐队列、设计效应 ----
        m09 = read_csv(os.path.join(self.q, "09_分结局对照.csv"))
        want = {"CHARLS": 0.895, "ELSA": 0.896, "HRS": 0.858, "SHARE": 0.825}
        for r in m09:
            if r["结局"] != "跌倒（任意）":
                continue
            for k, v in want.items():
                if r["队列"].startswith(k):
                    self.check("逐队列", f"{k} OR", v, "Results", f(r["OR"]), "09_分结局对照.csv",
                               0.0015)

        deff = read_csv(os.path.join(self.q, "40_四队列设计效应汇总.csv"))
        want_deff = {"CHARLS": (1.42, 3327), "ELSA": (1.25, 7162), "HRS": (1.71, 7648),
                     "SHARE": (3.71, 23488)}
        seen = set()
        for r in deff:
            if r["结局"] != "跌倒（任意）":
                continue
            for k, (dv, nv) in want_deff.items():
                if r["队列"].startswith(k) and k not in seen:
                    seen.add(k)
                    self.check("设计效应", f"{k} DEFF", dv, "Results", f(r["DEFF_总体"]),
                               "40_…设计效应汇总.csv", 0.01)
                    self.check("设计效应", f"{k} 有效样本量", nv, "Results", f(r["有效样本量"]),
                               "40_…设计效应汇总.csv", 0.6)

        # ---- 队列：逐窗口、敏感性、阳性对照 ----
        win = read_csv(os.path.join(self.q, "51_逐窗口.csv"))
        ors = sorted(f(r["OR"]) for r in win)
        n_sig = sum(1 for r in win if f(r["p"]) < 0.05)
        mid = (ors[len(ors) // 2 - 1] + ors[len(ors) // 2]) / 2 if len(ors) % 2 == 0 \
            else ors[len(ors) // 2]
        self.check("稳健性", "窗口数（16）", 16, "Results", len(win), "51_逐窗口.csv", 0)
        self.check("稳健性", "窗口最小值（0.758）", 0.758, "Results", ors[0], "51_逐窗口.csv", 0.001)
        self.check("稳健性", "窗口最大值（0.954）", 0.954, "Results", ors[-1], "51_逐窗口.csv", 0.001)
        self.check("稳健性", "窗口中位数（0.862）", 0.862, "Results", mid, "51_逐窗口.csv", 0.001)
        self.check("稳健性", "单独显著窗口数（7）", 7, "Results", n_sig, "51_逐窗口.csv", 0)

        sens = read_csv(os.path.join(self.q, "52_关键敏感性.csv"))
        so = sorted(f(r["OR"]) for r in sens)
        self.check("稳健性", "敏感性分析数（11）", 11, "Results", len(sens), "52_关键敏感性.csv", 0)
        self.check("稳健性", "敏感性最小值（0.778）", 0.778, "Results", so[0], "52_关键敏感性.csv", 0.001)
        self.check("稳健性", "敏感性最大值（0.951）", 0.951, "Results", so[-1], "52_关键敏感性.csv", 0.001)

        pc = [r for r in read_csv(os.path.join(self.q, "50_阳性对照.csv"))
              if r["队列"].startswith("CHARLS") and "BMI" in r["对照"]]
        if pc:
            self.check("阳性对照", "CHARLS BMI OR（1.000）", 1.000, "Results", f(pc[0]["OR"]),
                       "50_阳性对照.csv", 0.001)
            self.check("阳性对照", "CHARLS BMI P（0.994）", 0.994, "Results", f(pc[0]["p"]),
                       "50_阳性对照.csv", 0.001)

        # ---- 队列：竞争风险、死亡构成 ----
        cr = read_csv(os.path.join(self.q, "72_竞争风险_合并.csv"))
        want_cr = [0.833, 0.857, 0.808, 0.892, 0.861]
        self.check("竞争风险", "口径数（5）", 5, "Results", len(cr), "72_…合并.csv", 0)
        for w, r in zip(want_cr, cr):
            self.check("竞争风险", f"{r['口径'][:1]} 口径 OR", w, "Results", f(r["合并OR"]),
                       "72_…合并.csv", 0.001)

        comp = read_csv(os.path.join(self.q, "71_竞争风险_人数构成.csv"))
        deaths = {r["队列"].split()[0]: f(r["结局前死亡"]) for r in comp}
        elig = {r["队列"].split()[0]: f(r["符合条件"]) for r in comp}
        total = sum(deaths.values())
        for k, (d, n) in {"CHARLS": (199, 7658), "ELSA": (200, 15762), "HRS": (841, 25049)}.items():
            self.check("死亡构成", f"{k} 结局前死亡", d, "Results", deaths.get(k),
                       "71_…人数构成.csv", 0)
            self.check("死亡构成", f"{k} 符合条件", n, "Results", elig.get(k),
                       "71_…人数构成.csv", 0)
        self.check("死亡构成", "死亡合计（4743）", 4743, "Results", total, "71_…人数构成.csv", 0)

        # ---- 队列：剂量反应 ----
        rcs = read_csv(os.path.join(self.q, "84_RCS_非线性检验.csv"))
        pooled = [r for r in rcs if "合并" in r["队列"] and "z" in r["尺度"]]
        if pooled:
            self.check("剂量反应", "合并非线性 P（0.0047）", 0.0047, "Results",
                       f(pooled[0]["P非线性"]), "84_RCS_非线性检验.csv", 1e-5)
        els = [r for r in rcs if r["队列"].startswith("ELSA") and "kg" in r["尺度"]]
        if els:
            self.check("剂量反应", "ELSA kg 尺度 P（0.011）", 0.011, "Results",
                       f(els[0]["P非线性"]), "84_RCS_非线性检验.csv", 5e-4)

        thr = read_csv(os.path.join(self.q, "81_RCS_上穿1的握力.csv"))
        want_thr = {"CHARLS": 15.5, "ELSA": 14.2, "HRS": 13.3, "SHARE": 18.6}
        for r in thr:
            for k, v in want_thr.items():
                if r["队列"].startswith(k):
                    self.check("剂量反应", f"{k} 阈值 kg", v, "Results",
                               f(r["主阈值kg_OR达1.25"]), "81_RCS_上穿1的握力.csv", 0.05)

        # ---- 队列：协变量口径分解（本地复算 DL 合并，正文 0.865→0.834）----
        cov = read_csv(os.path.join(self.q, "64_协变量口径分解.csv"))
        prim = [r for r in cov if r["模型"] == "主模型（本队列协变量）"]
        core = [r for r in cov if r["模型"] == "统一核心协变量"]
        if len(prim) == 4 and len(core) == 4:
            p1, _, i1 = dl_pool([f(r["OR"]) for r in prim], [f(r["下限"]) for r in prim],
                                [f(r["上限"]) for r in prim])
            p2, _, i2 = dl_pool([f(r["OR"]) for r in core], [f(r["下限"]) for r in core],
                                [f(r["上限"]) for r in core])
            self.check("协变量分解", "主模型合并 OR（0.865）", 0.865, "Results", p1,
                       "64（本脚本 DL 复算）", 0.001)
            self.check("协变量分解", "统一核心合并 OR（0.834）", 0.834, "Results", p2,
                       "64（本脚本 DL 复算）", 0.001)
            self.check("协变量分解", "主模型 I²（30.7%）", 30.7, "Results", i1,
                       "64（本脚本 DL 复算）", 0.5)
            self.check("协变量分解", "统一核心 I²（79.8%）", 79.8, "Results", i2,
                       "64（本脚本 DL 复算）", 0.5)

        # ---- MR ----
        summ = {f(r["序号"]): r for r in read_csv(os.path.join(self.mr, "09_汇总表.csv"))}
        lbl = {2: ("重叠 左手→UKB 跌倒", 0.940, 0.918, 0.963),
               3: ("重叠 左手→falling risk", 0.826, 0.756, 0.903),
               4: ("独立 左手→FinnGen", 0.971, 0.873, 1.080),
               5: ("独立 右手→FinnGen", 1.035, 0.937, 1.144),
               6: ("独立 低握力→FinnGen", 1.074, 1.003, 1.150),
               7: ("阳性对照 BMI", 1.119, 1.077, 1.162),
               8: ("阴性对照 LDL", 1.009, 0.983, 1.035)}
        for k, (name, o, lo, hi) in lbl.items():
            r = summ[k]
            self.check("MR", f"{name} OR", o, "Results", f(r["OR"]), "09_汇总表.csv", 0.001)
            self.check("MR", f"{name} 下限", lo, "Results", f(r["下限"]), "09_汇总表.csv", 0.001)
            self.check("MR", f"{name} 上限", hi, "Results", f(r["上限"]), "09_汇总表.csv", 0.001)
        w_over = f(summ[2]["上限"]) - f(summ[2]["下限"])
        w_ind = f(summ[4]["上限"]) - f(summ[4]["下限"])
        wl = math.log(f(summ[2]["上限"])) - math.log(f(summ[2]["下限"]))
        wi = math.log(f(summ[4]["上限"])) - math.log(f(summ[4]["下限"]))
        # 正文的 0.048 / 0.213 与比值 4.4 是对数 OR 尺度上的宽度
        self.check("MR", "CI 宽度 对数尺度 重叠（0.048）", 0.048, "Results", wl,
                   "09_汇总表.csv（复算）", 0.0006)
        self.check("MR", "CI 宽度 对数尺度 独立（0.213）", 0.213, "Results", wi,
                   "09_汇总表.csv（复算）", 0.0006)
        self.check("MR", "CI 宽度比 对数尺度（4.4）", 4.4, "Results", wi / wl,
                   "09_汇总表.csv（复算）", 0.05)
        self.note(f"CI 宽度 0.048／0.213 与比值 4.4 是对数 OR 尺度（{wl:.4f}／{wi:.4f}）。"
                  f"同一比较在 OR 尺度上是 {w_over:.4f}／{w_ind:.4f}，比值 {w_ind / w_over:.2f}；"
                  f"正文与图 4 面板 B 均未标注尺度，审稿人按报告的四舍五入 OR 复算会得到约 4.6 而非 4.4。")
        self.check("MR", "重叠点估计低于独立下限的差值（0.011）", 0.011, "Results",
                   f(summ[4]["下限"]) - 0.862, "09_汇总表.csv（复算）", 0.001)

        mr1 = read_csv(os.path.join(self.mr, "01_FinnGen_MR_主结果.csv"))
        nsnp = {r["暴露"]: f(r["n_snp"]) for r in mr1
                if r["方法"] == "Inverse variance weighted"}
        for key, v in [("左手", 145), ("右手", 161)]:
            hit = [n for k, n in nsnp.items() if key in k]
            if hit:
                self.check("MR", f"{key}工具变量数", v, "Methods", hit[0],
                           "01_FinnGen_MR_主结果.csv", 0)

        # ---- 反向 MR ----
        rev = read_csv(os.path.join(self.rev, "01_反向MR_主结果.csv"))
        ivw = [r for r in rev if r["方法"] == "Inverse variance weighted"]
        b_right = [r for r in ivw if "右手" in r["结局"]]
        b_left = [r for r in ivw if "左手" in r["结局"]]
        if b_right:
            self.check("反向MR", "右手 b（−0.016）", -0.016, "Results", f(b_right[0]["b"]),
                       "01_反向MR_主结果.csv", 0.001)
            self.check("反向MR", "右手 P（0.61）", 0.61, "Results", f(b_right[0]["pval"]),
                       "01_反向MR_主结果.csv", 0.01)
        if b_left:
            self.check("反向MR", "左手 b（−0.021）", -0.021, "Results", f(b_left[0]["b"]),
                       "01_反向MR_主结果.csv", 0.001)
            self.check("反向MR", "左手 P（0.48）", 0.48, "Results", f(b_left[0]["pval"]),
                       "01_反向MR_主结果.csv", 0.01)

        # ---- 工具变量强度：从原始工具变量文件现算 F ----
        fmin = None
        for fn in ("inst_gripL.csv", "inst_gripR.csv", "inst_lowGrip.csv"):
            path = os.path.join(self.mr, fn)
            if not os.path.exists(path):
                continue
            fs = [(f(r["beta.exposure"]) / f(r["se.exposure"])) ** 2 for r in read_csv(path)]
            if fs:
                fmin = min(fs) if fmin is None else min(fmin, min(fs))
        if fmin is not None:
            self.check("MR", "三个暴露集最小 F（29.6）", 29.6, "Methods", fmin,
                       "inst_*.csv（现算 F=beta²/se²）", 0.1)
            self.note(f"三个暴露集最小 F 的精确值为 {fmin:.2f}（低握力集），正文写 29.6 属四舍五入；"
                      f"全部工具变量 F 均 > 10。")

        pw = {r["项"]: f(r["值"]) for r in read_csv(os.path.join(self.rev, "17_功效换算.csv"))}
        self.check("反向MR", "80% 功效 MDE（0.088 SD）", 0.088, "Results",
                   pw.get("80%功效_MDE_SD"), "17_功效换算.csv", 0.001)
        if pw.get("80%功效_MDE_SD"):
            self.check("反向MR", "MDE/隐含效应 比值（3.05）", 3.05, "Results",
                       pw["80%功效_MDE_SD"] / 0.029, "17_功效换算.csv（复算）", 0.05)
        return self.rows


def main():
    pkg = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_PKG
    a = Audit(pkg)
    rows = a.run()
    width = max(len(r[1]) for r in rows) + 2
    print(f"{'组':<8}{'核对项':<{width}}{'正文':>12}{'结果文件':>12}  判定")
    print("-" * (width + 44))
    fails = 0
    for g, label, claimed, actual, csrc, asrc, ok in rows:
        shown = "None" if actual is None else f"{actual:.6g}"
        flag = "PASS" if ok else "FAIL"
        if not ok:
            fails += 1
        print(f"{g:<8}{label:<{width}}{claimed:>12.6g}{shown:>12}  {flag}")
    print("-" * (width + 44))
    print(f"合计 {len(rows)} 项，通过 {len(rows) - fails} 项，未通过 {fails} 项")
    if a.notes:
        print("\n附注（不影响判定，投稿前值得处理）：")
        for n in a.notes:
            print(f"  · {n}")
    if fails:
        print("\n未通过明细（正文值 vs 结果文件值）：")
        for g, label, claimed, actual, csrc, asrc, ok in rows:
            if not ok:
                print(f"  [{g}] {label}：正文 {claimed:g}｜文件 {actual}｜"
                      f"正文来源 {csrc}｜文件来源 {asrc}")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
