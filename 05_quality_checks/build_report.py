#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""把机检结果机械写入审计报告的表格与附录，避免手抄。

用法：python3 build_report.py
读取同目录 check_manuscript_numbers.py 的 Audit 结果，替换报告里的 <!-- TABLE_B --> 占位。
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from check_manuscript_numbers import Audit, DEFAULT_PKG  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
REPORT = os.path.join(HERE, "握力与跌倒-STROBE逐条审计（证据锚定版）.md")


def fmt(x):
    if x is None:
        return "—"
    ax = abs(x)
    if ax != 0 and (ax < 0.001 or ax >= 100000):
        return f"{x:.3g}"
    return f"{x:g}"


def main():
    rows = Audit(DEFAULT_PKG).run()
    a = Audit(DEFAULT_PKG)
    rows = a.run()

    groups = {}
    for g, *_rest in rows:
        groups.setdefault(g, [0, 0])
        groups[g][0] += 1
    for g, _l, _c, _v, _cs, _asrc, ok in rows:
        if ok:
            groups[g][1] += 1

    out = []
    out.append("| 分组 | 核对项数 | 通过 | 结果 |")
    out.append("| --- | --- | --- | --- |")
    for g, (n, p) in groups.items():
        out.append(f"| {g} | {n} | {p} | {'全部通过' if n == p else f'{n - p} 项未通过'} |")
    out.append(f"| **合计** | **{len(rows)}** | **{sum(1 for r in rows if r[6])}** | "
               f"{'**全部通过**' if all(r[6] for r in rows) else '**存在未通过**'} |")
    out.append("")
    out.append("### 附录：机检逐项输出")
    out.append("")
    out.append("「正文值」为稿件中报告的数字，「文件值」为从 `04-关键结果数据/` 复算或读取的值。")
    out.append("")
    out.append("| 组 | 核对项 | 正文值 | 文件值 | 判定 |")
    out.append("| --- | --- | --- | --- | --- |")
    for g, label, claimed, actual, _cs, _asrc, ok in rows:
        out.append(f"| {g} | {label} | {fmt(claimed)} | {fmt(actual)} | "
                   f"{'通过' if ok else '**不一致**'} |")
    block = "\n".join(out)

    txt = open(REPORT, encoding="utf-8").read()

    # 幂等：先把上一次生成的正文块整体摘掉（正文块到下一个 "## " 之前），
    # 再重新按"摘要表留在第二节、逐项附录挪到文末"的版式拼回去。
    lines = txt.split("\n")
    kept, carry_on = [], False
    for ln in lines:
        if ln.startswith("| 分组 | 核对项数"):
            carry_on = True
            continue
        if carry_on and ln.startswith("## "):
            carry_on = False
        if not carry_on:
            kept.append(ln)
    txt = "\n".join(kept)

    summary, appendix = block.split("### 附录：机检逐项输出", 1)
    appendix = "### 附录：机检逐项输出" + appendix

    if "<!-- TABLE_B -->" in txt:
        txt = txt.replace("<!-- TABLE_B -->", summary.strip())
    else:
        anchor = "## 二、表 B：关键数字回归核对（机检 96 项）"
        idx = txt.find(anchor)
        nl = txt.find("\n", idx) + 1
        txt = txt[:nl] + "\n" + summary.strip() + "\n" + txt[nl:]

    txt = txt.rstrip() + "\n\n" + appendix.strip() + "\n"
    open(REPORT, "w", encoding="utf-8").write(txt)
    print(f"已写入 {len(rows)} 行机检结果 → {os.path.basename(REPORT)}")


if __name__ == "__main__":
    main()
