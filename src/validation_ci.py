import math

import pandas as pd

d = pd.read_csv("data/processed/validation_results.csv")


def wilson(k, n, z=1.96):
    p = k / n
    den = 1 + z * z / n
    c = p + z * z / (2 * n)
    a = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return (c - a) / den, (c + a) / den


d["bin"] = pd.cut(d["max_sim"], [0, 0.4, 0.7, 1.01], labels=["<0.4", "0.4-0.7", ">=0.7"])
lines = []
for name, sub in [("overall", d), ("0.4-0.7", d[d["bin"] == "0.4-0.7"]), (">=0.7", d[d["bin"] == ">=0.7"])]:
    n = len(sub)
    for col in ["top1", "top5", "top10", "baseline_top10"]:
        k = int(sub[col].sum())
        lo, hi = wilson(k, n)
        lines.append(f"{name:8} {col:15} n={n:3}  {k / n:.3f}  (95% CI {lo:.2f}-{hi:.2f})")
text = "\n".join(lines)
print(text)
open("docs/validation_ci.txt", "w").write(text + "\n")