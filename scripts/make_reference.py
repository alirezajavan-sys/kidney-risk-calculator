"""Generate test/reference_values.json from the tested Python implementation (ckdtools).

    pip install git+https://github.com/alirezajavan-sys/ckd-toolkit.git
    python scripts/make_reference.py

The JavaScript calculator must reproduce every value (see test/parity.test.mjs).
"""
import itertools
import json
import warnings
from pathlib import Path

from ckdtools import albuminuria_category, egfr_creatinine, egfr_creatinine_cystatin, egfr_cystatin, gfr_category, kdigo_risk, kfre

warnings.simplefilter("ignore")
cases = []
for age, sex, scr, scys in itertools.product([18, 35, 52, 67, 84], ["female", "male"], [0.5, 0.7, 0.9, 1.2, 2.1, 4.5],
                                              [0.6, 0.8, 1.4, 2.5]):
    cases.append({
        "kind": "egfr", "age": age, "sex": sex, "scr": scr, "scys": scys,
        "egfr_cr": egfr_creatinine(scr, age, sex),
        "egfr_cys": egfr_cystatin(scys, age, sex),
        "egfr_cr_cys": egfr_creatinine_cystatin(scr, scys, age, sex),
    })
for age, sex, egfr, acr, region in itertools.product([45, 60, 75], ["female", "male"], [12, 25, 40, 55, 70],
                                                     [5, 30, 120, 450, 2000], ["non_north_american", "north_american"]):
    cases.append({
        "kind": "kfre", "age": age, "sex": sex, "egfr": egfr, "acr": acr, "region": region,
        "kfre2": kfre(age, sex, egfr, acr, years=2, region=region),
        "kfre5": kfre(age, sex, egfr, acr, years=5, region=region),
    })
for egfr, acr in itertools.product([95, 90, 89.9, 60, 59.9, 45, 44.9, 30, 29.9, 15, 14.9, 5], [10, 29.9, 30, 300, 301]):
    g, a = gfr_category(egfr), albuminuria_category(acr)
    cases.append({"kind": "kdigo", "egfr": egfr, "acr": acr, "g": g, "a": a, "risk": kdigo_risk(g, a)})

out = Path(__file__).resolve().parents[1] / "test" / "reference_values.json"
out.write_text(json.dumps(cases, indent=1))
print(f"{len(cases)} reference cases written to {out}")
