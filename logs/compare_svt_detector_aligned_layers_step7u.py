#!/usr/bin/env python3
"""Compare pre/post-Step-7u overlap identities, ignoring order and samples.

ROOT overlap identities are deterministic and must be identical. Geant4 uses
random surface samples, so its named warning set is reported for review rather
than required to match exactly.
"""

import collections
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BASELINE = ROOT / "logs/check_svt_opening_separation_step7n.txt"
CURRENT = ROOT / "logs/check_svt_detector_aligned_layers_step7u.txt"


def root_overlaps(text):
    return collections.Counter(
        re.findall(r"^ = Overlap ov\d+: (.*?) ovlp=[0-9.eE+-]+$", text, re.MULTILINE)
    )


def geant4_overlaps(text):
    return collections.Counter(
        re.findall(
            r"Overlap is detected for volume (.*?) \([^\n]+\) with its mother volume "
            r"(.*?) \([^\n]+\)",
            text,
        )
    )


baseline_text = BASELINE.read_text()
current_text = CURRENT.read_text()
baseline_root = root_overlaps(baseline_text)
current_root = root_overlaps(current_text)
baseline_g4 = geant4_overlaps(baseline_text)
current_g4 = geant4_overlaps(current_text)

print(f"ROOT baseline/current: {sum(baseline_root.values())}/{sum(current_root.values())}")
print(f"ROOT only-before/only-after: {sum((baseline_root-current_root).values())}/"
      f"{sum((current_root-baseline_root).values())}")
assert baseline_root == current_root

print(f"Geant4 baseline/current: {sum(baseline_g4.values())}/{sum(current_g4.values())}")
print(f"Geant4 only-before/only-after: {sum((baseline_g4-current_g4).values())}/"
      f"{sum((current_g4-baseline_g4).values())}")
print("PASS: deterministic ROOT overlap identities are unchanged")
print("NOTE: Geant4 pair differences are allowed because surface sampling is stochastic")
