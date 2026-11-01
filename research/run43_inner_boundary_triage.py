"""Read-only triage of the original Run43 inner-band acceptance deficits.

This is a diagnostic, NOT a replacement for the frozen audit or a new
physical solve. Input values are the previously reported October 22 audit
orders; update them only by supplying verified original physical audit output.
"""
from __future__ import annotations

THRESHOLD = 1.5
PAIRS = ("8->16", "16->32", "32->64")
INNER = {
    "H": (1.4888, 1.7354, 1.8671),
    "M2": (0.5082, 0.8524, 1.0922),
}

def main() -> None:
    print("Original physical inner-boundary audit: triage only")
    failures = 0
    for metric, orders in INNER.items():
        for pair, order in zip(PAIRS, orders):
            deficit = max(0.0, THRESHOLD - order)
            failures += deficit > 0
            print(f"{metric:>2} {pair:>6}: p={order:.4f}  deficit={deficit:.4f}  "
                  f"{'FAIL' if deficit else 'PASS'}")
    print(f"Inner-band failed gates: {failures}/6")
    print("Next: inspect boundary stencil, tensor-divergence consistency, and "
          "source-coupled residuals on identical frozen physical bands.")
    print("No equations, thresholds, sources, checkpoints, or audit files changed.")

if __name__ == "__main__":
    main()
