#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path


ORIG_REL_HE3DB = {
    1: 11, 2: 11, 3: 11, 4: 10, 5: 31, 6: 31, 7: 31, 8: 31,
    9: 51, 10: 114, 11: 112, 12: 113, 13: 115, 14: 188, 15: 189,
    16: 193, 17: 194, 18: 194, 19: 266, 20: 266, 21: 266,
    22: 266, 23: 271, 24: 334, 25: 339, 26: 339, 27: 340,
    28: 340, 29: 411, 30: 412, 31: 411, 32: 411,
}

ORIG_REL_ETHMSB = {
    1: 9.6, 2: 9.58, 3: 9.76, 4: 9.59, 5: 19.55, 6: 20.00,
    7: 19.31, 8: 19.44, 9: 51.02, 10: 63.26, 11: 63.30,
    12: 63.61, 13: 63.61, 14: 100.27, 15: 100.22, 16: 100.36,
    17: 100.34, 18: 100.18, 19: 137.51, 20: 137.29, 21: 137.50,
    22: 136.95, 23: 136.66, 24: 175.57, 25: 175.94, 26: 175.33,
    27: 175.54, 28: 175.84, 29: 215.91, 30: 215.41, 31: 215.60,
    32: 215.52,
}

ORIG_EQ_HE3DB = {
    1: 30, 2: 31, 3: 30, 4: 30, 5: 71, 6: 74, 7: 71, 8: 81,
    9: 114, 10: 242, 11: 233, 12: 237, 13: 233, 14: 383, 15: 383,
    16: 382, 17: 381, 18: 380, 19: 526, 20: 526, 21: 526,
    22: 527, 23: 529, 24: 676, 25: 669, 26: 671, 27: 676,
    28: 671, 29: 813, 30: 817, 31: 813, 32: 819,
}

ORIG_EQ_ETHMSB = {
    1: 30.12, 2: 30.32, 3: 30.20, 4: 30.07, 5: 50.69, 6: 50.52,
    7: 50.29, 8: 50.72, 9: 112.71, 10: 138.28, 11: 138.36,
    12: 138.58, 13: 138.83, 14: 214.84, 15: 214.88, 16: 214.94,
    17: 214.75, 18: 214.96, 19: 404.54, 20: 404.26, 21: 403.93,
    22: 404.88, 23: 404.50, 24: 521.61, 25: 522.09, 26: 521.27,
    27: 521.67, 28: 520.78, 29: 640.35, 30: 640.70, 31: 640.25,
    32: 640.59,
}

ORIG_TPCH = {
    ("Q6", 1024): (65.531, 61.65),
    ("Q6", 2048): (72.596, 61.04),
    ("Q6", 4096): (89.351, 71.815),
    ("Q6", 8192): (128.37, 94.553),
    ("Q6", 16384): (202.435, 132.45),
    ("Q14", 1024): (130.027, 126.646),
    ("Q14", 2048): (151.46, 136.052),
    ("Q14", 4096): (197.038, 161.808),
    ("Q14", 8192): (301.265, 220.094),
    ("Q14", 16384): (574.243, 338.078),
}


def read_bit_csv(path):
    rows = {}
    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            bit = int(row["bit"])
            rows[bit] = (
                float(row["he3db_ms"]),
                float(row["ethmsb_ms"]),
                float(row["speedup"]),
            )
    return rows


def pct(value, target):
    return (value - target) / target * 100.0


def print_bit_report(label, rows, orig_he, orig_eth, key_bits):
    bits = sorted(set(rows) & set(orig_he) & set(orig_eth))
    he_errors = [abs(pct(rows[b][0], orig_he[b])) for b in bits]
    eth_errors = [abs(pct(rows[b][1], orig_eth[b])) for b in bits]
    print(f"[{label}]")
    print(
        f"  HE3DB avg_abs_err={sum(he_errors) / len(he_errors):.2f}% "
        f"max_abs_err={max(he_errors):.2f}%"
    )
    print(
        f"  ETHMSB avg_abs_err={sum(eth_errors) / len(eth_errors):.2f}% "
        f"max_abs_err={max(eth_errors):.2f}%"
    )
    for bit in key_bits:
        if bit not in rows:
            continue
        he, eth, speedup = rows[bit]
        print(
            f"  bit={bit:2d} "
            f"HE3DB={he:8.3f}/{orig_he[bit]:8.3f} ({pct(he, orig_he[bit]):+6.1f}%) "
            f"ETHMSB={eth:8.3f}/{orig_eth[bit]:8.3f} ({pct(eth, orig_eth[bit]):+6.1f}%) "
            f"speedup={speedup:.3f}"
        )
    print()


def read_tpch_csv(path):
    rows = {}
    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            query = row["query"]
            impl = row["impl"]
            row_count = int(row["rows"])
            rows[(query, impl, row_count)] = float(row["total_ms"]) / 1000.0
    return rows


def print_tpch_report(path):
    if not path.exists():
        return
    rows = read_tpch_csv(path)
    print("[TPC-H]")
    for query in ("Q6", "Q14"):
        for row_count in (1024, 2048, 4096, 8192, 16384):
            target = ORIG_TPCH[(query, row_count)]
            he = rows.get((query, "HE3DB", row_count))
            eth = rows.get((query, "ETHMSB", row_count))
            if he is None or eth is None:
                continue
            print(
                f"  {query} rows={row_count:5d} "
                f"HE3DB={he:8.3f}/{target[0]:8.3f} ({pct(he, target[0]):+6.1f}%) "
                f"ETHMSB={eth:8.3f}/{target[1]:8.3f} ({pct(eth, target[1]):+6.1f}%)"
            )
    print()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-dir", default="results")
    parser.add_argument("--tpch-csv", default="results/tpch_q6_q14_fullcore_release_20260626_102801.csv")
    args = parser.parse_args()

    results = Path(args.results_dir)
    key_bits = [4, 8, 9, 16, 24, 32]
    bit_inputs = [
        ("rel/gt", "greater_than_input_bit_timing_optimized_1_32.csv", ORIG_REL_HE3DB, ORIG_REL_ETHMSB),
        ("rel/ge", "greater_equal_input_bit_timing_optimized_1_32.csv", ORIG_REL_HE3DB, ORIG_REL_ETHMSB),
        ("eq", "equal_input_bit_timing_optimized_1_32.csv", ORIG_EQ_HE3DB, ORIG_EQ_ETHMSB),
        ("ne", "not_equal_input_bit_timing_optimized_1_32.csv", ORIG_EQ_HE3DB, ORIG_EQ_ETHMSB),
    ]
    for label, name, orig_he, orig_eth in bit_inputs:
        path = results / name
        if path.exists():
            print_bit_report(label, read_bit_csv(path), orig_he, orig_eth, key_bits)
    print_tpch_report(Path(args.tpch_csv))


if __name__ == "__main__":
    main()
