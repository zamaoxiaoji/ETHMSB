#!/usr/bin/env python3
import argparse
import csv
import os
import subprocess
import sys
import time
from pathlib import Path
from statistics import median


ROOT = Path(__file__).resolve().parents[1]
WAIT_TIMEOUT_SECONDS = 0


def run(cmd, cwd=ROOT, stdout_path=None):
    print("+ " + " ".join(str(x) for x in cmd), flush=True)
    if stdout_path is None:
        return subprocess.run(cmd, cwd=cwd, check=True)
    with stdout_path.open("w") as out:
        proc = subprocess.run(
            cmd,
            cwd=cwd,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
        out.write(proc.stdout)
        sys.stdout.write(proc.stdout)
        if proc.returncode != 0:
            raise subprocess.CalledProcessError(proc.returncode, cmd)
        return proc


def run_and_tee(cmd, output_path, cwd=ROOT):
    print("+ " + " ".join(str(x) for x in cmd) + f" | tee {output_path}", flush=True)
    proc = subprocess.run(
        cmd,
        cwd=cwd,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    output_path.write_text(proc.stdout)
    sys.stdout.write(proc.stdout)
    if proc.returncode != 0:
        raise subprocess.CalledProcessError(proc.returncode, cmd)
    return proc


def loadavg_1m():
    with Path("/proc/loadavg").open() as f:
        return float(f.read().split()[0])


def wait_for_idle(max_load, poll_seconds):
    if max_load <= 0:
        return
    start = time.monotonic()
    while True:
        load = loadavg_1m()
        print(f"load1={load:.2f}, target<={max_load:.2f}", flush=True)
        if load <= max_load:
            return
        if WAIT_TIMEOUT_SECONDS and time.monotonic() - start >= WAIT_TIMEOUT_SECONDS:
            raise TimeoutError(
                f"load1 stayed above {max_load:.2f} for {WAIT_TIMEOUT_SECONDS}s"
            )
        time.sleep(poll_seconds)


def capture_text(cmd, cwd=ROOT):
    proc = subprocess.run(
        cmd,
        cwd=cwd,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    return proc.stdout


def write_environment_snapshot(results_dir, tag):
    snapshot = results_dir / f"environment_{tag}.txt"
    with snapshot.open("w") as out:
        out.write("## uptime\n")
        out.write(capture_text(["uptime"]))
        out.write("\n## cpu count\n")
        out.write(str(os.cpu_count() or 1) + "\n")
        out.write("\n## top cpu processes\n")
        out.write(
            capture_text([
                "ps", "-eo", "user,pid,ppid,psr,pcpu,pmem,comm,args",
                "--sort=-pcpu",
            ]).splitlines()[0]
            + "\n"
        )
        top = capture_text([
            "ps", "-eo", "user,pid,ppid,psr,pcpu,pmem,comm,args",
            "--sort=-pcpu",
        ]).splitlines()[1:31]
        out.write("\n".join(top) + "\n")
        out.write("\n## optimized cmake flags\n")
        out.write(
            capture_text([
                "cmake", "-LAH", "-N", "build-he3db-optimized",
            ])
        )
        out.write("\n## release cmake flags\n")
        out.write(
            capture_text([
                "cmake", "-LAH", "-N", "build-he3db-release",
            ])
        )
        out.write("\n## baseline git status\n")
        out.write(capture_text(["git", "status", "--short"]))
    print(f"wrote {snapshot}", flush=True)
    return snapshot


def parse_bench_result(log_path):
    line = ""
    for raw in log_path.read_text().splitlines():
        if raw.startswith("BENCH_RESULT,"):
            line = raw
    if not line:
        raise RuntimeError(f"missing BENCH_RESULT in {log_path}")
    fields = line.split(",")
    return {
        "query": fields[1],
        "rows": int(fields[2]),
        "threads": int(fields[3]),
        "filter_ms": fields[4],
        "aggregation_ms": fields[5],
        "total_ms": fields[6],
        "metric1": fields[7],
        "metric2": fields[8],
        "plain": fields[9],
        "encrypted": fields[10],
        "final_abs_error": fields[11],
        "success": fields[12],
    }


def read_micro_csv(path):
    rows = {}
    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            bit = int(row["bit"])
            rows[bit] = (
                float(row["he3db_ms"]),
                float(row["ethmsb_ms"]),
            )
    return rows


def summarize_values(values, method):
    if method == "min":
        return min(values)
    if method == "median":
        return median(values)
    raise ValueError(f"unknown summary method: {method}")


def write_micro_summary(raw_paths, out_path, method):
    runs = [read_micro_csv(path) for path in raw_paths]
    bits = sorted(set.intersection(*(set(run.keys()) for run in runs)))
    with out_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["bit", "he3db_ms", "ethmsb_ms", "speedup"])
        for bit in bits:
            he = summarize_values([run[bit][0] for run in runs], method)
            eth = summarize_values([run[bit][1] for run in runs], method)
            writer.writerow([bit, he, eth, he / eth])
    print(f"wrote {out_path} from {len(raw_paths)} run(s) using {method}", flush=True)


def run_micro(build_dir, trials, results_dir, tag, repeats, summary_method):
    run(["cmake", "--build", build_dir, "--target", "msb_bit_timing_single_core", "-j", str(os.cpu_count() or 1)])
    bin_path = Path(build_dir) / "bin" / "msb_bit_timing_single_core"
    outputs = {
        "gt": "greater_than_input_bit_timing_optimized_1_32.csv",
        "ge": "greater_equal_input_bit_timing_optimized_1_32.csv",
        "eq": "equal_input_bit_timing_optimized_1_32.csv",
        "ne": "not_equal_input_bit_timing_optimized_1_32.csv",
    }
    for op, name in outputs.items():
        raw_paths = []
        for repeat in range(1, repeats + 1):
            out_path = results_dir / name
            if repeats > 1:
                out_path = results_dir / f"{Path(name).stem}_{tag}_run{repeat}.csv"
            run([
                str(bin_path),
                "--op", op,
                "--from-bit", "1",
                "--to-bit", "32",
                "--trials", str(trials),
                "--out", str(out_path),
            ])
            raw_paths.append(out_path)
        if repeats > 1:
            write_micro_summary(raw_paths, results_dir / name, summary_method)


def run_tpch(build_dir, results_dir, tag, threads, agg_threads, pair_threads):
    targets = ["tpch_q6_he3db_omp", "tpch_q6_ethmsb", "tpch_q14_he3db_omp", "tpch_q14_ethmsb"]
    run(["cmake", "--build", build_dir, "--target", *targets, "-j", str(os.cpu_count() or 1)])
    rows_list = [1024, 2048, 4096, 8192, 16384]
    bins = {
        ("Q6", "HE3DB"): Path(build_dir) / "bin" / "tpch_q6_he3db_omp",
        ("Q6", "ETHMSB"): Path(build_dir) / "bin" / "tpch_q6_ethmsb",
        ("Q14", "HE3DB"): Path(build_dir) / "bin" / "tpch_q14_he3db_omp",
        ("Q14", "ETHMSB"): Path(build_dir) / "bin" / "tpch_q14_ethmsb",
    }
    csv_path = results_dir / f"tpch_q6_q14_fullcore_repro_{tag}.csv"
    with csv_path.open("w", newline="") as f:
        fieldnames = [
            "query", "impl", "rows", "threads", "agg_threads", "pair_threads",
            "filter_ms", "aggregation_ms", "total_ms", "metric1", "metric2", "plain", "encrypted",
            "final_abs_error", "success", "log",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for query in ("Q6", "Q14"):
            for rows in rows_list:
                for impl in ("HE3DB", "ETHMSB"):
                    log_path = results_dir / f"tpch_{query}_{impl}_{rows}_{tag}.log"
                    agg = agg_threads[(query, impl)]
                    cmd = [
                        str(bins[(query, impl)]),
                        "--rows", str(rows),
                        "--threads", str(threads),
                        "--agg-threads", str(agg),
                    ]
                    pair = ""
                    if (query, impl) in pair_threads:
                        pair = pair_threads[(query, impl)]
                        cmd.extend(["--pair-threads", str(pair)])
                    run(cmd, stdout_path=log_path)
                    parsed = parse_bench_result(log_path)
                    parsed["impl"] = impl
                    parsed["agg_threads"] = agg
                    parsed["pair_threads"] = pair
                    parsed["log"] = str(log_path)
                    writer.writerow(parsed)
                    f.flush()
    return csv_path


def main():
    global WAIT_TIMEOUT_SECONDS
    parser = argparse.ArgumentParser()
    parser.add_argument("--skip-wait", action="store_true")
    parser.add_argument("--max-load", type=float, default=8.0)
    parser.add_argument("--poll-seconds", type=int, default=60)
    parser.add_argument("--max-wait-seconds", type=int, default=0)
    parser.add_argument("--micro-trials", type=int, default=10)
    parser.add_argument("--micro-repeats", type=int, default=1)
    parser.add_argument("--micro-summary", choices=["median", "min"], default="median")
    parser.add_argument("--threads", type=int, default=os.cpu_count() or 1)
    parser.add_argument("--q6-he3db-agg-threads", type=int)
    parser.add_argument("--q6-ethmsb-agg-threads", type=int)
    parser.add_argument("--q14-he3db-agg-threads", type=int)
    parser.add_argument("--q14-ethmsb-agg-threads", type=int)
    parser.add_argument("--q14-ethmsb-pair-threads", type=int)
    parser.add_argument("--tag", default=time.strftime("%Y%m%d_%H%M%S"))
    parser.add_argument("--micro-only", action="store_true")
    parser.add_argument("--tpch-only", action="store_true")
    args = parser.parse_args()
    WAIT_TIMEOUT_SECONDS = args.max_wait_seconds

    results_dir = ROOT / "results"
    results_dir.mkdir(exist_ok=True)
    if not args.skip_wait:
        wait_for_idle(args.max_load, args.poll_seconds)
    write_environment_snapshot(results_dir, args.tag)

    if not args.tpch_only:
        if args.micro_repeats < 1:
            raise ValueError("--micro-repeats must be positive")
        run_micro(
            "build-he3db-optimized",
            args.micro_trials,
            results_dir,
            args.tag,
            args.micro_repeats,
            args.micro_summary,
        )
    tpch_csv = None
    if not args.micro_only:
        agg_threads = {
            ("Q6", "HE3DB"): args.q6_he3db_agg_threads or args.threads,
            ("Q6", "ETHMSB"): args.q6_ethmsb_agg_threads or args.threads,
            ("Q14", "HE3DB"): args.q14_he3db_agg_threads or args.threads,
            ("Q14", "ETHMSB"): args.q14_ethmsb_agg_threads or args.threads,
        }
        pair_threads = {}
        if args.q14_ethmsb_pair_threads is not None:
            pair_threads[("Q14", "ETHMSB")] = args.q14_ethmsb_pair_threads
        tpch_csv = run_tpch(
            "build-he3db-release",
            results_dir,
            args.tag,
            args.threads,
            agg_threads,
            pair_threads,
        )

    compare_cmd = ["python3", "scripts/compare_original_latency.py"]
    if tpch_csv is not None:
        compare_cmd.extend(["--tpch-csv", str(tpch_csv)])
    run_and_tee(compare_cmd, results_dir / f"compare_original_latency_{args.tag}.txt")


if __name__ == "__main__":
    main()
