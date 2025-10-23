#!/usr/bin/env python3
"""
compute_convergence_dmstag.py

命令行工具：对给定的可执行命令模板（包含占位符 {N} 表示网格大小），按多个网格尺寸运行可执行程序，解析输出中的 L2 误差并计算收敛阶。

用法示例：
  python3 tools/compute_convergence_dmstag.py "./build/src/poisson/poisson_dmstag_tests_cli --grid {N} --pair poly2 --poisson_check_error" 16 32 64 128

要求：可执行程序需要在输出中包含一行类似：
  ||u - u_exact||_2 = 2.51004e-05
脚本默认通过该模式提取 L2 误差。如果输出格式不同，可以使用 --pattern 自定义正则表达式（第一个捕获组为误差值）。
"""

import argparse
import shlex
import subprocess
import sys
import re
import math
from statistics import mean

DEFAULT_PATTERN = r"\|\|u - u_exact\|\|_2\s*=\s*([0-9Ee+\-\.]+)"


def run_and_extract(cmd, pattern, timeout):
    try:
        proc = subprocess.run(shlex.split(cmd), capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        print(f"ERROR: command timed out: {cmd}", file=sys.stderr)
        return None, None
    out = proc.stdout + "\n" + proc.stderr
    m = re.search(pattern, out)
    if m:
        try:
            val = float(m.group(1))
            return val, out
        except ValueError:
            print(f"ERROR: cannot parse numeric value from '{m.group(1)}'", file=sys.stderr)
            return None, out
    else:
        return None, out


def compute_rates(ns, errs):
    # compute local rates between consecutive pairs and global least-squares fit
    rates = []
    for i in range(1, len(ns)):
        e0 = errs[i-1]
        e1 = errs[i]
        n0 = ns[i-1]
        n1 = ns[i]
        if e0 is None or e1 is None or e0 <= 0 or e1 <= 0:
            rates.append(None)
            continue
        p = math.log(e0 / e1) / math.log(n1 / n0)
        rates.append(p)
    # global fit: linear regression of log(error) vs log(h) where h=1/N
    xs = []
    ys = []
    for n, e in zip(ns, errs):
        if e is None or e <= 0:
            continue
        xs.append(math.log(1.0 / n))
        ys.append(math.log(e))
    global_p = None
    if len(xs) >= 2:
        xm = mean(xs)
        ym = mean(ys)
        num = sum((x - xm) * (y - ym) for x, y in zip(xs, ys))
        den = sum((x - xm) ** 2 for x in xs)
        if den != 0:
            slope = num / den
            global_p = slope  # slope is p because y = p * log(h) + C
    return rates, global_p


def main():
    parser = argparse.ArgumentParser(description='Run executable at multiple grids and compute convergence order')
    parser.add_argument('cmd_template', help='Command template containing {N} placeholder for grid size (e.g. "./exe --grid {N} --poisson_check_error")')
    parser.add_argument('Ns', nargs='+', type=int, help='Grid sizes (e.g. 16 32 64 128)')
    parser.add_argument('--pattern', default=DEFAULT_PATTERN, help='Regex to extract L2 error; first capture group must be the numeric error')
    parser.add_argument('--timeout', type=int, default=300, help='Timeout in seconds for each run')
    parser.add_argument('--show-output', action='store_true', help='Show program output for each run')
    args = parser.parse_args()

    ns = args.Ns
    errs = []
    outputs = []
    for n in ns:
        cmd = args.cmd_template.replace('{N}', str(n))
        print(f"Running N={n}: {cmd}")
        val, out = run_and_extract(cmd, args.pattern, args.timeout)
        if val is None:
            print(f"WARNING: failed to extract error for N={n}")
        else:
            print(f"  extracted L2 error = {val:e}")
        errs.append(val)
        outputs.append(out)
        if args.show_output and out is not None:
            print('--- program output start ---')
            print(out)
            print('--- program output end ---')

    rates, global_p = compute_rates(ns, errs)

    print('\nConvergence summary:')
    print(f"{'N':>8} {'L2 error':>16} {'rate (vs prev)':>18}")
    for i, n in enumerate(ns):
        e = errs[i]
        r = rates[i-1] if i >= 1 else None
        e_str = f"{e:.5e}" if e is not None else 'N/A'
        r_str = f"{r:.4f}" if r is not None else 'N/A'
        print(f"{n:8d} {e_str:16s} {r_str:18s}")

    if global_p is not None:
        print(f"\nEstimated global convergence order p (least squares fit): {global_p:.4f}")
    else:
        print('\nNot enough valid data to estimate global convergence order.')

    # exit code: 0 if at least two valid errors found
    valid_count = sum(1 for e in errs if e is not None and e > 0)
    if valid_count >= 2:
        sys.exit(0)
    else:
        sys.exit(2)


if __name__ == '__main__':
    main()
