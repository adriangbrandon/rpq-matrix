import sys
import os
import math

TIME_OUT = 60000000000
NQ = 21

if len(sys.argv) < 4:
    print(f"Usage: {sys.argv[0]} <query-results-file> <query-types-file> <selected-queries-file>")
    sys.exit(1)

line = ""
n_line = 0
q = 0

M = {
    "v /* c": 1,
    "v * c": 2,
    "v + c": 3,
    "c * v": 4,
    "c /* v": 5,
    "v / c": 6,
    "v / v": 7,
    "v */* c": 8,
    "v |* c": 9,
    "v | v": 10,
    "v ^ v": 11,
    "v /* v": 12,
    "v */*/*/*/* c": 13,
    "v * v": 14,
    "v /? c": 15,
    "v + v": 16,
    "v /+ c": 17,
    "v ^/ c": 18,
    "v || v": 19,
    "v /^ v": 20,
    "v | c": 21
}

fp = [None] + [open(sys.argv[1] + "-" + str(i), "w") for i in range(1, NQ + 3)]
for i in range(1, NQ + 3):
    fp[i].write(f"\\begin{{filecontents}}{{q{i}-{sys.argv[1]}.dat}}\n")

fp_all = open(sys.argv[1] + "-all", "w")
fp_v_to_v = open(sys.argv[1] + "-all_v_to_v", "w")
fp_c_to_v = open(sys.argv[1] + "-all_c_to_v", "w")

m_excluded = set()
with open(sys.argv[3], "r") as ifs_excluded:
    for line in ifs_excluded:
        eq = int(line.strip())
        m_excluded.add(eq)
        print(f"Including query {eq}...")

nq_to_qtype = {}
with open(sys.argv[2], "r") as ifs_types:
    for i in range(1, 2111):
        line = ifs_types.readline().strip()
        nq_to_qtype[i] = line

with open(sys.argv[1], "r") as ifs_results:
    for i in range(1, 2111):
        nQ, semicolon, nResults, semicolon, qTime = map(int, ifs_results.readline().strip().split())

        if i not in m_excluded:
            continue

        if nResults == -1 or nResults == -2 or qTime > TIME_OUT:
            qTime = TIME_OUT

        fp_all.write(f"{qTime / 1000000000:.9f}\n")

        if nq_to_qtype[nQ][0] == "v" and nq_to_qtype[nQ][-1] == "v":
            fp_v_to_v.write(f"{qTime / 1000000000:.9f}\n")
        else:
            fp_c_to_v.write(f"{qTime / 1000000000:.9f}\n")

        if nq_to_qtype[nQ] in M:
            fp[M[nq_to_qtype[nQ]]].write(f"{qTime / 1000000000:.9f}\n")

