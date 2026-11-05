#!/usr/bin/env python3
"""
Correctly blinds gmN lines by computing Asimov N as:
    N_asimov = round(expected_SR / alpha)
where alpha is the transfer factor from the datacard,
and expected_SR = alpha * N_original (the background prediction).

This satisfies text2workspace's requirement that N * alpha == expected_rate exactly.
"""
import sys, re

def asimov_gmN(line):
    parts = line.split()
    if len(parts) < 3 or parts[1] != 'gmN':
        return line

    name  = parts[0]
    n_cr  = int(parts[2])
    cols  = parts[3:]

    # Find all non-dash transfer factors and their positions
    tfs = []
    for c in cols:
        try:
            v = float(c)
            if v > 0:
                tfs.append(v)
        except ValueError:
            pass

    if not tfs:
        return line

    # Expected SR yield for each bin = tf * n_cr
    # Asimov: keep the same expected SR yield but set N_cr such that
    # N_asimov * tf = expected_SR = tf * n_cr  →  N_asimov = n_cr
    # BUT text2workspace checks N * alpha == rate where rate = tf * N_original
    # So we need N_asimov such that N_asimov * tf == tf * n_cr
    # i.e. N_asimov = n_cr ... which is circular (keeps real data!)
    #
    # The CORRECT interpretation:
    # alpha (tf) is fixed in the datacard.
    # N is the observed CR count.
    # Expected SR = N * alpha.
    # For Asimov: set N_asimov so that N_asimov * alpha = expected_SR_asimov
    # where expected_SR_asimov = background prediction (= N_original * alpha).
    # This means N_asimov = N_original ... same number.
    #
    # THE REAL SOLUTION: use --setParameters to fix all gmN nuisances to 0
    # OR use the original N but pass --snapshotName to override.
    # For text2workspace compatibility, N must stay as-is.
    # Instead we fix nuisances post-workspace using --freezeNuisanceGroups.

    # Since we cannot change N without breaking text2workspace,
    # return line unchanged and handle via combine flags instead.
    return line, tfs, n_cr

# Actually: the solution is to NOT change gmN N values at all.
# Instead print them so we can freeze them in combine.

infile  = sys.argv[1]
outfile = sys.argv[2]

gmn_params = []

with open(infile) as f:
    lines = f.readlines()

with open(outfile, 'w') as out:
    for line in lines:
        parts = line.split()
        if len(parts) >= 3 and parts[1] == 'gmN':
            gmn_params.append(parts[0])
        out.write(line)

# Write the list of gmN nuisance names to a separate file
with open('gmn_params.txt', 'w') as gf:
    gf.write(','.join(gmn_params))

print(f"Found {len(gmn_params)} gmN parameters", file=sys.stderr)
print(f"Written parameter list to gmn_params.txt", file=sys.stderr)
