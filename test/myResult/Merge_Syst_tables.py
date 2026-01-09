import os

# Directory containing input tables
directory = "logs"

# Input files and labels (ordered)
files = [
    ("eventcount_2016preVFP.txt", "2016 preVFP"),
    ("eventcount_2016postVFP.txt", "2016 postVFP"),
    ("eventcount_2017.txt", "2017"),
    ("eventcount_2018.txt", "2018"),
]

# --------------------------------------
# Extract rows inside tabular environment
# --------------------------------------
def extract_table_rows(text):
    start = text.find(r"\begin{tabular")
    end = text.find(r"\end{tabular}")

    if start == -1 or end == -1:
        return []

    block = text[start:end]
    lines = block.splitlines()

    cleaned = []
    for ln in lines:
        ln = ln.strip()
        if not ln:
            continue
        if ln.startswith(r"\begin{tabular"):
            continue
        if ln.startswith(r"\end{tabular"):
            continue
        if "Source of Uncertainties" in ln:
            continue
        if ln.startswith(r"\hline"):
            continue
        cleaned.append(ln)

    return cleaned


# --------------------------------------
# Read tables
# --------------------------------------
year_tables = []

for fname, label in files:
    path = os.path.join(directory, fname)
    with open(path, "r") as f:
        rows = extract_table_rows(f.read())

    if not rows:
        raise RuntimeError(f"❌ No table extracted from {fname}")

    year_tables.append((label, rows))


# --------------------------------------
# Build merged LaTeX table
# --------------------------------------
merged = r"""
\begin{table}[h]
\centering
\resizebox{\linewidth}{!}{
\begin{tabular}{|c||c|c|c|c|c|}
\hline
\multicolumn{6}{|c|}{\textbf{Systematic uncertainties of the SM backgrounds}} \\
\hline
\textbf{Source of Uncertainties} &
$e\rightarrow\gamma$ fakes &
$\text{jet}\rightarrow\gamma$ fakes &
$\text{jet}\rightarrow l$ fakes &
$V+\gamma$ &
rare EWK \\
\hline
"""

for label, rows in year_tables:
    merged += rf"\multicolumn{{6}}{{|c|}}{{\textbf{{{label}}}}} \\"
    merged += "\n\\hline\n"
    for r in rows:
        merged += f"  {r}\n"
    merged += "\\hline\n"

merged += r"""
\end{tabular}
}
\caption{Systematic uncertainties of the SM backgrounds for different years.}
\label{table:ch4-systematic}
\end{table}
"""

# --------------------------------------
# Write output
# --------------------------------------
with open("merged_table.tex", "w") as f:
    f.write(merged)

print("✅ merged_table.tex created successfully")
