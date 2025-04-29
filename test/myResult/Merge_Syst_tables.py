import os

# Define the directory containing the text files
directory = 'logs/'

# Define the filenames
filenames = [
    'eventcount_2016postVFP.txt',
    'eventcount_2016preVFP.txt',
    'eventcount_2017.txt',
    'eventcount_2018.txt'
]

# Initialize an empty list to store the tables
tables = []

# Read each file and extract the table content
for filename in filenames:
    with open(os.path.join(directory, filename), 'r') as file:
        content = file.read()
        start = content.find('\\begin{tabular}')
        end = content.find('\\end{tabular}') + len('\\end{tabular}')
        table_content = content[start:end]
        tables.append(table_content)

# Merge the tables into a single table
merged_table = """
\\begin{table}[h]
\\centering
\\resizebox{\\linewidth}{!}{
\\begin{tabular}{|c||c|c|c|c|c|}
\\hline
\\multicolumn{6}{|c|}{\\textbf{Systematic uncertainties of the SM backgrounds (2016 preVFP)}} \\\\
\\hline
\\hline
  Source of Uncertainties & $e\\rightarrow\\gamma$ fakes & $\\text{jet}\\rightarrow\\gamma$ fakes & $\\text{jet}\\rightarrow l$ fakes & $V+\\gamma$ & rare EWK \\\\
  \\hline
""" + tables[1][len('\\begin{tabular}{|c|c|c|c|c|c|}'):] + """
\\multicolumn{6}{|c|}{\\textbf{Systematic uncertainties of the SM backgrounds (2016 postVFP)}} \\\\
\\hline
""" + tables[0][len('\\begin{tabular}{|c|c|c|c|c|c|}'):] + """
\\multicolumn{6}{|c|}{\\textbf{Systematic uncertainties of the SM backgrounds (2017)}} \\\\
\\hline
""" + tables[2][len('\\begin{tabular}{|c|c|c|c|c|c|}'):] + """
\\multicolumn{6}{|c|}{\\textbf{Systematic uncertainties of the SM backgrounds (2018)}} \\\\
\\hline
""" + tables[3][len('\\begin{tabular}{|c|c|c|c|c|c|}'):] + """
\\end{tabular}
}
\\caption{Systematic uncertainties of the SM backgrounds for different years.}
\\label{table:ch4-systematic}
\\end{table}
"""

# Save the merged table to a new file
with open('merged_table.tex', 'w') as file:
    file.write(merged_table)

print("The tables have been successfully merged into merged_table.tex.")
