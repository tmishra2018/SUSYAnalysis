#input_file = "input.txt"
#output_file = "ULCrossSectionT5WG.txt"

input_file = "input_T6Wg.txt"
output_file = "ULCrossSectionT6WG.txt"

def process_file(input_file, output_file):
    try:
        with open(input_file, "r", encoding="utf-8") as infile, open(output_file, "w", encoding="utf-8") as outfile:
            for line in infile:
                line = line.strip()
                if not line:
                    continue
                parts = line.split()
                if len(parts) == 3:
                    num1 = int(parts[0])
                    num2 = float(parts[1])*1000
                    percent = float(parts[2]) / 100  
                    result = num2 * percent
                    formatted_line = f"{num1} {num2:.3E} {result:.4f}\n"
                    outfile.write(formatted_line)
    except Exception as e:
        print(f"Error: {e}")
process_file(input_file, output_file)
