import csv
import itertools

# Define a polynomial ring with 6 variables (x0 to x5)
# to generate exactly 56 cubic monomials
R = PolynomialRing(QQ, 6, 'x')
x = R.gens()

print("Generating 56 cubic monomials...")
monomials = [prod(c) for c in itertools.combinations_with_replacement(x, 3)]
idx_to_mono = {i: m for i, m in enumerate(monomials)}

# Define the output file name
csv_filename = "cubic_monomials.csv"

# Write the indexed monomials to a 2-column CSV file
with open(csv_filename, mode='w', newline='') as file:
    writer = csv.writer(file)

    # Write the header row
    writer.writerow(["Index", "Monomial"])

    # Write the data rows
    for idx, mono in idx_to_mono.items():
        writer.writerow([idx, str(mono)])

print(f"Successfully saved {len(idx_to_mono)} monomials to '{csv_filename}'.")









