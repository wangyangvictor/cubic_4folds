import itertools
import csv
from sage.all import *

# ==========================================
# 0. USER CONFIGURATION
# ==========================================
K_VAL = 3                # Change this to choose k monomials
START_INDEX = 0          # Change this to resume from a specific index (e.g., 500)
csv_filename = f"cubic_4fold_k{K_VAL}_data.csv"

# ==========================================
# 1. SETUP RINGS AND SPACES
# ==========================================
R = PolynomialRing(QQ, names=['x0','x1','x2','x3','x4','x5'], order='degrevlex')
x = R.gens()
P_X = ProjectiveSpace(QQ, 5, 'x')

BigRing = PolynomialRing(QQ, names=['x0','x1','x2','x3','x4','x5', 'y0','y1','y2','y3','y4','y5'], order='degrevlex')
X_vars = BigRing.gens()[:6]
Y_vars = BigRing.gens()[6:]

Y_ring = PolynomialRing(QQ, names=['y0','y1','y2','y3','y4','y5'], order='degrevlex')
P_Y = ProjectiveSpace(QQ, 5, 'y')

# ==========================================
# 2. LOAD PRE-COMPUTED ORBIT REPRESENTATIVES
# ==========================================
# We still need the monomial list to translate indices back to polynomials
print("Generating 56 cubic monomials...")
monomials = [prod(c) for c in itertools.combinations_with_replacement(x, 3)]
idx_to_mono = {i: m for i, m in enumerate(monomials)}

orbit_file = f"orbit_reps_k{K_VAL}.txt"
orbit_reps = []

print(f"Loading pre-computed orbits from {orbit_file}...")
with open(orbit_file, 'r') as f:
    for line in f:
        if line.strip():
            # Convert "0,1,2,3,4" into a tuple of integers (0, 1, 2, 3, 4)
            rep = tuple(map(int, line.strip().split(',')))
            orbit_reps.append(rep)

print(f"Loaded {len(orbit_reps)} distinct orbits for k={K_VAL}.")

# ==========================================
# 3. GEOMETRY CHECK FUNCTIONS
# ==========================================
def is_non_conical(partials):
    basis_deg2 = [prod(c) for c in itertools.combinations_with_replacement(x, 2)]
    M = matrix(QQ, 6, len(basis_deg2))
    for i, p in enumerate(partials):
        for j, b in enumerate(basis_deg2):
            M[i, j] = p.monomial_coefficient(b)
    return M.rank() == 6

def compute_invariants(F):
    results = {}

    # --- (b) & (c) Reduced and Irreducible ---
    factors = list(F.factor())
    results['reduced'] = all(mult == 1 for f, mult in factors)
    results['irreducible'] = (len(factors) == 1 and factors[0][1] == 1)

    # EARLY EXIT 1
    if not (results['reduced'] and results['irreducible']):
        return None

    # --- Partials ---
    partials = [derivative(F, v) for v in x]

    # --- (a) Non-conical ---
    results['non_conical'] = is_non_conical(partials)

    # EARLY EXIT 2
    if not results['non_conical']:
        return None

    # --- (1) & (2) Singular Subscheme Dimension & Degree ---
    J = R.ideal(partials)
    sing_scheme = P_X.subscheme(J)
    sing_dim = sing_scheme.dimension()
    results['sing_dim'] = sing_dim

    if sing_dim >= 0:
        results['sing_deg'] = sing_scheme.degree()
    else:
        return None

    # --- (3) & (4) Dual Variety Computation ---
    F_big = BigRing(F)
    partials_big = [derivative(F_big, xv) for xv in X_vars]

    minors = []
    for i in range(6):
        for j in range(i+1, 6):
            minors.append(Y_vars[i] * partials_big[j] - Y_vars[j] * partials_big[i])

    incidence = sum([X_vars[i] * Y_vars[i] for i in range(6)])
    I_graph = BigRing.ideal(minors + [F_big, incidence])

    J_big = BigRing.ideal(partials_big)

    I_sat = I_graph.saturation(J_big)[0]
    Dual_Ideal_Big = I_sat.elimination_ideal(X_vars)

    dual_gens_Y = [Y_ring(g) for g in Dual_Ideal_Big.gens()]
    Dual_Ideal = Y_ring.ideal(dual_gens_Y)

    dual_scheme = P_Y.subscheme(Dual_Ideal)
    dual_dim = dual_scheme.dimension()
    results['dual_dim'] = dual_dim

    if dual_dim >= 0:
        results['dual_deg'] = dual_scheme.degree()
    else:
        results['dual_deg'] = 0

    results['is_hypersurface'] = bool(dual_dim == 4)

    return results

# ==========================================
# 4. MAIN EXECUTION LOOP
# ==========================================
# If starting > 0, append to the file so we don't overwrite previous work.
file_mode = 'a' if START_INDEX > 0 else 'w'

print(f"\nProcessing forms starting from index {START_INDEX}... (Skipping conical, reducible, non-reduced)\n")

with open(csv_filename, mode=file_mode, newline='') as file:
    writer = csv.writer(file)

    # Only write the header if we are starting a fresh file
    if START_INDEX == 0:
        writer.writerow(['Index', 'Cubic Form F', 'Non-Conical', 'Reduced', 'Irreducible', 'Singular Dim', 'Singular Deg', 'Dual Deg', 'Dual Is Hypersurface'])

    try:
        for i, rep in enumerate(orbit_reps):
            # Skip iterations until we hit the desired starting index
            if i < START_INDEX:
                continue

            F = sum([idx_to_mono[idx] for idx in rep])

            inv = compute_invariants(F)

            if inv is None:
                continue

            writer.writerow([
                i,
                str(F),
                inv['non_conical'],
                inv['reduced'],
                inv['irreducible'],
                inv['sing_dim'],
                inv['sing_deg'],
                inv['dual_deg'],
                inv['is_hypersurface']
            ])

            file.flush()

    except KeyboardInterrupt:
        print(f"\nProcess manually interrupted. The data processed so far is safely saved in {csv_filename}.")

print(f"\nFinished! Your filtered data is available in {csv_filename}.")



