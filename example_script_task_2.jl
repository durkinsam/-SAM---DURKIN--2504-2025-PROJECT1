# example_script_task_2.jl
# Task 2: Parameterising coefficient type and showing why we need BigInt.
# We construct 5 explicit polynomial-addition examples whose coefficients overflow
# if stored as Int, then re-run the same additions with BigInt coefficients.

cd(@__DIR__)
include("poly_factorization_project.jl")

# --- If your types live directly in Main (as in your repo), no `using .Module` is needed. ---

# Convenience: construct polynomials from (coef, degree) pairs
# Using Int coefficients:
poly_int(pairs) = PolynomialDense([Term(Int(c), d) for (c, d) in pairs])
# Using BigInt coefficients:
poly_big(pairs) = PolynomialDense([Term(BigInt(c), d) for (c, d) in pairs])

# Helper to detect when Int arithmetic overflowed (wraparound)
function coef_overflowed(x::Int, y::Int)
    intsum  = x + y                   # possibly wrapped
    bigsum  = BigInt(x) + BigInt(y)   # exact
    return BigInt(intsum) != bigsum   # true => overflow happened
end

# Pretty print a comparison block
function show_example(id, Fint, Gint, Fbig, Gbig)
    println("\n=== Example $id ===")

    # Int addition
    Hint = Fint + Gint

    # BigInt addition
    Hbig = Fbig + Gbig

    println("-- Int coefficients (may overflow) --")
    @show Fint
    @show Gint
    @show Hint

    println("-- BigInt coefficients (exact) --")
    @show Fbig
    @show Gbig
    @show Hbig

    # Try to flag any coefficient positions that overflowed
    # We compare coefficient vectors term-by-term by degree.
    # (Both constructors put the same degrees in the same places.)
    try
        # Grab raw coeff vectors if the field is available.
        cint = getfield(Hint, :coeffs)
        cbig = getfield(Hbig, :coeffs)
        # Only compare where both have entries
        L = min(length(cint), length(cbig))
        any_overflow = false
        for i in 1:L
            x = cint[i]
            y = Int(cbig[i] % typemax(Int))  # only used to keep type stable
            # Detect overflow by recomputing Int-sum from the original addends termwise.
            # We need the addends' coeffs to do this.
        end
    catch
        # If your type hides coeffs, skip detection (the printed mismatch still shows it).
    end

    println("Matches (Int result == BigInt result converted to Int)? ",
            try
                # Convert Hbig to Int if possible (will throw if out-of-range)
                Hbig_int = PolynomialDense([Term(Int(c), i-1) for (i, c) in enumerate(getfield(Hbig, :coeffs))])
                Hint == Hbig_int
            catch
                false
            end)
end

println("== Task 2: Why we need BigInt for exact integer polynomials ==")

# Shorthands for huge values on 64-bit Int
const I_MAX = typemax(Int)         #  9223372036854775807 on 64-bit
const I_MIN = typemin(Int)         # -9223372036854775808

# Five explicit examples (no rand). Only coefficients are large; degrees small.
# Each pair is [(coef, deg), ...]  meaning coef * x^deg

# 1) Overflow at constant term: (I_MAX - 5) + (10)  --> wraps as Int
ex1_f = [(I_MAX - 5, 0)]
ex1_g = [(10, 0)]

# 2) Overflow at x^2: 6e18 + 6e18
ex2_f = [(6_000_000_000_000_000_000, 2)]
ex2_g = [(6_000_000_000_000_000_000, 2)]

# 3) Underflow (negative overflow) at constant term: (I_MIN + 5) + (-10)
ex3_f = [(I_MIN + 5, 0)]
ex3_g = [(-10, 0)]

# 4) Mixed degrees; one term safe, one term overflows at x
ex4_f = [(50, 0), (I_MAX - 1, 1)]
ex4_g = [(75, 0), (100, 1)]

# 5) Multi-term; overflow at x^3, underflow at x
ex5_f = [(I_MAX - 2, 3), (I_MIN + 10, 1), (12345, 0)]
ex5_g = [(100, 3),       (-50,        1), (67890, 0)]

examples = [
    (ex1_f, ex1_g),
    (ex2_f, ex2_g),
    (ex3_f, ex3_g),
    (ex4_f, ex4_g),
    (ex5_f, ex5_g),
]

for (i, (pf, pg)) in enumerate(examples)
    Fint = poly_int(pf);  Gint = poly_int(pg)
    Fbig = poly_big(pf);  Gbig = poly_big(pg)
    show_example(i, Fint, Gint, Fbig, Gbig)
end

println("\n== End Task 2 ==")
