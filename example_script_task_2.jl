# example_script_task_2.jl
# Task 2: Parameterising coefficient type and showing why we need BigInt.
# We construct 5 explicit polynomial-addition examples whose coefficients overflow
# if stored as Int, then re-run the same additions with BigInt coefficients.

# Load your project code (adjust this path/name if yours differs)
cd(@__DIR__)
include("poly_factorization_project.jl")

# ──────────────────────────────────────────────────────────────────────────────
# Utilities (generic; use only the public interface: Term, PolynomialDense, +)
# ──────────────────────────────────────────────────────────────────────────────

# Build a PolynomialDense{C,D} from (coef,deg) pairs (deg are plain Ints)
makeP(::Type{C}, ::Type{D}, pairs::Vector{Tuple{<:Integer,<:Integer}}) where {C,D} =
    PolynomialDense{C,D}(Term{C,D}[ Term{C,D}(convert(C,c), convert(D,d)) for (c,d) in pairs ])

# Materialize a coefficient vector (a₀, a₁, …, a_deg) from any Polynomial{C,D}
function coeffvec(p)::Vector
    # Determine length from leading degree
    ld = leading(p)                  # Term{C,D}
    D  = typeof(ld.degree)
    C  = typeof(ld.coeff)
    deg = Int(ld.degree)
    v = [zero(C) for _ in 0:deg]
    # Fill from the iterator of non-zero terms
    for t in p
        v[Int(t.degree)+1] = t.coeff
    end
    return v
end

# Lift a coefficient vector to BigInt for safe comparison
bigvec(v) = BigInt[ big(x) for x in v ]

# Pretty printer for one example
function show_example(label::String, Pint, Qint, Pbig, Qbig)
    println("\n", "─"^80)
    println(label)
    println("─"^80)

    # Compute sums
    Rint = Pint + Qint
    Rbig = Pbig + Qbig

    # Show polynomials
    println("[Int]  p(x) = ", Pint)
    println("[Int]  q(x) = ", Qint)
    println("[Int]  p(x)+q(x) = ", Rint, "   (⚠ may overflow)")

    println("[Big]  p(x) = ", Pbig)
    println("[Big]  q(x) = ", Qbig)
    println("[Big]  p(x)+q(x) = ", Rbig, "   (✓ exact)")

    # Compare coefficient vectors safely (lift Int→BigInt)
    cint = coeffvec(Rint)
    cbig = coeffvec(Rbig)
    cint_big = bigvec(cint)

    # Align lengths
    L = max(length(cint_big), length(cbig))
    resize!(cint_big, L; init=big(0))
    resize!(cbig, L; init=zero(eltype(cbig)))

    # Detect mismatch == evidence of overflow/wrap
    mismatch_positions = Int[]
    for i in 1:L
        if cint_big[i] != cbig[i]
            push!(mismatch_positions, i-1)  # degree = index-1
        end
    end

    println("[Int] coeffs: ", cint)
    println("[Big] coeffs: ", cbig)
    println("Overflow detected at degrees: ",
            isempty(mismatch_positions) ? "none" : string(mismatch_positions))
end

# ──────────────────────────────────────────────────────────────────────────────
# Examples (explicit, no randomness). Only coefficients are large; degrees small
# ──────────────────────────────────────────────────────────────────────────────

const I_MAX = typemax(Int)      #  9223372036854775807 on 64-bit
const I_MIN = typemin(Int)      # -9223372036854775808

# Each example is a pair of vectors of (coef, degree)
ex1_f = [(I_MAX - 5, 0)];                   ex1_g = [(10, 0)]                       # overflow at constant
ex2_f = [(6_000_000_000_000_000_000, 2)];   ex2_g = [(6_000_000_000_000_000_000, 2)]# overflow at x^2
ex3_f = [(I_MIN + 5, 0)];                   ex3_g = [(-10, 0)]                      # negative overflow
ex4_f = [(50, 0), (I_MAX - 1, 1)];          ex4_g = [(75, 0), (100, 1)]             # mix; overflow at x
ex5_f = [(I_MAX - 2, 3), (I_MIN + 10, 1), (12345, 0)]
ex5_g = [(100,       3), (-50,        1), (67890, 0)]                                # overflow at x^3 & x

examples = [
    ("Example 1 — positive overflow at constant term", ex1_f, ex1_g),
    ("Example 2 — large equal terms overflow at x^2",  ex2_f, ex2_g),
    ("Example 3 — negative overflow at constant",      ex3_f, ex3_g),
    ("Example 4 — mixed degrees; overflow at x",       ex4_f, ex4_g),
    ("Example 5 — multi-term; wraps at x^3 and x",     ex5_f, ex5_g),
]

println("\n=== Task 2: Int vs BigInt overflow demonstration ===")

for (label, pf, pg) in examples
    # Int,Int polynomials (will wrap)
    Fint = makeP(Int,    Int, pf)
    Gint = makeP(Int,    Int, pg)
    # BigInt,Int polynomials (exact)
    Fbig = makeP(BigInt, Int, pf)
    Gbig = makeP(BigInt, Int, pg)
    show_example(label, Fint, Gint, Fbig, Gbig)
end

# A pure BigInt scale example (orders of magnitude beyond Int)
println("\n", "─"^80)
println("Example 6 — truly huge BigInt coefficients (no Int counterpart)")
println("─"^80)
Pbig = makeP(BigInt, Int, [(big"10"^40,0), (-big"10"^35,1), (big"10"^30,2)])
Qbig = makeP(BigInt, Int, [(big"10"^40,0), ( big"10"^35,1), (-big"10"^25,2)])
Rbig = Pbig + Qbig
println("[Big] P(x)+Q(x) = ", Rbig)
println("[Big] coeffs: ", coeffvec(Rbig))

println("\n=== End Task 2 ===\n")

# Tip to capture output for PDF:
# julia --project=. example_script_task_2.jl | tee overflow_demo_output.txt
