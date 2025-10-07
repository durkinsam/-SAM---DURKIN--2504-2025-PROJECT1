# example_script_task_3.jl
# Task 3: Check sparse works + show Sparse outperforms Dense (time & memory).
# Uses only non-overflowing cases (BigInt coeffs, Int degrees).

cd(@__DIR__)
include("poly_factorization_project.jl")  # <- whatever loads your package
# Append Task 3.8 code (your sparse ops overrides)
# If you split ops, include them all; otherwise include the single ops file.
include("src/basic_polynomial_operations/sparse/polynomial_sparse_ops.jl")

using BenchmarkTools, Printf

# ── builders ────────────────────────────────────────────────────────────────
makeDense(::Type{C}, ::Type{D}, pairs) where {C,D} =
    PolynomialDense{C,D}(Term{C,D}[Term{C,D}(convert(C,c), convert(D,d)) for (c,d) in pairs])

makeSparse(::Type{C}, ::Type{D}, pairs) where {C,D} =
    PolynomialSparse{C,D}(Term{C,D}[Term{C,D}(convert(C,c), convert(D,d)) for (c,d) in pairs])

coeffvec(p) = begin
    ld = leading(p); deg = Int(ld.degree); C = typeof(ld.coeff)
    v = [zero(C) for _ in 0:deg]
    for t in p; v[Int(t.degree)+1] = t.coeff; end
    v
end

# ── pretty print a benchmark comparison ────────────────────────────────────
function compare!(label, f_dense, f_sparse)
    bd = @benchmark $f_dense()
    bs = @benchmark $f_sparse()
    td = minimum(bd).time / 1e6;  md = minimum(bd).memory
    ts = minimum(bs).time / 1e6;  ms = minimum(bs).memory
    @printf("\n%s\n", label)
    @printf("Dense : %8.3f ms  %10d bytes\n", td, md)
    @printf("Sparse: %8.3f ms  %10d bytes\n", ts, ms)
    @printf("Speedup (Dense/Sparse): %.2fx   Mem ratio: %.2fx\n", td/ts, md/max(ms,1))
end

println("\n=== Task 3: PolynomialSparse vs PolynomialDense (no overflow cases) ===")

C = BigInt; D = Int

# Choose a high degree where dense is disadvantaged but still reasonable to run.
const N = 100_000

# x^N - 1 and x^N + 1 (sparse: 2 terms; dense: N+1 terms)
p_dense = makeDense(C,D, [(1,N), (-1,0)])
q_dense = makeDense(C,D, [(1,N), ( 1,0)])
p_sparse = makeSparse(C,D, [(1,N), (-1,0)])
q_sparse = makeSparse(C,D, [(1,N), ( 1,0)])

# Term for term-multiplication (shift by k)
t_dense = Term{C,D}(one(C), D(N÷2))
t_sparse = Term{C,D}(one(C), D(N÷2))

# Sanity: results match exactly
@assert coeffvec(p_dense + q_dense) == coeffvec(p_sparse + q_sparse)
@assert coeffvec(t_dense * p_dense) == coeffvec(t_sparse * p_sparse)
@assert coeffvec(derivative(p_dense)) == coeffvec(derivative(p_sparse))
@assert coeffvec(mod(p_dense, 101)) == coeffvec(mod(p_sparse, 101))

# 1) Addition: (x^N - 1) + (x^N + 1) = 2 x^N
compare!("Addition: (x^N - 1) + (x^N + 1)",
    () -> (p_dense + q_dense),
    () -> (p_sparse + q_sparse)
)

# 2) Term multiplication: x^(N÷2) * (x^N - 1)
compare!("Term multiplication: x^(N÷2) * (x^N - 1)",
    () -> (t_dense * p_dense),
    () -> (t_sparse * p_sparse)
)

# 3) Derivative: d/dx (x^N - 1) = N x^(N-1)
compare!("Derivative: d/dx (x^N - 1)",
    () -> derivative(p_dense),
    () -> derivative(p_sparse)
)

# 4) Mod prime: coefficients mod 101
compare!("mod 101: reduce coefficients (structure preserved)",
    () -> mod(p_dense, 101),
    () -> mod(p_sparse, 101)
)

println("\n=== End Task 3 ===\n")

# Hints for your report:
# Run and capture:
#   julia --project=. example_script_task_3.jl | tee sparse_vs_dense_output.txt
# Add a small table in your PDF:
#   Pros of Sparse: memory for high-degree with few terms; faster on +/t* p; storage matches structure.
#   Cons of Sparse: worse for very dense polys; random access by degree is O(log n)/O(n); more merging logic.
