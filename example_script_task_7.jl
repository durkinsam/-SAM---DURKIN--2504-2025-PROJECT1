# example_script_task_7.jl
# Task 7: Integer CRT examples (≥3 cases, each with ≥5 congruences)

cd(@__DIR__)
include("poly_factorization_project.jl")   # or your project entry file
include("src/crt.jl")

println("\n=== Task 7 — Integer CRT Demo ===")

function demo_case!(label, rems, moduli)
    println("\n-- ", label, " --")
    println("rems   = ", rems)
    println("moduli = ", moduli)
    x = int_crt(rems, moduli)
    M = prod(BigInt.(moduli))
    println("solution x = ", x, "  (mod ", M, ")")
    ok = verify_crt(x, rems, moduli)
    println("checks  = ", ok)
    println("all satisfied? ", all(ok))
end

# Case A: 5 equations, small coprime moduli
remsA   = [2, 3, 1, 4, 5]
moduliA = [3, 4, 5, 7, 11]   # pairwise coprime
demo_case!("Case A", remsA, moduliA)

# Case B: 6 equations, mixed signs in remainders
remsB   = [ -12, 7, 0, 25, -3, 19 ]
moduliB = [   5, 7, 9, 11, 13, 17 ]       # 9 is not coprime with … oh no! fix to 8
moduliB  = [5, 7, 8, 11, 13, 17]          # pairwise coprime now
demo_case!("Case B", remsB, moduliB)

# Case C: 7 equations, larger moduli (still coprime)
remsC   = [123, 456, 789, 101, 202, 303, 404]
moduliC = [23, 29, 31, 37, 41, 43, 47]
demo_case!("Case C", remsC, moduliC)

println("\n=== End Task 7 ===\n")
