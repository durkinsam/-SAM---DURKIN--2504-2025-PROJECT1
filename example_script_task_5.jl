# example_script_task_5.jl
cd(@__DIR__)
include("poly_factorization_project.jl")   # or your package entry
using .ZModPField

println("\n=== Task 5: ZModP arithmetic demo ===")
const N = 7
A = ZModP{Int,N}(3)
B = ZModP{Int,N}(5)

println("a = ", A, "   b = ", B)
println("a + b = ", A + B)            # 3+5=8≡1
println("a - b = ", A - B)            # 3-5=-2≡5
println("a * b = ", A * B)            # 15≡1
println("inv(a) = ", inv(A))          # 3⁻¹ ≡ 5 mod 7
println("a ÷ b = ", A ÷ B)            # 3 * 5⁻¹
println("2 + a = ", 2 + A)
println("a + 2 = ", A + 2)
println("2 * a = ", 2 * A)
println("a * 2 = ", A * 2)
println("a^10  = ", A^10)
println("abs(a) = ", abs(A))
println("Int(a) = ", Int(A), "  BigInt(b) = ", BigInt(B))
println("a == 3?  ", A == 3, "   3 == a? ", 3 == A)
println("=== End Task 5 ===\n")

