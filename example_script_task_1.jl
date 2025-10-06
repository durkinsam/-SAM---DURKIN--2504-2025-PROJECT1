cd(@__DIR__)                     # ensures relative paths work
include("poly_factorization_project.jl")  # loads the main module
# example_script_task_1.jl
# Demonstrates construction, basic ops, derivative product rule,
# modular division, and gcd over F_p.


println("== Task 1 example script ==")

# 1) Construct at least 3 polynomials with the given constructor.
#    In most templates, PolynomialDense takes a coefficient vector
# f(x) = 1 - 2x + 3x^2
f = PolynomialDense([Term(1,0), Term(-2,1), Term(3,2)])

# g(x) = 5x + x^2
g = PolynomialDense([Term(5,1), Term(1,2)])

# h(x) = 7 + x^3
h = PolynomialDense([Term(7,0), Term(1,3)])


@show f
@show g
@show h

println("\n-- Basic operations --")
@show f + g
@show f * h
@show g * h

println("\n-- Derivative product rule check: d(fg) = f'g + fg' --")
lhs = derivative(f * g)
rhs = derivative(f) * g + f * derivative(g)
@show lhs
@show rhs
@show lhs == rhs

println("\n-- Modular division: (f*h) ÷ h over F_p --")
# Your repo almost certainly has a function called `div_mod_p` that returns (q, r).
# If your example_script.jl shows a different call (e.g. `div_mod_p!` or a method `÷(a,b,p)`),
# use that exact call here instead.
for p in (5, 17, 101)
    q, r = div_mod_p(f * h, h, p)   # quotient, remainder over F_p[x]
    println("p = $p")
    @show q
    @show r
    # Confirm the quotient equals f reduced mod p:
    @show q == mod(f, p)
end

println("\n-- GCD over F_p --")
for p in (5, 11, 13)
    g1 = gcd_mod_p(f * h, g * h, p)
    println("p = $p")
    @show g1
end

println("\n== Done ==")
