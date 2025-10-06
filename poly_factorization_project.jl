#############################################################################
#############################################################################
#
# This is the main project file for polynomial factorization
#                                                                               
#############################################################################
#############################################################################

using Distributions, StatsBase, Random, Primes, StaticArrays

import Base: %, gcd
import Base: push!, pop!, iszero, show, isless, map, map!, iterate, length, last, isempty
import Base: +, -, *, mod, %, ÷, ==, ^, rand, div, rem, zero, one

include("src/term.jl")

# Utilities
include("src/utils/general_alg.jl")
include("src/utils/sample_primes.jl")
include("src/utils/heap.jl")
include("src/utils/misc.jl")

# Polynomial
include("src/polynomial_definitions/polynomial.jl")
    include("src/basic_polynomial_operations/abstract/polynomial_addition.jl")
    include("src/basic_polynomial_operations/abstract/polynomial_multiplication.jl")
    include("src/basic_polynomial_operations/abstract/polynomial_division.jl")
    include("src/basic_polynomial_operations/abstract/polynomial_gcd.jl")
    include("src/basic_polynomial_operations/abstract/polynomial_pseudo_operations.jl")
    include("src/basic_polynomial_operations/abstract/polynomial_field_ops.jl")

# Dense
include("src/polynomial_definitions/polynomial_dense.jl")
    include("src/basic_polynomial_operations/dense/polynomial_addition.jl")
    include("src/basic_polynomial_operations/dense/polynomial_multiplication.jl")
    include("src/basic_polynomial_operations/dense/polynomial_division.jl")
    include("src/basic_polynomial_operations/dense/polynomial_gcd.jl")

include("src/polynomial_factorization/factor.jl")
include("src/utils/pretty_print.jl")

nothing
