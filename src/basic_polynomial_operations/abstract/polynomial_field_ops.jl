#############################################################################
#############################################################################
#
# This file implements operations for polynomials over a field, e.g. Zp[x].
#                                                                               
#############################################################################
#############################################################################

#= TASK LIST:
NOTE - the following functions will NOT work when C is an integer until you 
finish implementing all tasks (and even then, some functions such as gcd simply
cannot be implemented over the ring of integers - see bonus task 3 for details).

TODO (Task 2) done
    In this task, simply edit this file such that the first argument (::Type{C})
    to each function is also the type of the coefficients of the polynomials.
    I.e., we want the type signatures to contain:
        {C, D, P <: Polynomial{C, D}}

TODO (Task 6)
    In this task, override the unimplemented functions (after duplicating this 
    file as per the instructions) for C <: ZModP. When overriding, use the
    corresponding _mod_p function as a guide.

TODO (Task 7)
    Here, (in the duplicated file as per the instructions) override the factor
    function for C <: Integer and implement the Chinese Remainder Theorem (CRT).

=#


"""
Returns the factors of `f` in an array of tuples (g, n).
NOTE: Override this in Task 6 for Zp[x] and in Task 7 for Z[x].
"""
function factor(::Type{C}, f::P)::Vector{Tuple{P, Integer}} where {C,D,P<:Polynomial{C,D}}
    not_implemented_error(f, "factor")
end

"""
Returns (q, r) for num ÷ den.
NOTE: Override this in Task 6 for Zp[x].
"""
function div_rem(::Type{C}, num::P, den::P)::Tuple{P,P} where {C,D,P<:Polynomial{C,D}}
    not_implemented_error(num, "div_rem")
end

"""
Distinct degree factorization.
NOTE: Override this in Task 6 for Zp[x].
"""
function dd_factor(::Type{C}, f::P)::Array{P} where {C,D,P<:Polynomial{C,D}}
    not_implemented_error(f, "dd_factor")
end

"""
Distinct degree split.
NOTE: Override this in Task 6 for Zp[x].
"""
function dd_split(::Type{C}, f::P, d::Integer)::Vector{P} where {C,D,P<:Polynomial{C,D}}
    not_implemented_error(f, "dd_split")
end

"""
Quotient only.
"""
div(::Type{C}, num::P, den::P) where {C,D,P<:Polynomial{C,D}} = first(div_rem(C, num, den))

"""
Remainder only.
"""
rem(::Type{C}, num::P, den::P) where {C,D,P<:Polynomial{C,D}} = last(div_rem(C, num, den))

"""
Extended Euclid (generic over a field C).
"""
function extended_euclid_alg(::Type{C}, f::P, g::P) where {C,D,P<:Polynomial{C,D}}
    return ext_euclid_alg(f, g, rem, div)
end

"""
Greatest common divisor.
"""
gcd(::Type{C}, f::P, g::P) where {C,D,P<:Polynomial{C,D}} =
    extended_euclid_alg(C, f, g) |> first

"""
Yun's algorithm: square-free part over a perfect field C.
"""
function square_free(::Type{C}, f::P) where {C,D,P<:Polynomial{C,D}}
    # Remove minimum degree (in case char(C) != 0)
    min_deg = last(f).degree
    vt = filter(t -> !iszero(t), collect(f))
    # Ensure we rebuild with Term{C,D} so degrees/coeffs keep their param types
    shifted_terms = map(t -> Term{C,D}(t.coeff, t.degree - min_deg), vt)
    f = P(shifted_terms)

    # Compute the gcd of f, f'
    der_f = derivative(f)
    sqr_part = gcd(C, f, der_f)

    iszero(sqr_part) && return f * (min_deg > zero(min_deg) ? x_poly(P) : one(P))

    # Remove factors with multiplicity > 1
    sqr_free = div(f, sqr_part)

    # Add one factor of x back in if necessary
    if min_deg > zero(min_deg)
        sqr_free *= x_poly(P)
    end

    return sqr_free
end

"""
Multiplicity of g in f.
"""
function multiplicity(::Type{C}, f::P, g::P)::Integer where {C,D,P<:Polynomial{C,D}}
    degree(gcd(C, f, g)) == 0 && return 0
    return 1 + multiplicity(C, div(f, g), g)
end
