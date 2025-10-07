#############################################################################
#############################################################################
#
# This file implements polynomial addition for dense polynomials.
#
#############################################################################
#############################################################################

# Assumes:
#   mutable struct PolynomialDense{C,D} <: Polynomial{C,D}
#       terms :: Vector{Term{C,D}}      # degree = index-1 (ascending), canonical zero as [0*x^0]
#   end
# and trim!(::Polynomial{C,D}) is available from the abstract layer.

"""
Add a dense polynomial and a term (functional).
Merges at degree = t.degree and trims trailing zeros.
"""
function Base.:+(p::PolynomialDense{C,D}, t::Term{C,D}) where {C,D}
    iszero(t) && return p
    q = deepcopy(p)
    idx = Int(t.degree) + 1                # array indices are Int
    if idx > length(q.terms)
        # grow with zero terms up to idx
        zt = Term{C,D}(zero(C), zero(D))
        append!(q.terms, fill(zt, idx - length(q.terms) - 0))
        push!(q.terms, t)
    else
        # merge at that slot
        if !iszero(q.terms[idx])
            q.terms[idx] = q.terms[idx] + t
        else
            q.terms[idx] = t
        end
    end
    return trim!(q)
end

# Term + PolynomialDense (commutative)
Base.:+(t::Term{C,D}, p::PolynomialDense{C,D}) where {C,D} = p + t

"""
In-place add a term (mutating).
"""
function Base.:+=(p::PolynomialDense{C,D}, t::Term{C,D}) where {C,D}
    iszero(t) && return p
    idx = Int(t.degree) + 1
    if idx > length(p.terms)
        zt = Term{C,D}(zero(C), zero(D))
        append!(p.terms, fill(zt, idx - length(p.terms) - 0))
        push!(p.terms, t)
    else
        p.terms[idx] = iszero(p.terms[idx]) ? t : (p.terms[idx] + t)
    end
    trim!(p)
end

"""
Add an integer constant to a dense polynomial.
(Converts n to coefficient type C and uses degree zero of type D.)
"""
Base.:+(p::PolynomialDense{C,D}, n::Integer) where {C,D} =
    p + Term{C,D}(convert(C, n), zero(D))
Base.:+(n::Integer, p::PolynomialDense{C,D}) where {C,D} =
    p + n

# (Optional) If you want a fast dense+dense path (otherwise abstract fallback is fine):
# function Base.:+(p::PolynomialDense{C,D}, q::PolynomialDense{C,D})::PolynomialDense{C,D} where {C,D}
#     np, nq = length(p.terms), length(q.terms)
#     L = max(np, nq)
#     zt = Term{C,D}(zero(C), zero(D))
#     vt = Vector{Term{C,D}}(undef, L)
#     @inbounds for i in 1:L
#         tp = (i <= np) ? p.terms[i] : zt
#         tq = (i <= nq) ? q.terms[i] : zt
#         vt[i] = tp + tq
#     end
#     return trim!(PolynomialDense{C,D}(vt))
# end
