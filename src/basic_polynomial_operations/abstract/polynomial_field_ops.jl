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
using .ZModPField: ZModP

# ---------- helpers: rebuild same concrete poly type with new coeff type ----
@generated function _rebuild_like(p::P, ::Type{C2}, terms::Vector{Term{C2,D}}) where {C2,D,P<:Polynomial{<:Any,D}}
    head = :( $(Expr(:quote, :(typeof(p).name.wrapper))) )
    :( $head{C2,D}(terms) )
end

# Converters between rings (used when we reuse existing *_mod_p code)
"Int/BigInt → ZModP coefficients"
function to_zp_poly(p::P, ::Type{Z}) where {C<:Integer,D,P<:Polynomial{C,D},T<:Integer,N,Z<:ZModP{T,N}}
    vt = Term{Z,D}[ Term{Z,D}(Z(t.coeff), t.degree) for t in p ]
    _rebuild_like(p, Z, vt)
end

"ZModP → Int/BigInt coefficients (choose target Int type T2)"
function to_int_poly(p::P, ::Type{T2}) where {T2<:Integer,D,T<:Integer,N,P<:Polynomial{ZModP{T,N},D}}
    vt = Term{T2,D}[ Term{T2,D}(T2(t.coeff.val), t.degree) for t in p ]
    _rebuild_like(p, T2, vt)
end

# extract the modulus N from the coefficient type
_modulus(::Type{ZModP{T,N}}) where {T,N} = N

# ======================= REQUIRED API (Task 6) ==============================

"""
div_rem(::Type{ZModP{T,N}}, num, den) over Zp[x].
Returns (q,r) with num = q*den + r and deg(r) < deg(den) or r=0.
"""
function div_rem(::Type{Cz}, num::P, den::P)::Tuple{P,P} where {T<:Integer,N,D,Cz<:ZModP{T,N},P<:Polynomial{Cz,D}}
    iszero(den) && throw(DivideError())
    q = _rebuild_like(num, Cz, Term{Cz,D}[])
    r = deepcopy(num)
    while !iszero(r) && degree(r) >= degree(den)
        lt_r = leading(r); lt_d = leading(den)
        t = Term{Cz,D}(lt_r.coeff ÷ lt_d.coeff, lt_r.degree - lt_d.degree)  # exact in Zp
        q = q + t
        r = r - (t * den)
    end
    (trim!(q), trim!(r))
end

# We can reuse your existing *_mod_p implementations by converting
# Zp[x] -> BigInt polynomials, calling the old function, then converting back.
const _Carrier = BigInt

"""
Returns the factors of f as (g, multiplicity) over Zp[x].
"""
function factor(::Type{Cz}, f::P)::Vector{Tuple{P,Int}} where {T<:Integer,N,D,Cz<:ZModP{T,N},P<:Polynomial{Cz,D}}
    Nmod = _modulus(Cz)
    f_int = to_int_poly(f, _Carrier)
    parts_int = factor_mod_p(f_int, Nmod)      # existing code
    [ (to_zp_poly(fi, Cz), m) for (fi,m) in parts_int ]
end

"""
Distinct degree factorization over Zp[x].
"""
function dd_factor(::Type{Cz}, f::P)::Array{P} where {T<:Integer,N,D,Cz<:ZModP{T,N},P<:Polynomial{Cz,D}}
    Nmod = _modulus(Cz)
    f_int = to_int_poly(f, _Carrier)
    arr_int = dd_factor_mod_p(f_int, Nmod)     # existing code
    [ to_zp_poly(g, Cz) for g in arr_int ]
end

"""
Distinct degree split over Zp[x].
"""
function dd_split(::Type{Cz}, f::P, d::Integer)::Vector{P} where {T<:Integer,N,D,Cz<:ZModP{T,N},P<:Polynomial{Cz,D}}
    Nmod = _modulus(Cz)
    f_int = to_int_poly(f, _Carrier)
    vec_int = dd_split_mod_p(f_int, d, Nmod)   # existing code
    [ to_zp_poly(g, Cz) for g in vec_int ]
end

# The following generic helpers from the abstract file will now work automatically:
#   div(::Type{C}, num, den) = first(div_rem(C, num, den))
#   rem(::Type{C}, num, den) = last(div_rem(C, num, den))
#   extended_euclid_alg, gcd, square_free, multiplicity
# because they only require field ops and div/rem defined above.

