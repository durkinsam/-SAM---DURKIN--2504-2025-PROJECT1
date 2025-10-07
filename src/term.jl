#############################################################################
#############################################################################
#
# This file defines the Term type with several operations 
#                                                                               
#############################################################################
#############################################################################

##################################
# Term type and utility function #
##################################

"""
A term.

For Term{C, D},
    C: The type of the coefficient of the term.
    D: The type of the degree of the term.
"""
struct Term{C <: Integer, D <: Integer} #structs are immutable by default
    coeff::C
    degree::D

    function Term{C, D}(coeff::C, degree::D) where {C, D}
        degree < 0 && error("Degree must be non-negative")
        new(coeff, degree)
    end
end

######################
# Outer Constructors #
######################

""" Convenience outer constructor to infer coefficient/degree types. """
function Term(coeff::C, degree::D) where {C, D}
    Term{C, D}(coeff, degree)
end

"""
Creates the zero term.
"""
function Base.zero(::Type{Term{C, D}})::Term{C, D} where {C, D} 
    Term(zero(C), zero(D))
end
Base.zero(t::Term) = zero(typeof(t))

"""
Creates the unit term.
"""
function Base.one(::Type{Term{C, D}})::Term{C, D} where {C, D}
    Term(one(C), zero(D))
end
Base.one(t::Term) = one(typeof(t))

###########
# Display #
###########

"""
Show a term.
"""
Base.show(io::IO, t::Term) = print(io, "$(t.coeff)⋅x^$(t.degree)")

########################
# Queries about a term #
########################

"""
Check if a term is 0.
"""
Base.iszero(t::Term)::Bool = iszero(t.coeff)

"""
Compare two terms (with the same coefficient/degree types).
"""
function Base.isless(t1::Term{C, D}, t2::Term{C, D})::Bool  where {C, D}
    t1.degree == t2.degree ? (t1.coeff < t2.coeff) : (t1.degree < t2.degree)  
end

function Base.:(==)(t1::Term{C, D}, t2::Term{C, D}) where {C, D}
    t1.coeff == t2.coeff && t1.degree == t2.degree
end

"""
Evaluate a term at a point x.
Note: if C is ZModP, you should also pass x as ZModP for field evaluation.
"""
function evaluate(t::Term, x::S) where {S <: Number} 
    t.coeff * x^t.degree
end

##########################
# Operations with a term #
##########################

"""
Add two terms of the same degree (with the same coefficient/degree types).
"""
function Base.:+(t1::Term{C, D}, t2::Term{C, D})::Term{C, D} where {C, D}
    @assert t1.degree == t2.degree
    Term(t1.coeff + t2.coeff, t1.degree)
end

"""
Negate a term.
"""
function Base.:-(t::Term{C, D},)::Term{C, D} where {C, D} 
    Term(-t.coeff,t.degree)  
end

"""
Subtract two terms with the same degree (and the same coefficient/degree types).
"""
function Base.:-(t1::Term{C, D}, t2::Term{C, D})::Term{C, D} where {C, D}
    t1 + (-t2) 
end

"""
Multiply two terms (with the same coefficient/degree types).
"""
function Base.:*(t1::Term{C, D}, t2::Term{C, D})::Term{C, D} where {C, D}
    Term(t1.coeff * t2.coeff, t1.degree + t2.degree)
end

"""
Multiply a term by a constant.
"""
function Base.:*(t::Term{C, D}, n::S)::Term{C, D} where {C, D, S <: Integer}
    Term(t.coeff * n, t.degree)
end
function Base.:*(n::S, t::Term{C, D})::Term{C, D} where {C, D, S <: Integer}
    t * n
end

"""
Power of a term.
"""
function Base.:^(t::Term{C, D}, n::S)::Term{C, D} where {C, D, S <: Integer}
    Term(t.coeff^n, t.degree*D(n))
end

"""
Compute the symmetric mod of a term with an integer.
(Only for integer coefficient types; ZModP should not use this.)
"""
function Base.mod(t::Term{C, D}, p::Integer)::Term{C, D} where {C <: Integer, D}
    Term(mod(t.coeff, p), t.degree)
end

"""
Compute the derivative of a term.
"""
function derivative(t::Term{C, D})::Term{C, D} where {C, D}  
    Term{C, D}(t.coeff*C(t.degree), max(t.degree-one(D), zero(D)))
end

"""
Exact division when the coefficient type `C` is a field (e.g. Zp).

Default placeholder (will be overridden for ZModP in Task 5).
"""
function Base.div(t1::Term{C, D}, t2::Term{C, D}) where {C, D} # Base.:÷
    not_yet_implemented_error(t1, "div")
end

# ---------------------------------------------------------------------------
# --- Task 5: ZModP support (exact division + div_mod_p via conversion) -----
# ---------------------------------------------------------------------------

# Bring the ZModP type into scope (assumes src/z_mod_p.jl defines module ZModPField)
import .ZModPField: ZModP

"""
Exact division for terms with ZModP coefficients.
"""
function Base.div(t1::Term{ZModP{T,N},D}, t2::Term{ZModP{T,N},D}) where {T<:Integer,N,D}
    iszero(t2.coeff) && throw(DivideError())
    Term{ZModP{T,N},D}(t1.coeff ÷ t2.coeff, t1.degree - t2.degree)
end

"""
Refactored div_mod_p: convert to ZModP, divide exactly, convert back.
"""
function div_mod_p(t1::Term{C,D}, t2::Term{C,D}, prime::Integer) where {C<:Integer,D}
    N = Int(prime)
    Tz = ZModP{C,N}
    a = Term{Tz,D}(Tz(t1.coeff), t1.degree)
    b = Term{Tz,D}(Tz(t2.coeff), t2.degree)
    r = div(a, b)
    return Term{C,D}(C(r.coeff.val), r.degree)
end

"""
Integer divide a term by an integer (helper forwards to the above).
"""
function div_mod_p(t::Term{C, D}, n::Integer, prime::Integer) where {C <: Integer, D} 
    return div_mod_p(t, Term(C(n), zero(D)), prime)
end

#############################
# Vectorization with a term #
#############################

""" Enable broadcasting an operation on a term """
Base.broadcastable(t::Term) = Ref(t)
