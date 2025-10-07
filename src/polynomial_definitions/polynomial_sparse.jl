#############################################################################
# Sparse polynomial operations for PolynomialSparse{C,D}
# Optimised for:
#   • p + q         : linear merge in degree order
#   • p + t (term)  : linear insert/merge
#   • t * p         : linear map, preserves sorted unique degrees
# Also provides sparse-aware: mod(·, prime), derivative(·)
#############################################################################

# ========== helpers (local; rely only on the sparse representation) ==========
@inline _iszeroC(::Type{C}, x) where {C} = x == zero(C)

# Merge two ascending-degree term streams into a new term vector
function _merge_add(ps::Vector{Term{C,D}}, qs::Vector{Term{C,D}}) where {C,D}
    vt = Term{C,D}[]
    sizehint!(vt, length(ps) + length(qs))
    i = 1; j = 1
    while i <= length(ps) || j <= length(qs)
        if i > length(ps)
            push!(vt, qs[j]); j += 1
        elseif j > length(qs)
            push!(vt, ps[i]); i += 1
        else
            ti = ps[i]; tj = qs[j]
            if ti.degree == tj.degree
                c = ti.coeff + tj.coeff
                iszero(c) || push!(vt, Term{C,D}(c, ti.degree))
                i += 1; j += 1
            elseif ti.degree < tj.degree
                push!(vt, ti); i += 1
            else
                push!(vt, tj); j += 1
            end
        end
    end
    return vt
end

# Insert/merge a single term t (ascending degree) into p.terms (ascending)
function _insert_add(ps::Vector{Term{C,D}}, t::Term{C,D}) where {C,D}
    iszero(t) && return copy(ps)
    vt = Term{C,D}[]
    sizehint!(vt, length(ps) + 1)
    inserted = false
    for u in ps
        if !inserted && t.degree < u.degree
            push!(vt, t); push!(vt, u); inserted = true
        elseif !inserted && t.degree == u.degree
            c = u.coeff + t.coeff
            iszero(c) || push!(vt, Term{C,D}(c, u.degree))
            inserted = true
        else
            push!(vt, u)
        end
    end
    inserted || push!(vt, t)
    return vt
end

# ========== Addition (optimised) ============================================

# p + q  (PolynomialSparse + PolynomialSparse) — linear merge
function Base.:+(p::PolynomialSparse{C,D}, q::PolynomialSparse{C,D})::PolynomialSparse{C,D} where {C,D}
    vt = _merge_add(p.terms, q.terms)
    return PolynomialSparse{C,D}(vt)  # already sparse/ordered
end

# p + t  (PolynomialSparse + Term) — linear insert/merge
function Base.:+(p::PolynomialSparse{C,D}, t::Term{C,D})::PolynomialSparse{C,D} where {C,D}
    vt = _insert_add(p.terms, t)
    return PolynomialSparse{C,D}(vt)
end
Base.:+(t::Term{C,D}, p::PolynomialSparse{C,D}) where {C,D} = p + t

# p + n  (add constant; use degree 0 of type D, coefficient type C)
Base.:+(p::PolynomialSparse{C,D}, n::Integer) where {C,D} =
    p + Term{C,D}(convert(C,n), zero(D))
Base.:+(n::Integer, p::PolynomialSparse{C,D}) where {C,D} = p + n

# ========== Term multiplication (optimised) =================================
# t * p — map once, keep order (degrees strictly increase so no merging needed)
function Base.:*(t::Term{C,D}, p::PolynomialSparse{C,D})::PolynomialSparse{C,D} where {C,D}
    iszero(t) && return PolynomialSparse{C,D}()
    vt = Term{C,D}[]
    sizehint!(vt, length(p.terms))
    @inbounds for u in p.terms
        iszero(u) && continue
        c = t.coeff * u.coeff
        iszero(c) && continue
        push!(vt, Term{C,D}(c, u.degree + t.degree))
    end
    return PolynomialSparse{C,D}(vt)
end
Base.:*(p::PolynomialSparse{C,D}, t::Term{C,D}) where {C,D} = t * p

# (optional) scalar-on-either-side
Base.:*(α::Number, p::PolynomialSparse{C,D}) where {C,D} =
    (iszero(α) ? PolynomialSparse{C,D}() :
     PolynomialSparse{C,D}([Term{C,D}(convert(C,α)*u.coeff, u.degree) for u in p.terms]))
Base.:*(p::PolynomialSparse{C,D}, α::Number) where {C,D} = α * p

# ========== Derivative (sparse-aware) =======================================
function derivative(p::PolynomialSparse{C,D})::PolynomialSparse{C,D} where {C,D}
    vt = Term{C,D}[]
    sizehint!(vt, max(0, length(p.terms)-1))
    @inbounds for u in p.terms
        iszero(u.degree) && continue
        # multiply coefficient by degree (cast degree -> C), then drop degree by 1
        c = u.coeff * convert(C, u.degree)
        iszero(c) && continue
        push!(vt, Term{C,D}(c, u.degree - one(D)))
    end
    return PolynomialSparse{C,D}(vt)
end

# ========== Mod prime (sparse-aware; no indexing) ============================
# Reduce coefficients modulo a prime (keep non-zero residues only)
function Base.mod(f::PolynomialSparse{C,D}, p::Int)::PolynomialSparse{C,D} where {C,D}
    vt = Term{C,D}[]
    sizehint!(vt, length(f.terms))
    @inbounds for t in f.terms
        r = mod(t.coeff, p)
        iszero(r) || push!(vt, Term{C,D}(convert(C, r), t.degree))
    end
    return PolynomialSparse{C,D}(vt)
end

# ========== Equality fast path (vector compare is fine) ======================
Base.:(==)(p::PolynomialSparse{C,D}, q::PolynomialSparse{C,D}) where {C,D} =
    p.terms == q.terms
