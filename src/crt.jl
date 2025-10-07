###############################
# Task 7 — Chinese Remainder #
###############################

# Normalize to least nonnegative representative
@inline _norm_mod(x, m) = ((x % m) + m) % m

"""
    int_crt(rems::Vector{T}, moduli::Vector{T}) where {T<:Union{Int,BigInt}}

Solve the system of congruences:
    x ≡ rems[i] (mod moduli[i])  for i=1..k,
assuming ALL moduli are pairwise coprime.

Returns a `BigInt` solution x in the range 0 ≤ x < M where M = ∏ moduli[i].

Errors if:
- lengths mismatch
- any modulus ≤ 0
- moduli are not pairwise coprime
"""
function int_crt(rems::Vector{T}, moduli::Vector{T}) where {T<:Union{Int,BigInt}}
    length(rems) == length(moduli) || throw(ArgumentError("rems and moduli must have same length"))
    k = length(moduli)
    k > 0 || throw(ArgumentError("need at least one congruence"))

    # promote to BigInt to avoid overflow
    r = BigInt.(_norm_mod.(rems, moduli))
    m = BigInt.(moduli)

    any(x -> x ≤ 0, m) && throw(ArgumentError("all moduli must be positive"))
    # pairwise coprime check (O(k^2) but k is small here)
    for i in 1:k, j in i+1:k
        gcd(m[i], m[j]) == 1 || throw(ArgumentError("moduli must be pairwise coprime"))
    end

    # iterative CRT: combine one modulus at a time
    x = r[1]
    M = m[1]
    for i in 2:k
        # solve x ≡ x (mod M) and x ≡ r[i] (mod m[i])
        # find t ≡ (r[i] - x) * inv(M, m[i]) (mod m[i])
        Mi = M % m[i]
        invMi = invmod(Mi, m[i])                     # exists because gcd(M, m[i]) = 1
        t = _norm_mod((r[i] - (x % m[i])) * invMi, m[i])
        x += M * t
        M *= m[i]
        x = _norm_mod(x, M)                          # keep x small
    end
    return x  # BigInt in [0, M-1]
end

"""
    verify_crt(x, rems, moduli)

Return a Bool vector indicating which congruences x satisfies.
"""
function verify_crt(x::Integer, rems::Vector, moduli::Vector)
    [ ((x % m) + m) % m == ((r % m) + m) % m for (r, m) in zip(rems, moduli) ]
end
