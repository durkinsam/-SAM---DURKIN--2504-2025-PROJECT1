# src/utils/pretty_print.jl
# Try common ways to fetch coefficients (a0, a1, a2, ... low to high)
function coeffvec(p)
    if hasproperty(p, :coeffs)
        return getproperty(p, :coeffs)     # p.coeffs
    else
        try
            return coefficients(p)         # if a method exists
        catch
            error("pretty_print: can't locate coefficients for $(typeof(p)).")
        end
    end
end

# Build a display string that follows all Task 1 pretty-print rules
function _poly_str(p)
    cs = coeffvec(p)
    n  = length(cs) - 1
    terms = String[]

    for k in n:-1:0
        a = cs[k+1]
        a == 0 && continue

        sgn = a < 0 ? " - " : (isempty(terms) ? "" : " + ")
        c   = abs(a)

        t = if k == 0
            string(c)                        # constant: no x^0
        elseif k == 1
            c == 1 ? "x" : string(c, "x")    # degree 1: no ^1; hide unit coeff
        else
            c == 1 ? "x^$k" : string(c, "x^", k)  # hide unit coeff
        end

        push!(terms, sgn * t)
    end

    isempty(terms) && return "0"
    s = join(terms, "")
    startswith(s, " - ") ? s[2:end] : s
end

import Base: show
function show(io::IO, ::MIME"text/plain", p)
    print(io, _poly_str(p))
end
function show(io::IO, p)
    print(io, _poly_str(p))
end
