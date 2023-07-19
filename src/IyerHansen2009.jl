module IyerHansen2009

using ArbNumerics

# utility functions of commonly repeated equations
𝒜(h, hsc) = √((1 - 2h / hsc) * (1 + 6h / hsc))
ℬ(h, hsc, c) = 1 - c * h / hsc - 𝒜(h, hsc)
𝒞(ω₀, s) = 1 + s * √(1 - ω₀)

# alias around elliptical integrals
function Π(n::T, k::T)::T where {T}
    T(elliptic_pi(ArbFloat(n), ArbFloat(k)))
end
function Π(n::T, ϕ::T, k::T)::T where {T} 
    T(elliptic_pi(ArbFloat(n), ArbFloat(ϕ), ArbFloat(k)))
end

"""
    r₀(b, a, M, s)

Calculate Eq. (20) for impact parameter `b`, spin `a`,
black hole mass `M` and sign `s`.
"""
function r₀(b, a, M, s)
    numerator = -3 * √3 * M * (1 - a / (s * b))^2
    denominator = b * (1 - (a^2 / b^2))^(3 / 2)

    ang = (1 / 3) * acos(numerator / denominator)

    (2b / √3) * √(1 - (a^2 / b^2)) * cos(ang)
end

"""
    r₀_Q(h, hsc)

Eq. (40)
"""
r₀_Q(h, hsc) = inv(hsc * 𝒜(h, hsc))

"""
    k²(h, hsc)

Eq. (41)
"""
k²(h, hsc) = -ℬ(h, hsc, 6) / (2 * 𝒜(h, hsc))

"""
    ψ(h, hsc)

Eq. (42)
"""
ψ(h, hsc) = asin(√(ℬ(h, hsc, 2) / ℬ(h, hsc, 6)))

"""
    Ω₊₋(h, hsc, ω₀, ωs, s)

Eq. (43)
"""
function Ω₊₋(h, hsc, ω₀, ωs, s)
    numerator = s * 𝒞(ω₀, s) * (1 - ωs) - s * ω₀ / 2
    denominator = √(1 - ω₀) * (𝒞(ω₀, s) - (ω₀ * hsc / 4) * ℬ(h, hsc, 2))
    numerator / denominator
end
Ω₊(h, hsc, ω₀, ωs) = Ω₊₋(h, hsc, ω₀, ωs, 1)
Ω₋(h, hsc, ω₀, ωs) = Ω₊₋(h, hsc, ω₀, ωs, -1)

"""
    n₊₋(h, hsc, ω₀)

Eq. (44)
"""
function n₊₋(h, hsc, ω₀, s)
    numerator = ℬ(h, hsc, 6)
    denominator = ℬ(h, hsc, 2) - (4 / (ω₀ * hsc)) * 𝒞(ω₀, s)
    numerator / denominator
end
n₊(h, hsc, ω₀) = n₊₋(h, hsc, ω₀, 1)
n₋(h, hsc, ω₀) = n₊₋(h, hsc, ω₀, -1)

# impact parameter limits for the shadow
b₊₋(M, a, s) = -a + s * 6M * cos((1/3) * acos(-s * a / M))
b₊(M, a) = b₊₋(M, a, 1)
b₋(M, a) = b₊₋(M, a, -1)

"""
    deflection_angle(M, a, b)

Calculate the analytic deflection angle for a geodesic to an asymptotic 
observer in the equitorial plane for a black hole with mass `M`, spin `a`, 
at impact parameter `b`.

Uses Eq. (34) to compute the semi-analytic result using elliptical integrals.
"""
function deflection_angle(M, a, b)
    # make sure impact parameter does not correspond
    # to geodesic that falls into the event horizon
    if (b₋(M, a) ≤ b ≤ b₊(M, a))
        return NaN
    end
    _deflection_angle(M, a, b)
end

function _deflection_angle(M, a, b)
    r0 = r₀(abs(b), a, M, sign(b))

    ωs = a / b
    ω₀ = a^2 / M^2
    h = M / r0
    hsc = (1 + ωs) / (1 - ωs)

    k2 = k²(h, hsc)
    ϕ = ψ(h, hsc)

    np = n₊(h, hsc, ω₀)
    nn = n₋(h, hsc, ω₀)

    t1 = Ω₊(h, hsc, ω₀, ωs) * (Π(np, k2) - Π(np, ϕ, k2))
    t2 = Ω₋(h, hsc, ω₀, ωs) * (Π(nn, k2) - Π(nn, ϕ, k2))

    -π + (4 / (1 - ωs)) * √(r₀_Q(h, hsc)) * (t1 + t2)
end

export deflection_angle

end # module IyerHansen2009
