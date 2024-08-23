using LinearAlgebra
using SpecialFunctions
using Roots

logbinomial(N, d) = logabsbinomial(N, d)[1]
function binlogpdf(N, k, p)::Float64
    if k < 0 || k > N
        return -Inf
    end
    if N == 0 || (p ≤ 0 && k == 0) || (p ≥ 1 && k == N)
        return 0.0
    end
    logbinomial(N, k) + log(p) * k + log(1 - p) * (N - k)
end
binpdf(N, k, p) = exp(binlogpdf(N, k, p))

#-------------------------------------------------------------------------------
# Full consistency

function floyd_strong_risk(d, β, N)::Float64
    if β ≥ 1 || d ≤ 0
        return 0.0
    end
    if β ≤ 0 || d ≥ N
        return 1.0
    end
    return 1 - exp((log(β) - logbinomial(N, d)) / (N - d))
end

function campi_strong_risk(d, β, N)::Float64
    if β ≥ 1 || d ≤ 0
        return 0.0
    end
    if β ≤ 0 || d ≥ N
        return 1.0
    end
    g(i, v) = binlogpdf(N, i, v) - log(β)
    f(v) = 1 - sum(i -> exp(g(i, v)), 0:(d - 1))
    return find_zero(f, (0, 1))
end

function berger_strong_risk(d, β, N)::Float64
    if β ≥ 1 || d ≤ 0
        return 0.0
    end
    if β ≤ 0 || d ≥ N
        return 1.0
    end
    g1(m, v) = log(1 - v) * (N - m) - logbinomial(m, d)
    g2(v) = minimum(m -> g1(m, v), d:N)
    g3(v) = v ≥ 1 ? -Inf : g2(v)
    f(v) = 1 - exp(logbinomial(N, d) + g3(v) - log(β))
    return find_zero(f, (0, 1))
end

#-------------------------------------------------------------------------------
# Wait and Judge

function garatti_posteriori_risk(k, β, N)::Float64
    if β ≥ 1
        return 0.0
    end
    if β ≤ 0 || k ≥ N
        return 1.0
    end
    g1(m, v) = logbinomial(m, k) - logbinomial(N, k) +
               log(1 - v) * (m - N) + log(β) - log(N)
    g2(v) = sum(m -> exp(g1(m, v)), k:(N - 1))
    f(v) = 1 - g2(v)
    return find_zero(f, (0, 1))
end

#-------------------------------------------------------------------------------
# Sample and Discard

function margellos_weak_risk(d, β, N, r)::Float64
    if β ≥ 1 || d ≤ 0
        return 0.0
    end
    if β ≤ 0 || d + r ≥ N
        return 1.0
    end
    g(i, v) = binlogpdf(N - d, i, v) + logbinomial(N, d) - log(β)
    f(v) = 1 - sum(i -> exp(g(i, v)), 0:r)
    return find_zero(f, (0, 1))
end

function campi_weak_risk(d, β, N, r)::Float64
    if β ≥ 1 || d ≤ 0
        return 0.0
    end
    if β ≤ 0 || d + r ≥ N
        return 1.0
    end
    g(i, v) = binlogpdf(N, i, v) + logbinomial(r + d - 1, r) - log(β)
    f(v) = 1 - sum(i -> exp(g(i, v)), 0:(r + d - 1))
    return find_zero(f, (0, 1))
end

function romao_weak_risk(d, β, N, r)::Float64
    if β ≥ 1 || d ≤ 0
        return 0.0
    end
    if β ≤ 0 || d + r ≥ N
        return 1.0
    end
    g(i, v) = binlogpdf(N, i, v) - log(β)
    f(v) = 1 - sum(i -> exp(g(i, v)), 0:(r + d - 1))
    return find_zero(f, (0, 1))
end

function berger_weak_risk(d, β, N, r)::Float64
    if β ≥ 1 || d ≤ 0
        return 0.0
    end
    if β ≤ 0 || d + r ≥ N
        return 1.0
    end
    g1(i, v) = binlogpdf(N - d, i, v) + logbinomial(N, d) - log(β)
    f1(v) = sum(i -> exp(g1(i, v)), 0:r)
    g2(m, v) = log(1 - v) * (N - r - m) - logbinomial(m, d)
    g3(v) = minimum(m -> g2(m, v), d:(N - r))
    g4(v) = v ≥ 1 ? -Inf : g3(v)
    f2(v) = exp(logbinomial(N, r) + logbinomial(N - r, d) + g4(v) - log(β))
    f(v) = 1 - min(f1(v), f2(v))
    return find_zero(f, (0, 1))
end