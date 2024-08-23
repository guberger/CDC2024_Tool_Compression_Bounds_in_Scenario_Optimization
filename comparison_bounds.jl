module PlotDomain

using LaTeXStrings
using Plots

include("bounds.jl")
colors = palette(:default)
bounds_with_info = [
    ((k, β, N, _) -> floyd_strong_risk(k, β, N), "Theorem 1", colors[1]),
    ((k, β, N, _) -> campi_strong_risk(k, β, N), "Theorem 2", colors[2]),
    ((k, β, N, _) -> garatti_posteriori_risk(k, β, N), "Theorem 3", colors[3]),
    ((k, β, N, _) -> berger_strong_risk(k, β, N), "Theorem 4", colors[4]),
    (campi_weak_risk, "Theorem 5", colors[5]),
    (margellos_weak_risk, "Theorem 6", colors[6]),
    (romao_weak_risk, "Theorem 7", colors[7]),
    (berger_weak_risk, "Theorem 8", :black),
]
bounds_strong_with_info = bounds_with_info[1:4]
bounds_weak_with_info = bounds_with_info[5:8]
β = 0.05
N = 500
r = 50

#-------------------------------------------------------------------------------
plt = plot(xlabel=L"d", ylabel=L"\epsilon", legend=:bottomright,
           title="\$N=$(N)\$, \$q(N,\\epsilon)=$(β)\$", dpi=1000)
ks = collect(0:N)
R = fill(NaN, length(ks), length(bounds_strong_with_info) + 1)
R[:, 1] = ks

for (i, (bound, name, color)) in enumerate(bounds_strong_with_info)
    risk(k) = bound(k, β, N, Inf)
    # R[:, i + 1] = risk.(ks)
    # plot!(plt, ks, R[:, i + 1], c=color, label=name)
end

savefig(plt, "figures/compare_strong.png")
file = open("data_strong.txt", "w")
for row in eachrow(R)
    for val in row
        print(file, val, " ")
    end
    println(file, "")
end
close(file)

#-------------------------------------------------------------------------------
plt = plot(xlabel=L"d", ylabel=L"\epsilon", legend=:bottomright,
           title="\$N=$(N)\$, \$r=$(r)\$, \$q(N,\\epsilon)=$(β)\$", dpi=1000)
ks = collect(0:N)
R = fill(NaN, length(ks), length(bounds_weak_with_info) + 1)
R[:, 1] = ks

for (i, (bound, name, color)) in enumerate(bounds_weak_with_info)
    risk(k) = bound(k, β, N, r)
    # R[:, i + 1] = risk.(ks)
    # plot!(plt, ks, R[:, i + 1], c=color, label=name)
end

savefig(plt, "figures/compare_weak.png")
file = open("data_weak.txt", "w")
for row in eachrow(R)
    for val in row
        print(file, val, " ")
    end
    println(file, "")
end
close(file)

#-------------------------------------------------------------------------------
β = 0.005

plt = plot(xlabel=L"d", ylabel=L"\epsilon", legend=:bottomright,
           title="\$N=$(N)\$, \$r=$(r)\$, \$q(N,\\epsilon)=$(β)\$", dpi=1000)
ks = collect(0:N)
R = fill(NaN, length(ks), length(bounds_with_info) + 1)
R[:, 1] = ks

for (i, (bound, name, color)) in enumerate(bounds_with_info)
    risk(k) = bound(k, β, N, r)
    # R[:, i + 1] = risk.(ks)
    # plot!(plt, ks, R[:, i + 1], c=color, label=name)
end

savefig(plt, "figures/compare_beta.png")
file = open("data_beta.txt", "w")
for row in eachrow(R)
    for val in row
        print(file, val, " ")
    end
    println(file, "")
end
close(file)

β = 0.05

#-------------------------------------------------------------------------------
N = 5000
r = 500

plt = plot(xlabel=L"d", ylabel=L"\epsilon", legend=:bottomright,
           title="\$N=$(N)\$, \$r=$(r)\$, \$q(N,\\epsilon)=$(β)\$", dpi=1000)
ks = collect(0:10:N)
R = fill(NaN, length(ks), length(bounds_with_info) + 1)
R[:, 1] = ks

for (i, (bound, name, color)) in enumerate(bounds_with_info)
    risk(k) = bound(k, β, N, r)
    # R[:, i + 1] = risk.(ks)
    # plot!(plt, ks, R[:, i + 1], c=color, label=name)
end

savefig(plt, "figures/compare_N.png")
file = open("data_N.txt", "w")
for row in eachrow(R)
    for val in row
        print(file, val, " ")
    end
    println(file, "")
end
close(file)

N = 500
r = 50

#-------------------------------------------------------------------------------
r = 25

plt = plot(xlabel=L"d", ylabel=L"\epsilon", legend=:bottomright,
           title="\$N=$(N)\$, \$r=$(r)\$, \$q(N,\\epsilon)=$(β)\$", dpi=1000)
ks = collect(0:N)
R = fill(NaN, length(ks), length(bounds_weak_with_info) + 1)
R[:, 1] = ks

for (i, (bound, name, color)) in enumerate(bounds_weak_with_info)
    risk(k) = bound(k, β, N, r)
    R[:, i + 1] = risk.(ks)
    plot!(plt, ks, R[:, i + 1], c=color, label=name)
end

savefig(plt, "figures/compare_r.png")
file = open("data_r.txt", "w")
for row in eachrow(R)
    for val in row
        print(file, val, " ")
    end
    println(file, "")
end
close(file)

end # module