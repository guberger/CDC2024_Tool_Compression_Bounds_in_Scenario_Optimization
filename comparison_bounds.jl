module PlotDomain

using LaTeXStrings
using Plots

include("bounds.jl")
colors = palette(:default)
bounds_with_info = [
    ((k, β, N, _) -> floyd_full_risk(k, β, N), "Theorem 1", colors[1]),
    ((k, β, N, _) -> campi_full_risk(k, β, N), "Theorem 2", colors[2]),
    ((k, β, N, _) -> garatti_wait_risk(k, β, N), "Theorem 3", colors[3]),
    ((k, β, N, _) -> berger_full_risk(k, β, N), "Theorem 4", colors[4]),
    (campi_partial_risk, "Theorem 5", colors[5]),
    (margellos_partial_risk, "Theorem 6", colors[6]),
    (romao_partial_risk, "Theorem 7", colors[7]),
    (berger_partial_risk, "Theorem 8", :black),
]
bounds_full_with_info = bounds_with_info[1:4]
bounds_partial_with_info = bounds_with_info[5:8]
β = 0.05
N = 500
r = 50

#-------------------------------------------------------------------------------
plt = plot(xlabel=L"d", ylabel=L"\epsilon", legend=:bottomright,
           title="\$N=$(N)\$, \$q(N,\\epsilon)=$(β)\$", dpi=1000)
ks = collect(0:N)

for (bound, name, color) in bounds_full_with_info
    risk(k) = bound(k, β, N, Inf)
    plot!(plt, ks, risk.(ks), c=color, label=name)
end

savefig(plt, "figures/compare_full.png")

#-------------------------------------------------------------------------------
plt = plot(xlabel=L"d", ylabel=L"\epsilon", legend=:bottomright,
           title="\$N=$(N)\$, \$r=$(r)\$, \$q(N,\\epsilon)=$(β)\$", dpi=1000)
ks = collect(0:N)

for (bound, name, color) in bounds_partial_with_info
    risk(k) = bound(k, β, N, r)
    plot!(plt, ks, risk.(ks), c=color, label=name)
end

savefig(plt, "figures/compare_partial.png")

#-------------------------------------------------------------------------------
βs_with_info = [
    (0.05, :solid, "solid"),
    (0.005, :dash, "dash"),
]
title_β = join([
    "$(β)\$ ($(lsn))" for (β, _, lsn) in βs_with_info 
], ", \$", " or \$")
title = join(["\$N=$(N)\$, \$r=$(r)\$, \$q(N,\\epsilon)=", title_β])

plt = plot(xlabel=L"d", ylabel=L"\epsilon", legend=:bottomright,
           title=title, dpi=1000)
ks = collect(0:N)
for (β, ls, _) in βs_with_info
    for (bound, name, color) in bounds_with_info
        risk(k) = bound(k, β, N, r)
        plot!(plt, ks, risk.(ks), c=color, ls=ls, label=false)
    end
end
# labels
for (_, name, color) in bounds_with_info
    plot!(plt, [], [], c=color, ls=:solid, label=name)
end

savefig(plt, "figures/compare_beta.png")

#-------------------------------------------------------------------------------
Ntest = 5000
rtest = 500
kstest = collect(0:Ntest)

plt = plot(xlabel=L"d", ylabel=L"\epsilon", legend=:bottomright,
           title="\$N=$(Ntest)\$, \$r=$(rtest)\$, \$q(N,\\epsilon)=$(β)\$",
           dpi=1000)
ks = collect(0:N)

for (bound, name, color) in bounds_with_info
      risk(k) = bound(k, β, Ntest, rtest)
      plot!(plt, kstest, risk.(kstest), c=color, label=name)
end

savefig(plt, "figures/compare_N.png")

#-------------------------------------------------------------------------------
rs_with_info = [
    (50, :solid, "solid"),
    (25, :dash, "dash"),
]
title_r = join([
    "$(r)\$ ($(lsn))" for (r, _, lsn) in rs_with_info 
], ", \$", " or \$")
title = join(["\$N=$(N)\$, \$r=", title_r, ", \$q(N,\\epsilon)=$(β)\$"])

plt = plot(xlabel=L"d", ylabel=L"\epsilon", legend=:bottomright,
           title=title, dpi=1000)
ks = collect(0:N)
for (r, ls, _) in rs_with_info
    for (bound, name, color) in bounds_partial_with_info
        risk(k) = bound(k, β, N, r)
        plot!(plt, ks, risk.(ks), c=color, ls=ls, label=false)
    end
end
# labels
for (_, name, color) in bounds_partial_with_info
    plot!(plt, [], [], c=color, ls=:solid, label=name)
end

savefig(plt, "figures/compare_r.png")

end # module