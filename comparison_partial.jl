module PlotDomain

using LaTeXStrings
using Plots

include("bounds.jl")
colors = palette(:default)
bounds_with_info = [
      (margellos_partial_risk, "Margellos", colors[1]),
      (campi_partial_risk, "Campi", colors[2]),
      (romao_partial_risk, "Romao", colors[3]),
      (berger_partial_risk, "us", colors[4]),
]

#-------------------------------------------------------------------------------
β = 0.05
N = 500
r = 50

plt = plot(xlabel=L"d", ylabel=L"\epsilon",
           title="\$N=$(N)\$, \$r=$(r)\$, \$q(N,\\epsilon)=$(β)\$", dpi=1000)
ks = collect(0:N)

for (bound, name, color) in bounds_with_info
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

plt = plot(xlabel=L"d", ylabel=L"\epsilon", title=title, dpi=1000)
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

savefig(plt, "figures/compare_partial_beta.png")

#-------------------------------------------------------------------------------
rs_with_info = [
      (50, :solid, "solid"),
      (25, :dash, "dash"),
]
title_r = join([
    "$(r)\$ ($(lsn))" for (r, _, lsn) in rs_with_info 
], ", \$", " or \$")
title = join(["\$N=$(N)\$, \$r=", title_r, ", \$q(N,\\epsilon)=$(β)\$"])

plt = plot(xlabel=L"d", ylabel=L"\epsilon", title=title, dpi=1000)
ks = collect(0:N)
for (r, ls, _) in rs_with_info
    for (bound, name, color) in bounds_with_info
        risk(k) = bound(k, β, N, r)
        plot!(plt, ks, risk.(ks), c=color, ls=ls, label=false)
    end
end
# labels
for (_, name, color) in bounds_with_info
    plot!(plt, [], [], c=color, ls=:solid, label=name)
end

savefig(plt, "figures/compare_partial_r.png")

#-------------------------------------------------------------------------------
N = 5000
r = 500

plt = plot(xlabel=L"d", ylabel=L"\epsilon",
           title="\$N=$(N)\$, \$r=$(r)\$, \$q(N,\\epsilon)=$(β)\$", dpi=1000)
ks = collect(0:N)

for (bound, name, color) in bounds_with_info
      risk(k) = bound(k, β, N, r)
      plot!(plt, ks, risk.(ks), c=color, label=name)
end

savefig(plt, "figures/compare_partial_$(N).png")

end # module