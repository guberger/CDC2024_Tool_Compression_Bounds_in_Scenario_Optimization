module PlotDomain

using LaTeXStrings
using Plots

include("bounds.jl")
colors = palette(:default)
bounds_with_info = [
      (floyd_full_risk, "Floyd", colors[1]),
      (campi_full_risk, "Campi", colors[2]),
      (garatti_wait_risk, "Garatti", colors[3]),
      (berger_full_risk, "us", colors[4]),
]

#-------------------------------------------------------------------------------
β = 0.05
N = 500

plt = plot(xlabel=L"d", ylabel=L"\epsilon",
           title="\$N=$(N)\$, \$q(N,\\epsilon)=$(β)\$", dpi=1000)
ks = collect(0:N)

for (bound, name, color) in bounds_with_info
      risk(k) = bound(k, β, N)
      plot!(plt, ks, risk.(ks), c=color, label=name)
end

savefig(plt, "figures/compare_full.png")

#-------------------------------------------------------------------------------
βs_with_info = [
      (0.05, :solid, "solid"),
      (0.005, :dash, "dash"),
]
title_β = join([
    "$(β)\$ ($(lsn))" for (β, _, lsn) in βs_with_info 
], ", \$", " or \$")
title = join(["\$N=$(N)\$, \$q(N,\\epsilon)=", title_β])

plt = plot(xlabel=L"d", ylabel=L"\epsilon", title=title, dpi=1000)
ks = collect(0:N)
for (β, ls, _) in βs_with_info
    for (bound, name, color) in bounds_with_info
        risk(k) = bound(k, β, N)
        plot!(plt, ks, risk.(ks), c=color, ls=ls, label=false)
    end
end
# labels
for (_, name, color) in bounds_with_info
    plot!(plt, [], [], c=color, ls=:solid, label=name)
end

savefig(plt, "figures/compare_full_beta.png")

#-------------------------------------------------------------------------------
N = 5000

plt = plot(xlabel=L"d", ylabel=L"\epsilon",
           title="\$N=$(N)\$, \$q(N,\\epsilon)=$(β)\$", dpi=1000)
ks = collect(0:N)

for (bound, name, color) in bounds_with_info
      risk(k) = bound(k, β, N)
      plot!(plt, ks, risk.(ks), c=color, label=name)
end

savefig(plt, "figures/compare_full_$(N).png")

end # module