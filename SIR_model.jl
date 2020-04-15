using PyPlot

# Conditions initiales:

S = 99
I = 1
R = 0

# Paramètre de l'épidémie:

β = 1
γ = 0.25
pas_tps = 0.1
tps = 40
T = range(0, tps, length = Int(tps/pas_tps))

function epidemie(S, I, R, β, γ, pas_tps, tps)

    Susceptible = Float64[S]
    Infectious = Float64[I]
    Recovered = Float64[R]

    N = S + I + R

    for t in range(pas_tps, tps, length = Int(tps/pas_tps))
        dSdt = - β*I*S/N
        dIdt = β*I*S/N - γ*I
        dRdt = γ*I

        S += dSdt*pas_tps
        I += dIdt*pas_tps
        R += dRdt*pas_tps

        append!(Susceptible, S)
        append!(Infectious, I)
        append!(Recovered, R)
    end
    T = range(0, tps, length = Int(tps/pas_tps)+1)
    plot(T, Susceptible, color="#4a4c4d", linewidth=2.0, linestyle="-", label="Susceptible")
    plot(T, Infectious, color="#c13607", linewidth=2.0, linestyle="-", label="Infectious")
    plot(T, Recovered, color="#3c9f66", linewidth=2.0, linestyle="-", label="Recovered")
    title(L"Model \ SIR")
    xlabel(L"Temps")
    ylabel(L"Nombre \ de \ personne")
    legend()
    grid("on")
    show()
end

epidemie(S, I, R, β, γ, pas_tps, tps)