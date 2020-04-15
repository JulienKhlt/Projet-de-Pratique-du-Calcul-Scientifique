using PyPlot

# Conditions initiales:

S = 99
I = 1
R = 0

# Paramètre de l'épidémie:

β = 1
γ = 1
pas_tps = 0.1
tps = 100
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
    plot(T, Susceptible, "x")
    plot(T, Infectious, "x")
    plot(T, Recovered, "x")
    show()
end

epidemie(S, I, R, β, γ, pas_tps, tps)