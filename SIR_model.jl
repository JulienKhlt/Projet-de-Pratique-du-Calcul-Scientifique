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

function f(U, β, γ, N)
    dUdt = [- β*U[2]*U[1]/N, β*U[2]*U[1]/N - γ*U[2], γ*U[2]]
    return dUdt
end

function epidemie(S, I, R, β, γ, pas_tps, tps)

    Susceptible = Float64[S]
    Infectious = Float64[I]
    Recovered = Float64[R]

    N = S + I + R
    U = Float64[S, I, R]

    for t in range(pas_tps, tps, length = Int(tps/pas_tps))
        dUdt = f(U, β, γ, N)

        U += dUdt * pas_tps

        append!(Susceptible, U[1])
        append!(Infectious, U[2])
        append!(Recovered, U[3])
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