using PyPlot

Nombre_maille = 100
x = range(0, 1, length=Nombre_maille)
dx = 1/Nombre_maille
D = 1.0

U = [zeros(Nombre_maille) for _ in 1:3]
laplacien = [zeros(Nombre_maille) for _ in 1:3]

for i = 1:Nombre_maille
    U[1][i] = 1
    U[2][i] = 0
    U[3][i] = 0
end

U[2][50] = 1

function f(U, β, γ, N)
    dUdt = [- β*U[2].*U[1]/N, β*U[2].*U[1]/N - γ*U[2], γ*U[2]]
    return dUdt
end

# stable si Ddt/dx^2 =< 1/2

pas_tps = 3e-2
tps = 10
β = 1
γ = 0.25

function diffusion(U, β, γ, D, dx, pas_tps, tps)
    Epidemie = [U]
    N = Nombre_maille + 1
    for t in range(pas_tps, tps, length = trunc(Int, (tps/pas_tps)))
        for k in 2:Nombre_maille-1
            for i = 1:length(U)
            laplacien[i][k] = (U[i][k+1]-2*U[i][k]+U[i][k-1])/(dx^2)
            end
        end
        U .+= pas_tps*(D*laplacien + f(U, β, γ, N)) 
        append!(Epidemie, [U])
        println(U)
        
    end
    println(Epidemie)
end

diffusion(U, β, γ, D, dx, pas_tps, tps)

# f(U, β, γ, Nombre_maille+1)