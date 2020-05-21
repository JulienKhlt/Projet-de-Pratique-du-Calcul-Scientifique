using LinearAlgebra
using SparseArrays

function newton(f, fp, x0, ϵ=1e-10, max_iter=100)
    N = 0
    while N < max_iter && norm(f(x0)) > ϵ
        x0 = x0 - fp(x0) \ f(x0)
        N += 1
    end

    if norm(f(x0)) > ϵ
        error("Newton did not converge.")
    else
        x0
    end
end

# Simple variable function will not work here, vectorize everything.
function euler_implicite(f, fp, x0, δt, T)
    Nδt = trunc(Int, T / δt)
    x = zeros((Nδt, size(x0, 1)))
    x[1, :] = x0

    for i in 2:Nδt
        x[i, :] = newton(y -> y - x[i-1, :] - f(y) * δt, y -> I - fp(y) * δt, x[i-1, :])
    end

    x
end

function euler_explicite(f, x0, δt, T)
    Nδt = trunc(Int, T / δt)
    x = zeros((Nδt, size(x0, 1)))
    x[1, :] = x0

    for i in 2:Nδt
        x[i, :] = x[i - 1, :] + f(x[i - 1, :]) .* δt
    end

    x
end

function laplacian_matrix(Nx, δx)
    Δ = zeros(Nx, Nx)

    for k in 2:(Nx - 1)
        Δ[k, k-1] = 1
        Δ[k, k] = -2
        Δ[k, k+1] = 1
    end

    Δ / (δx^2)
end


# Résout et trace l'équation de la chaleur.
# On vérifie la solution numérique avec la solution analytique.
# On vérifie que la méthode d'euler explicite diverge pour δt trop grand. 
function chaleur()
    Nx = 100
    δx = 1 / Nx
    δt = 0.0001
    T = 1

    Nt = trunc(Int, T / δt)

    Δ = laplacian_matrix(Nx, δx)
    tempX = collect(0:δx:1 - δx)
    x0 = sin.(2 * π * tempX)

    res = euler_explicite(y -> Δ * y, x0, δt, T)
    res2 = euler_implicite(y -> Δ * y, y -> Δ, x0, δt, T)
    sol = zeros((Nt, Nx))

    for t in 1:Nt
        sol[t, :] = sin.(2 * π * tempX) .* exp(- 4 * π^2 * t * δt)
    end

    animation = @animate for i in 1:100
        p = plot()
        plot!(p, res[i, :])
        plot!(p, res2[i, :])
        plot!(p, sol[i, :])
    end

    gif(animation)
end

function SIR(x, β, γ)
    Nx = Int(size(x, 1) / 3)
    y = zeros(3*Nx)
    
    S = x[1:Nx]
    I = x[Nx+1:2*Nx]
    R = x[2*Nx+1:3*Nx]
    
    y[1:Nx] = -β .* S .* I
    y[Nx+1:2*Nx] = β .* S .* I
    y[2*Nx+1:3*Nx] = γ .* I

    y
end

function SIRp(x, β, γ)
    Nx = Int(size(x, 1) / 3)
    y = zeros((3*Nx, 3*Nx))
        
    S = x[1:Nx]
    I = x[Nx+1:2*Nx]
    R = x[2*Nx+1:3*Nx]

    y[1:Nx, 1:Nx] = Diagonal(-β .* I)
    y[Nx+1:2*Nx, Nx+1:2*Nx] = Diagonal(β .* S .- γ)
    y[1:Nx, Nx+1:2*Nx] = Diagonal(-β .* S)
    y[Nx+1:2*Nx, 1:Nx] = Diagonal(β .* S)
    y[2*Nx+1:3*Nx, Nx+1:2*Nx] = Diagonal(γ .* I)

    sparse(y)
end

function SIR_diffusion()
    # Paramètres de calculs numériques
    Nx = 20
    δx = 1 / Nx
    δt = 0.001
    T = 10
    x0 = ones(3 * Nx)
    # Pas d'infecté au début sauf
    x0[Nx+1:2*Nx] .= 0
    x0[Nx+1:Nx+5] .= 1
    
    β = 1
    γ = 1
    
    # On commence par construire la matrice représentant f
    A = zeros((3*Nx, 3*Nx))
    Δ = laplacian_matrix(Nx, δx)

    # On construit les 3 blocs de diffusion
    A[1:Nx, 1:Nx] = Δ
    A[Nx+1:2*Nx, Nx+1:2*Nx] = Δ
    A[2*Nx+1:3*Nx, 2*Nx+1:3*Nx] = Δ
   
    # On transforme la matrice A en trigonale pour la vitesse de calcul
    Atrig = Tridiagonal(A)

    

    res = euler_implicite(y -> A * y + SIR(y, β, γ), y -> A + SIRp(y, β, γ), x0, δt, T)
    res_verif = euler_explicite(y -> A * y + SIR(y, β, γ), x0, δt, T)

    animation = @animate for i in 1:100
        p = plot()
        plot!(p, res[i, 1:Nx])
        plot!(p, res[i, Nx+1:2*Nx])
        plot!(p, res[i, 2*Nx+1:3*Nx])
        plot!(p, res_verif[i, 1:Nx])
        plot!(p, res_verif[i, Nx+1:2*Nx])
        plot!(p, res_verif[i, 2*Nx+1:3*Nx])
    end
    gif(animation)
end
