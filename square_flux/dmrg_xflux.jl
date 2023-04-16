using ITensors, JLD, HDF5

function main(Nx::Int, t::Float64, flux_f::Int, p::String, m::Int)
    N = 4*Nx
    #J = t^2/U
    #J1 = t1^2/U
    sites = siteinds("tJ", N, conserve_qns=true)
    
    flux = exp(2.0im*π/flux_f)

    os = OpSum()
    for i in 1:Nx
        for j in 1:3
            r = 4*(i-1)+j
            os .+= -t*flux^i, "c†↑", r, "c↑", r+1
            os .+= -t*conj(flux^i), "c†↑", r+1, "c↑", r
            os .+= -t*flux^i, "c†↓", r, "c↓", r+1
            os .+= -t*conj(flux^i), "c†↓", r+1, "c↓", r
        end
    end
    for i in 1:Nx-1
        for j in 1:4
            r = 4*(i-1)+j
            os .+= -t, "c†↑", r, "c↑", r+4
            os .+= -t, "c†↑", r+4, "c↑", r
            os .+= -t, "c†↓", r, "c↓", r+4
            os .+= -t, "c†↓", r+4, "c↓", r
        end
    end
    H = MPO(os,sites)

    state = ["0" for n=1:N]

    for n in 1:m
        if n % 2 == 1
            state[n] = "↑"
        else
            state[n] = (p == "F" ? "↑" : "↓")
        end
    end
    
    psi0 = randomMPS(sites,state,10)
    nsweeps = 35
    maxdim = [20,60,100,100,200,200,400,400,600,600,
              800,800,800,800,800,1000,1000,1000,1000,1000]
    #          1200,1200,1200,1200,1200,1400]
    cutoff = [1E-8]

    energy,psi = dmrg(H,psi0; nsweeps, maxdim, cutoff)

    f1 = h5open("data/data0/H_fx_$(Nx)_n_$(m)_p_$(p)_flux_$(flux_f)_chi_1000.h5","w")
    write(f1,"H",H)
    close(f1)

    f2 = h5open("data/data0/psi_fx_$(Nx)_n_$(m)_p_$(p)_flux_$(flux_f)_chi_1000.h5","w")
    write(f2,"psi",psi)
    close(f2)

    zz = correlation_matrix(psi,"Sz","Sz")
    pm = correlation_matrix(psi,"S+","S-")
    mp = correlation_matrix(psi,"S-","S+")
    Stot = real(sum(zz+pm/2+mp/2))

    return energy, Stot
end

Nx = parse(Int,ARGS[1])
flux_f = parse(Int,ARGS[2])

t = 1.0

E_l = Float64[]
S_l = Float64[]

for p in ["F", "A"]
    for m in 2:Nx*4-1
        E, S = main(Nx, t, flux_f, p, m)
        push!(E_l, E)
        push!(S_l, S)
    end
end

open("data/E_fx_N_$(Nx)_flux_$(flux_f)_chi_1000.dat","w") do io
    write(io, E_l)
end
open("data/S_fx_N_$(Nx)_flux_$(flux_f)_chi_1000.dat","w") do io
    write(io, S_l)
end

