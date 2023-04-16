using ITensors, JLD, HDF5

function main(Nx::Int, t::Float64, flux_f::Int, p::String, m::Int, chi0::Int,
        chi::Int)
    nsweeps = 10
    maxdim = [chi]
    cutoff = [1E-8]

    f1 = h5open("data/data0/H_fx_$(Nx)_n_$(m)_p_$(p)_flux_$(flux_f)_chi_$(chi0).h5","r")
    H = read(f1,"H",MPO)
    close(f1)

    f2 = h5open("data/data0/psi_fx_$(Nx)_n_$(m)_p_$(p)_flux_$(flux_f)_chi_$(chi0).h5","r")
    psi0 = read(f2,"psi",MPS)
    close(f2)

    energy,psi = dmrg(H,psi0; nsweeps, maxdim, cutoff)

    f1 = h5open("data/data0/H_fx_$(Nx)_n_$(m)_p_$(p)_flux_$(flux_f)_chi_$(chi).h5","w")
    write(f1,"H",H)
    close(f1)

    f2 = h5open("data/data0/psi_fx_$(Nx)_n_$(m)_p_$(p)_flux_$(flux_f)_chi_$(chi).h5","w")
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
chi0 = parse(Int,ARGS[3])
chi = parse(Int,ARGS[4])

t = 1.0

E_l = Float64[]
S_l = Float64[]

for p in ["F", "A"]
    for m in 2:Nx*4-1
        E, S = main(Nx, t, flux_f, p, m, chi0, chi)
        push!(E_l, E)
        push!(S_l, S)
    end
end

open("data/E_fx_N_$(Nx)_flux_$(flux_f)_chi_$(chi).dat","w") do io
    write(io, E_l)
end
open("data/S_fx_N_$(Nx)_flux_$(flux_f)_chi_$(chi).dat","w") do io
    write(io, S_l)
end

