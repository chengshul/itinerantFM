using ITensors, JLD, HDF5

function main(Nx::Int, t::Float64, t1::Float64, J::Float64, J1::Float64,
                p::String, m::Int, chi0::Int, chi::Int)
    nsweeps = 10
    maxdim = [chi]
    cutoff = [1E-8]

    f1 = h5open("data/data0/H_2l_$(Nx)_n_$(m)_p_$(p)_t1_$(t1)_chi_$(chi0).h5","r")
    H = read(f1,"H",MPO)
    close(f1)

    f2 = h5open("data/data0/psi_2l_$(Nx)_n_$(m)_p_$(p)_t1_$(t1)_chi_$(chi0).h5","r")
    psi0 = read(f2,"psi",MPS)
    close(f2)

    energy,psi = dmrg(H,psi0; nsweeps, maxdim, cutoff)

    f1 = h5open("data/data0/H_2l_$(Nx)_n_$(m)_p_$(p)_t1_$(t1)_chi_$(chi).h5","w")
    write(f1,"H",H)
    close(f1)

    f2 = h5open("data/data0/psi_2l_$(Nx)_n_$(m)_p_$(p)_t1_$(t1)_chi_$(chi).h5","w")
    write(f2,"psi",psi)
    close(f2)

    zz = correlation_matrix(psi,"Sz","Sz")
    pm = correlation_matrix(psi,"S+","S-")
    mp = correlation_matrix(psi,"S-","S+")
    Stot = sum(zz+pm/2+mp/2)

    return energy, Stot
end

Nx = parse(Int,ARGS[1])
m = parse(Int,ARGS[2])
chi0 = parse(Int,ARGS[3])
chi = parse(Int,ARGS[4])

t = 1.0
#t1_l = range(-3.0,3.0,31)
t1_l = [-3.0,-2.0,-1.5,-1.4,-1.3,-1.2,-1.1,-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,
        -0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.8,1.0,1.5,2.0,3.0]
J = 0.0
J1 = 0.0
E_l = zeros(length(t1_l))
S_l = zeros(length(t1_l))
for p in ["F", "A"]
    for i in 1:length(t1_l)
        E_l[i], S_l[i] = main(Nx, t, t1_l[i], J, J1, p, m, chi0, chi)
    end

    open("data/E_2l_N_$(Nx)_n_$(m)_p_$(p)_chi_$(chi).dat","w") do io
        write(io, E_l)
    end
    open("data/S_2l_N_$(Nx)_n_$(m)_p_$(p)_chi_$(chi).dat","w") do io
        write(io, S_l)
    end
end
