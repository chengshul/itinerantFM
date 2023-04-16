using ITensors, JLD, HDF5

function main(Nx::Int, t::Float64, t1::Float64, J::Float64, J1::Float64,
                p::String, m::Int)
    N = 4*Nx
    #J = t^2/U
    #J1 = t1^2/U
    sites = siteinds("tJ", N, conserve_qns=true)

    os = OpSum()
    for i in 1:Nx
        for j in 1:3
            r = 4*(i-1)+j
            os .+= -t, "c†↑", r, "c↑", r+1
            os .+= -t, "c†↑", r+1, "c↑", r
            os .+= -t, "c†↓", r, "c↓", r+1
            os .+= -t, "c†↓", r+1, "c↓", r
            os .+= -J, "Sz", r, "Sz", r+1
            os .+= -J/2, "S+", r, "S-", r+1
            os .+= -J/2, "S-", r, "S+", r+1
        end
    end
    for i in 1:Nx-1
        for j in 1:3
            r = 4*(i-1)+j
            os .+= -t1, "c†↑", r, "c↑", r+5
            os .+= -t1, "c†↑", r+5, "c↑", r
            os .+= -t1, "c†↓", r, "c↓", r+5
            os .+= -t1, "c†↓", r+5, "c↓", r
            os .+= -J1, "Sz", r, "Sz", r+5
            os .+= -J1/2, "S+", r, "S-", r+5
            os .+= -J1/2, "S-", r, "S+", r+5
        end
    end
    for i in 1:Nx-1
        for j in 1:4
            r = 4*(i-1)+j
            os .+= -t, "c†↑", r, "c↑", r+4
            os .+= -t, "c†↑", r+4, "c↑", r
            os .+= -t, "c†↓", r, "c↓", r+4
            os .+= -t, "c†↓", r+4, "c↓", r
            os .+= -J, "Sz", r, "Sz", r+4
            os .+= -J/2, "S+", r, "S-", r+4
            os .+= -J/2, "S-", r, "S+", r+4
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

    f1 = h5open("data/data0/H_$(Nx)_n_$(m)_p_$(p)_t1_$(t1)_chi_1000.h5","w")
    write(f1,"H",H)
    close(f1)

    f2 = h5open("data/data0/psi_$(Nx)_n_$(m)_p_$(p)_t1_$(t1)_chi_1000.h5","w")
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
        E_l[i], S_l[i] = main(Nx, t, t1_l[i], J, J1, p, m)
    end

    open("data/E_N_$(Nx)_n_$(m)_p_$(p)_chi_1000.dat","w") do io
        write(io, E_l)
    end
    open("data/S_N_$(Nx)_n_$(m)_p_$(p)_chi_1000.dat","w") do io
        write(io, S_l)
    end
end
