using ITensors, JLD, HDF5, Random

function main(Nx::Int, Ny::Int, t::Float64, t1::Float64, p::String, m::Int)
    N = 2*(Ny+1)*Nx+2*Ny
    #J = t^2/U
    #J1 = t1^2/U
    sites = siteinds("tJ", N, conserve_qns=true)

    os = OpSum()
    for i in 0:Nx-1
        for j in 0:2*Ny
            r = 2*(Ny+1)*i+j+1
            r1 = 2*(Ny+1)*i+j+2

            os .+= -t, "c†↑", r, "c↑", r1
            os .+= -t, "c†↑", r1, "c↑", r
            os .+= -t, "c†↓", r, "c↓", r1
            os .+= -t, "c†↓", r1, "c↓", r
        end
    end
   
    i = Nx
    for j in 0:2*Ny-2
        r = 2*(Ny+1)*i+j+1
        r1 = 2*(Ny+1)*i+j+2

        os .+= -t, "c†↑", r, "c↑", r1
        os .+= -t, "c†↑", r1, "c↑", r
        os .+= -t, "c†↓", r, "c↓", r1
        os .+= -t, "c†↓", r1, "c↓", r
    end

    for i in 0:Nx-2
        for j in 0:Ny÷2-1
            r = 2*(Ny+1)*i+4*j+1
            r1 = r+2*Ny+3
            r2 = 2*(Ny+1)*i+4*j+4
            r3 = r2+2*Ny+1

            os .+= -t, "c†↑", r, "c↑", r1
            os .+= -t, "c†↑", r1, "c↑", r
            os .+= -t, "c†↓", r, "c↓", r1
            os .+= -t, "c†↓", r1, "c↓", r
            os .+= -t, "c†↑", r2, "c↑", r3
            os .+= -t, "c†↑", r3, "c↑", r2
            os .+= -t, "c†↓", r2, "c↓", r3
            os .+= -t, "c†↓", r3, "c↓", r2
        end
        r = 2*(Ny+1)*i+2*Ny+2
        r1 = r+2*Ny+1

        os .+= -t, "c†↑", r, "c↑", r1
        os .+= -t, "c†↑", r1, "c↑", r
        os .+= -t, "c†↓", r, "c↓", r1
        os .+= -t, "c†↓", r1, "c↓", r
    end
   
    i = Nx-1
    for j in 0:Ny÷2-1
        r = 2*(Ny+1)*i+4*j+1
        r1 = r+2*Ny+2
        r2 = 2*(Ny+1)*i+4*j+4
        r3 = r2+2*Ny

        os .+= -t, "c†↑", r, "c↑", r1
        os .+= -t, "c†↑", r1, "c↑", r
        os .+= -t, "c†↓", r, "c↓", r1
        os .+= -t, "c†↓", r1, "c↓", r
        os .+= -t, "c†↑", r2, "c↑", r3
        os .+= -t, "c†↑", r3, "c↑", r2
        os .+= -t, "c†↓", r2, "c↓", r3
        os .+= -t, "c†↓", r3, "c↓", r2
    end
    r = 2*(Ny+1)*i+2*Ny+2
    r1 = r+2*Ny

    os .+= -t, "c†↑", r, "c↑", r1
    os .+= -t, "c†↑", r1, "c↑", r
    os .+= -t, "c†↓", r, "c↓", r1
    os .+= -t, "c†↓", r1, "c↓", r
   
    for i in 0:Nx-2
        for j in 1:2*Ny
            r = 2*(Ny+1)*i+j+1
            r2 = 2*(Ny+1)*i+j+1+2*(Ny+1)
            
            os .+= -t1, "c†↑", r, "c↑", r2
            os .+= -t1, "c†↑", r2, "c↑", r
            os .+= -t1, "c†↓", r, "c↓", r2
            os .+= -t1, "c†↓", r2, "c↓", r
        end
    end

    i = Nx-1
    for j in 1:2*Ny
        r = 2*(Ny+1)*i+j+1
        r2 = 2*(Ny+1)*i+j+1+2*Ny+1
        
        os .+= -t1, "c†↑", r, "c↑", r2
        os .+= -t1, "c†↑", r2, "c↑", r
        os .+= -t1, "c†↓", r, "c↓", r2
        os .+= -t1, "c†↓", r2, "c↓", r
    end

    H = MPO(os,sites)

    state = ["0" for n=1:N]

    n_r = randperm(N)
    for n in 1:m
        if n % 2 == 1
            state[n_r[n]] = "↑"
        else
            state[n_r[n]] = (p == "F" ? "↑" : "↓")
        end
    end
    
    psi0 = randomMPS(sites,state,10)
    nsweeps = 35
    maxdim = [20,60,100,100,200,200,400,400,600,600,
              800,800,800,800,800,1000,1000,1000,1000,1000]
    #          1200,1200,1200,1200,1200,1400]
    cutoff = [1E-8]
    noise = [1e-6,1e-6,1e-8,1e-8,0.0]

    energy,psi = dmrg(H,psi0; nsweeps, maxdim, cutoff, noise)

    f1 = h5open("data/data0/H_2l_obc_$(Nx)_n_$(m)_p_$(p)_t1_$(t1)_chi_1000.h5","w")
    write(f1,"H",H)
    close(f1)

    f2 = h5open("data/data0/psi_2l_obc_$(Nx)_n_$(m)_p_$(p)_t1_$(t1)_chi_1000.h5","w")
    write(f2,"psi",psi)
    close(f2)

    zz = correlation_matrix(psi,"Sz","Sz")
    pm = correlation_matrix(psi,"S+","S-")
    mp = correlation_matrix(psi,"S-","S+")
    Stot = real(sum(zz+pm/2+mp/2))

    return energy, Stot
end

Nx = parse(Int,ARGS[1])
m = parse(Int,ARGS[2])
Ny = 2

t = 1.0
t1_l = [-3.0,-2.0,-1.5,-1,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,
        0.0,0.1,0.2,0.3,0.4,0.5,0.8,1.0,1.5,2.0,3.0]

for p in ["F", "A"]
    E_l = Float64[]
    S_l = Float64[]

    for t1 in t1_l
        E, S = main(Nx, Ny, t, t1, p, m)
        push!(E_l, E)
        push!(S_l, S)
    end
    
    open("data/E_2l_obc_N_$(Nx)_n_$(m)_p_$(p)_chi_1000.dat","w") do io
        write(io, E_l)
    end
    open("data/S_2l_obc_N_$(Nx)_n_$(m)_p_$(p)_chi_1000.dat","w") do io
        write(io, S_l)
    end 
end

