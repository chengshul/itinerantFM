import numpy as np
import matplotlib.pyplot as plt
import sys

def nn(r,Lx,Ly):
    i = (r // 2) // Ly
    j = (r // 2) % Ly
    if r % 2 == 0:
        if i % 2 == 0 and i != Lx:
            return r+1, (i*Ly+(j-1)%Ly)*2+1, r+2*Ly
        elif i == Lx:
            return r+1, (i*Ly+(j-1)%Ly)*2+1
        else:
            return r+1, (i*Ly+(j-1)%Ly)*2+1, r-2*Ly
    else:
        if i == 0:
            return r-1, (i*Ly+(j+1)%Ly)*2
        elif i % 2 == 0:
            return r-1, (i*Ly+(j+1)%Ly)*2, r-2*Ly
        else:
            return r-1, (i*Ly+(j+1)%Ly)*2, r+2*Ly

def nnn(r,Lx,Ly):
    i = (r // 2) // Ly
    j = (r // 2) % Ly
    if r % 2 == 0:
        if i != Lx:
            return (i*Ly+(j+1)%Ly)*2, (i*Ly+(j-1)%Ly)*2
        else:
            return []
    else:
        if i != 0:
            return (i*Ly+(j+1)%Ly)*2+1, (i*Ly+(j-1)%Ly)*2+1
        else:
            return []

def H_k(t1, kx, ky):
    d1 = np.array([np.sqrt(3),0])
    d2 = np.array([np.sqrt(3)/2,3/2])
    d3 = d2 - d1
    #E = np.zeros((L,L))

    #for ix in range(L):
    #    kx = kx_l[ix]
    #    for iy in range(L):
    #        ky = ky_l[iy]
    k = np.array([kx,ky])
    H = np.array([[-t1*np.cos(k@d1),
                   -(1+np.exp(1j*k@d2)+np.exp(1j*k@d3))],
                  [-(1+np.exp(-1j*k@d2)+np.exp(-1j*k@d3)),
                   -t1*np.cos(k@d1)]])
    E = np.linalg.eigvalsh(H)[0]
    return E

def H0(t,t1,Lx,Ly):
    H = np.zeros((2*(Lx+1)*Ly,2*(Lx+1)*Ly))
    for r in range(2*(Lx+1)*Ly):
        for r1 in nn(r,Lx,Ly):
            H[r1,r] += -t
        for r1 in nnn(r,Lx,Ly):
            H[r1,r] += -t1
    return H

def H(t, t1, Lx, Ly):
    d1 = np.zeros((2*(Lx+1)*Ly, 2*(Lx+1)*Ly), dtype=int)
    d2 = np.zeros(((Lx+1)*Ly*(2*(Lx+1)*Ly-1), 2), dtype=int)
    c = 0
    for r1 in range(2*(Lx+1)*Ly):
        for r2 in range(r1+1,2*(Lx+1)*Ly):
            d1[r1,r2] = c
            d1[r2,r1] = c
            d2[c,:] = [r1,r2]
            c += 1

    Ht = np.zeros((c,c))
    Hs = np.zeros((c,c))
    for i0 in range(c):
        r1, r2 = d2[i0,:]
        for r1p in nn(r1, Lx, Ly):
            if r1p != r2:
                i1 = d1[r1p,r2]
                Ht[i1,i0] += -t*(-1)**int(r1p>r2)
                Hs[i1,i0] += -t
        for r2p in nn(r2, Lx, Ly):
            if r2p != r1:
                i1 = d1[r1,r2p]
                Ht[i1,i0] += -t*(-1)**int(r1>r2p)
                Hs[i1,i0] += -t
        for r1p in nnn(r1, Lx, Ly):
            if r1p != r2:
                i1 = d1[r1p,r2]
                Ht[i1,i0] += -t1*(-1)**int(r1p>r2)
                Hs[i1,i0] += -t1
        for r2p in nnn(r2, Lx, Ly):
            if r2p != r1:
                i1 = d1[r1,r2p]
                Ht[i1,i0] += -t1*(-1)**int(r1>r2p)
                Hs[i1,i0] += -t1
    return Ht, Hs

t = 1
t1_l = np.linspace(-1,0,51)
Lx = int(sys.argv[1])
Ly = int(sys.argv[2])
check = 1

E_all = []
E0_all = []
for t1 in t1_l:
    print(t1)
    Ht, Hs = H(t,t1,Lx,Ly)
    Et, Ut = np.linalg.eigh(Ht)
    Es, Us = np.linalg.eigh(Hs)
    if check == 1:
        H_0 = H0(t,t1,Lx,Ly)
        E0 = []
        e0, u0 = np.linalg.eigh(H_0)
        for i in range(2*(Lx+1)*Ly):
            for j in range(i+1,2*(Lx+1)*Ly):
                E0.append(e0[i]+e0[j])
        E0_all.append(np.sort(E0))
    E_all.append([Et,Es])
E_all = np.array(E_all)
#tc = t1_l[20+np.argmin(np.abs(E_all[20:,0,0]-E_all[20:,1,0]))]
#print("tc=",tc)

plt.plot(t1_l,E_all[:,0,0],label="$S=1$")
plt.plot(t1_l,E_all[:,1,0],label="$S=0$")
plt.legend()
plt.xlabel("$t_1$")
plt.ylabel("$E_0$")
plt.title(f"$N={Lx}\\times{Ly}$")
plt.savefig(f"plot/E0_honeycomb_{Lx}_{Ly}_cyl.pdf",format='pdf')

'''
L = 101
kx_l = np.linspace(-4*np.pi/3/np.sqrt(3),4*np.pi/3/np.sqrt(3),L)
ky_l = np.linspace(-2*np.pi/3,2*np.pi/3,L)

plt.figure()
t1_l = np.round(np.linspace(-1,0,11),4)
Ek_l = []
for t1 in t1_l:
    E = []
    for kx in kx_l:
        E.append(H_k(t1, kx, 0))
    Ek_l.append(E)

for i in range(len(t1_l)):
    plt.plot(kx_l,Ek_l[i],label=f"$t_1={t1_l[i]}$")
plt.legend()
'''

'''E = H_k(101,-0.8)
plt.imshow(E.T)
plt.colorbar()

plt.figure()
plt.plot(E[:,50])'''
plt.show()
