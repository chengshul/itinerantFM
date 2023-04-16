import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse.linalg import eigsh
import sys

def nn(r,Lx,Ly,bc):
    i = r // Ly
    j = r % Ly
    if bc == 'p' or (i < Lx-1 and j < Ly - 1 and i > 0 and j > 0):
        return (i-1)%Lx*Ly+j,(i+1)%Lx*Ly+j,i*Ly+(j-1)%Ly,i*Ly+(j+1)%Ly
    elif i == Lx-1 and j < Ly - 1 and j > 0:
        return (i-1)*Ly+j,i*Ly+(j-1),i*Ly+(j+1)
    elif i == 0 and j < Ly - 1 and j > 0:
        return (i+1)*Ly+j,i*Ly+(j-1),i*Ly+(j+1)
    elif j == Ly-1 and i < Lx - 1 and i > 0:
        return (i-1)*Ly+j,(i+1)*Ly+j,i*Ly+(j-1)
    elif j == 0 and i < Lx - 1 and i > 0:
        return (i-1)*Ly+j,(i+1)*Ly+j,i*Ly+(j+1)
    elif j == 0 and i == 0:
        return (i+1)*Ly+j,i*Ly+(j+1)
    elif j == 0 and i == Lx - 1:
        return (i-1)*Ly+j,i*Ly+(j+1)
    elif j == Ly-1 and i == 0:
        return (i+1)*Ly+j,i*Ly+(j-1)
    elif j == Ly-1 and i == Lx - 1:
        return (i-1)*Ly+j,i*Ly+(j-1)
    else:
        return []

def nnn(r,Lx,Ly,bc):
    i = r // Ly
    j = r % Ly
    if bc == 'p' or (i > 0 and i < Lx-1 and j > 0 and j < Ly-1):
        return (i-1)%Lx*Ly+(j-1)%Ly,(i+1)%Lx*Ly+(j+1)%Ly
    elif i < Lx-1 and j < Ly-1:
        return [(i+1)*Ly+(j+1)]
    elif i > 0 and j > 0:
        return [(i-1)*Ly+(j-1)]
    else:
        return []

def H0(t,t1,Lx,Ly,bc):
    H = np.zeros((Lx*Ly,Lx*Ly))
    for i in range(Lx):
        for j in range(Ly):
            r = i*Ly+j
            for r1 in nn(r,Lx,Ly,bc):
                H[r1,r] += -t
            for r1 in nnn(r,Lx,Ly,bc):
                H[r1,r] += -t1
    return H

def single_to_many(p1,p2,sgn):
    d = len(p1)
    D = d*(d-1)//2
    P = np.zeros(D)

    dic = np.zeros((d,d),dtype=int)
    c = 0
    for i in range(d):
        for j in range(i+1,d):
            dic[i,j] = c
            #dic[j,i] = c
            c += 1

    for i in range(d):
        for j in range(d):
            if i < j:
                P[dic[i,j]] += p1[i]*p2[j]
            elif i > j:
                P[dic[j,i]] += p1[i]*p2[j]*sgn
    P = P / np.sqrt(np.conj(P)@P)
    return P

def H(t,t1,J,J1,Lx,Ly,bc):
    #J = t**2/U
    #J1 = t1**2/U

    d1 = np.zeros((Lx*Ly,Lx*Ly),dtype=int)
    d2 = np.zeros((Lx*Ly*(Lx*Ly-1)//2,2),dtype=int)
    c = 0
    for i1 in range(Lx):
        for j1 in range(Ly):
            r1 = i1*Ly+j1
            for i2 in range(Lx):
                for j2 in range(Ly):
                    r2 = i2*Ly+j2
                    if r2 > r1:
                        d1[r1,r2] = c
                        d1[r2,r1] = c
                        d2[c,:] = [r1,r2]
                        c += 1

    Ht = np.zeros((c,c))
    Hs = np.zeros((c,c))
    for i0 in range(c):
        r1,r2 = d2[i0,:]
        
        for r1p in nn(r1,Lx,Ly,bc):
            if r1p != r2:
                i1 = d1[r1p,r2]
                Ht[i1,i0] += -t*(-1)**int(r1p>r2)
                Hs[i1,i0] += -t
            else:
                Ht[i0,i0] += J*3/4
                Hs[i0,i0] += -J/4
        for r2p in nn(r2,Lx,Ly,bc):
            if r2p != r1:
                i1 = d1[r1,r2p]
                Ht[i1,i0] += -t*(-1)**int(r1>r2p)
                Hs[i1,i0] += -t

        for r1p in nnn(r1,Lx,Ly,bc):
            if r1p != r2:
                i1 = d1[r1p,r2]
                Ht[i1,i0] += -t1*(-1)**int(r1p>r2)
                Hs[i1,i0] += -t1
            else:
                Ht[i0,i0] += J1*3/4
                Hs[i0,i0] += -J1/4
        for r2p in nnn(r2,Lx,Ly,bc):
            if r2p != r1:
                i1 = d1[r1,r2p]
                Ht[i1,i0] += -t1*(-1)**int(r1>r2p)
                Hs[i1,i0] += -t1
    return Ht, Hs

t = 1
#t1 = 1.2
#t1_l = [-1.5,1.5] #np.linspace(-3,3,61)
#t1_l = np.linspace(-3,3,121)
t1_l = np.linspace(-1.5,0,121)
#U = 10
J = 0
J1 = 0
Lx = int(sys.argv[1])
Ly = int(sys.argv[2])
check = 1
bc = sys.argv[3]

E_all = []
E0_all = []
E_gutz_all = []
psipsi_all = []
for t1 in t1_l:
    print(t1)
    Ht, Hs = H(t,t1,J,J1,Lx,Ly,bc)
    Et, Ut = np.linalg.eigh(Ht)
    Es, Us = np.linalg.eigh(Hs)
    #Et = eigsh(Ht,return_eigenvectors=False)
    #Es = eigsh(Hs,return_eigenvectors=False)
    if check == 1:
        H_0 = H0(t,t1,Lx,Ly,bc)
        E0 = []
        e0, u0 = np.linalg.eigh(H_0)
        for i in range(Lx*Ly):
            for j in range(i+1,Lx*Ly):
                E0.append(e0[i]+e0[j])
        E0_all.append(np.sort(E0))
        P_gutz = single_to_many(u0[:,0],u0[:,0],1)
        E_gutz_all.append(np.conj(P_gutz) @ Hs @ P_gutz)
        psipsi_all.append(np.abs(np.conj(P_gutz) @ Us[:,0]))
    #print(Et[:5])
    #print(Es[:5])
    E_all.append([Et,Es])
E_all = np.array(E_all)
tc = t1_l[np.argmin(np.abs(E_all[:,0,0]-E_all[:,1,0]))]
print("tc=",tc)
#tc_gutz = t1_l[200+np.argmin(np.abs(E_gutz_all[200:300]-E_all[200:300,0,0]))]
#print("tc_gutz=",tc_gutz)
plt.plot(t1_l,E_all[:,0,0],label="$S=1$",lw=3)
plt.plot(t1_l,E_all[:,1,0],label="$S=0$",lw=3)
#plt.plot(t1_l,E_gutz_all,label="Gutzwiller")
plt.legend(fontsize=30)
plt.xlabel("$t\'$",fontsize=40)
plt.ylabel("$E_0$",fontsize=40)
#plt.title(f"$N={Lx}\\times{Ly}$",fontsize=40)#,bc={bc}")
plt.xticks([-1.5,-1,-0.5,0],fontsize=30)
plt.yticks([-6,-5.5,-5,-4.5],fontsize=30)
plt.tight_layout()
plt.savefig(f"plot/E0_1_{Lx}_{Ly}_{bc}.pdf",format='pdf')

'''plt.figure()
plt.plot(t1_l,psipsi_all)
plt.xlabel("$t_1$")
plt.ylabel(r"$|\langle \psi_\mathrm{Gutz}|\psi_0\rangle|$")'''
plt.show()

