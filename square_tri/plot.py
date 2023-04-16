import numpy as np
import matplotlib.pyplot as plt
import sys

Nx = int(sys.argv[1])
m = int(sys.argv[2])
chi = 1400

if m % 2 == 1:
    Sz = "1/2"
else:
    Sz = "0"

#t_l = np.linspace(-3,3,31)
t_l = [-3.0,-2.0,-1.5,-1.4,-1.3,-1.2,-1.1,-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,
       -0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.8,1.0,1.5,2.0,3.0]

tc_dic = {(4,2):-0.55,(4,6):-0.85,(4,8):-0.55,(4,10):-0.35,
          (4,12):-0.25,(4,14):-0.45,(4,15):0.25,
          (6,2):-0.55,(6,5):-1.25,(6,6):-0.9,(6,7):-0.95,(6,10):-0.65,
          (6,12):-0.55,(6,16):-0.25,(6,18):-0.25,(6,20):-0.05,(6,23):0.25,
          (8,2):-0.55,(8,7):-1.45,(8,8):-1.25,(8,9):-1.35,(8,10):-1.35,
          (8,12):-0.75,(8,16):-0.5,(8,20):-0.25,(8,24):-0.15,(8,28):0.15,
          (8,31):0.25}

E_F = np.fromfile(f"data/new/E_N_{Nx}_n_{m}_p_F_chi_{chi}.dat")
E_A = np.fromfile(f"data/new/E_N_{Nx}_n_{m}_p_A_chi_{chi}.dat")
#E_A2 = np.fromfile(f"data/E_ex_N_{Nx}_n_{m}_p_A_chi_{chi}.dat")
S0_F = np.fromfile(f"data/new/S_N_{Nx}_n_{m}_p_F_chi_{chi}.dat")
S0_A = np.fromfile(f"data/new/S_N_{Nx}_n_{m}_p_A_chi_{chi}.dat")
#S0_A2 = np.fromfile(f"data/S_ex_N_{Nx}_n_{m}_p_A_chi_{chi}.dat")
S_F = (np.sqrt(1+4*S0_F)-1)/2
S_A = (np.sqrt(1+4*S0_A)-1)/2
#S_A2 = (np.sqrt(1+4*S0_A2)-1)/2

plt.plot(t_l,E_F,'.-',label="$S_z^{tot}=N/2$")
plt.plot(t_l,E_A,'.-',label=f"$S_z^{{tot}}={Sz}$")
#plt.plot(t_l,E_A2,'.-',label=f"$S_z^{{tot}}={Sz}$")
plt.legend()
plt.xlabel("$t_1$")
plt.ylabel("$E_0$")
plt.savefig(f"plot/E0_Nx_{Nx}_n_{m}_chi_{chi}.pdf",format='pdf')

plt.figure()
plt.plot(t_l[:],S_F[:],'.-',label="$S_z^{tot}=N/2$")
plt.plot(t_l[:],S_A[:],'.-',label=f"$S_z^{{tot}}={Sz}$")
if (Nx,m) in tc_dic:
    tc = tc_dic[(Nx,m)]
    plt.plot([tc,tc],[0,S_F[0]],'--',color='grey')
#plt.plot(t_l,S_A2,'.-',label=f"$S_z^{{tot}}={Sz}$")
plt.legend()
plt.title(f"$N=4\\times{Nx},N={m},\chi={chi}$")
plt.xlabel("$t\'$")
plt.ylabel("$S^{tot}$")
plt.savefig(f"plot/S0_Nx_{Nx}_n_{m}_chi_{chi}.pdf",format='pdf')
plt.show()

