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

tc_dic = {(6,2):-0.65,(6,4):-0.85,(6,6):-0.85,(6,8):-0.65,(6,10):0.1,(6,11):0.25,
          (8,2):-0.55,(8,4):-0.85,(8,6):-0.85,(8,8):-0.75,(8,10):-0.6,
          (8,12):-0.2,(8,14):0.25,(8,15):0.25,
          (10,2):-0.55,(10,4):-0.85,(10,6):-0.85,(10,8):-0.75,(10,10):-0.75,
          (10,12):-0.65,(10,14):-0.3,(10,16):0.1,(10,18):0.25,(10,19):0.25}

E_F = np.fromfile(f"data/new/E_2l_N_{Nx}_n_{m}_p_F_chi_{chi}.dat")
E_A = np.fromfile(f"data/new/E_2l_N_{Nx}_n_{m}_p_A_chi_{chi}.dat")
S0_F = np.fromfile(f"data/new/S_2l_N_{Nx}_n_{m}_p_F_chi_{chi}.dat")
S0_A = np.fromfile(f"data/new/S_2l_N_{Nx}_n_{m}_p_A_chi_{chi}.dat")
S_F = (np.sqrt(1+4*S0_F)-1)/2
S_A = (np.sqrt(1+4*S0_A)-1)/2

plt.plot(t_l,E_F,'.-',label="$S_z^{tot}=N/2$")
plt.plot(t_l,E_A,'.-',label=f"$S_z^{{tot}}={Sz}$")
plt.legend()
plt.xlabel("$t_1$")
plt.ylabel("$E_0$")
plt.savefig(f"plot/E0_2l_Nx_{Nx}_n_{m}_chi_{chi}.pdf",format='pdf')

plt.figure()
plt.plot(t_l,S_F,'.-',label="$S_z^{tot}=N/2$")
plt.plot(t_l,S_A,'.-',label=f"$S_z^{{tot}}={Sz}$")
if (Nx,m) in tc_dic:
    tc = tc_dic[(Nx,m)]
    plt.plot([tc,tc],[0,S_F[0]],'--',color='grey')
plt.legend()
plt.xlabel("$t_1$")
plt.ylabel("$S^{tot}$")
plt.savefig(f"plot/S0_2l_Nx_{Nx}_n_{m}_chi_{chi}.pdf",format='pdf')
plt.show()

