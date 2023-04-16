import numpy as np
import matplotlib.pyplot as plt
import sys

Nx = int(sys.argv[1])
m = int(sys.argv[2])
chi = 1000

if m % 2 == 1:
    Sz = "1/2"
else:
    Sz = "0"

#t_l = np.linspace(-3,3,31)
t_l = [-3.0,-2.0,-1.5,-1,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,
       0.0,0.1,0.2,0.3,0.4,0.5,0.8,1.0,1.5,2.0,3.0]

tc_dic = {(4,2):[-0.35,-2],(4,4):[-0.25,-0.45],(4,12):[-0.25,-0.55],
          (4,28):[0.05,-2],(4,29):[0.05,-2],
          (6,2):[-0.35,-2],(6,4):[-0.25,-0.35],(6,8):[-0.25,-0.45],
          (6,12):[-0.15,-0.45],(6,16):[-0.1,-0.4],
          (6,20):[-0.1,-0.4],(6,40):[0.05,-2],(6,41):[0.05,-2]}

E_F = np.fromfile(f"data/E_3l_N_{Nx}_n_{m}_p_F_chi_{chi}.dat")
E_A = np.fromfile(f"data/E_3l_N_{Nx}_n_{m}_p_A_chi_{chi}.dat")
#E_A2 = np.fromfile(f"data/E_ex_N_{Nx}_n_{m}_p_A_chi_{chi}.dat")
S0_F = np.fromfile(f"data/S_3l_N_{Nx}_n_{m}_p_F_chi_{chi}.dat")
S0_A = np.fromfile(f"data/S_3l_N_{Nx}_n_{m}_p_A_chi_{chi}.dat")
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
plt.savefig(f"plot/E0_3l_Nx_{Nx}_n_{m}_chi_{chi}.pdf",format='pdf')

plt.figure()
plt.plot(t_l,S_F,'.-',label="$S_z^{tot}=N/2$")
plt.plot(t_l,S_A,'.-',label=f"$S_z^{{tot}}={Sz}$")
#plt.plot(t_l,S_A2,'.-',label=f"$S_z^{{tot}}={Sz}$")
if (Nx,m) in tc_dic:
    t1,t2 = tc_dic[(Nx,m)]
    plt.plot([t1,t1],[0,S_F[0]],'--',color='grey')
    plt.plot([t2,t2],[0,S_F[0]],'--',color='grey')
plt.legend()
plt.xlabel("$t_1$")
plt.ylabel("$S^{tot}$")
plt.savefig(f"plot/S0_3l_Nx_{Nx}_n_{m}_chi_{chi}.pdf",format='pdf')
plt.show()

