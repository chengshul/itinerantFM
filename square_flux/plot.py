import numpy as np
import matplotlib.pyplot as plt

Nx = 10
m_l = np.arange(2,4*Nx)
n = m_l/4/Nx
chi = 1400

E_l = []
S_l = []

for flux in range(1,6):
    E = np.fromfile(f"data/E_fx_N_{Nx}_flux_{flux}_chi_{chi}.dat")
    S = np.fromfile(f"data/S_fx_N_{Nx}_flux_{flux}_chi_{chi}.dat")
    
    E = E.reshape((2,-1))
    S = S.reshape((2,-1))
    S = (np.sqrt(1+4*S)-1)/2
        
    E_l.append(E[1])
    S_l.append(S[1])

marker_list = ['.-']*5
label_list = ["$\phi=0$","$\phi=\pi$", "$\phi=2\pi/3$","$\phi=\pi/2$", "$\phi=2\pi/5$"]
fig, axes = plt.subplots(5, sharex=True) #, sharey=True)
for i in range(5):
    axes[i].plot(n, S_l[i]/(m_l/2), marker_list[i], label=label_list[i])
    axes[i].set_yticks([0,1],['',''],fontsize=30)
    #axes[i].legend(loc='center right',fontsize=30)
axes[2].set_ylabel("$S_{\mathrm{tot}}/S_{\mathrm{tot,max}}$",fontsize=40)
axes[0].set_yticks([0,1],['','1'],fontsize=30)
axes[4].set_yticks([0,1],[0,''],fontsize=30)
plt.xlabel("$n$",fontsize=40)
plt.xticks(fontsize=30)
#plt.yticks(fontsize=30)
plt.tight_layout()
plt.subplots_adjust(hspace=0)
plt.savefig(f"plot/S_square_all_chi_{chi}.pdf",format='pdf')

plt.show()
