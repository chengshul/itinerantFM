import numpy as np
import matplotlib.pyplot as plt

n_3_6 = np.array([2,4,8,12,16,20,40,41])/42

t1_3_6 = np.array([[-0.35,-2],[-0.25,-0.35],[-0.25,-0.45],[-0.15,-0.45],
                   [-0.1,-0.4],[-0.1,-0.4],[0.05,-2],[0.05,-2]])

for j in range(1,len(n_3_6)):
    n = n_3_6[j]
    t = t1_3_6[j]
    #if j in [1,5,6]:
    #    plt.plot([n,n],t,"o-",color='C0',lw=2)
    #else:
    plt.plot([n,n],t,"o",color='C0',lw=2)

#plt.plot(n_3_6[1:6],t1_3_6[1:6,0],color='C0',lw=2)
#plt.plot(n_3_6[1:6],t1_3_6[1:6,1],color='C0',lw=2)
#plt.plot(n_3_6[6:],t1_3_6[6:,0],color='C0',lw=2)

plt.fill(list(n_3_6[1:6])+list(n_3_6[5:0:-1]),list(t1_3_6[1:6,0])+list(t1_3_6[5:0:-1,1]),color='wheat')
plt.fill([n_3_6[6],1,1,n_3_6[6]],[t1_3_6[6,0],t1_3_6[6,0],-1,-1],color='wheat')

plt.plot([0,1],[-0.25,-0.25],'--',color='grey')
#plt.text(0.8,-0.3,"Lifshitz")
plt.xlabel("$n$",fontsize=24)
plt.ylabel("$t\'$",fontsize=24)
plt.xticks(fontsize=20)
plt.yticks([-1,-0.6,-0.2,0.2],fontsize=20)
plt.ylim([-1,0.2])
#plt.legend(fontsize=18,loc=[0.52,0.02])
plt.tight_layout()
plt.savefig("pd_honeycomb.pdf",format='pdf')

plt.show()
