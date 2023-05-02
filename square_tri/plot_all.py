import numpy as np
import matplotlib.pyplot as plt

n1 = 2/np.array([4,6,8,9,12,16,20,24,25,32,36])
t1 =      [-0.82,-0.73,-0.69,-0.65,-0.61,-0.59,-0.57,-0.56,-0.56,-0.53,-0.54]
n1p = 2/np.array([20,25,30,36,42,48])
t1p = [-0.86,-0.76,-0.71,-0.66,-0.64,-0.63]
Nd_4 = np.array([2,6,8,10,12,14,15])
Nd_6 = np.array([2,5,6,7,10,12,16,18,20,23])
Nd_8 = np.array([2,7,8,9,10,12,16,20,24,28,31])
Nd_6_2 = np.array([2,4,6,8,10,11])
Nd_8_2 = np.array([2,4,6,8,10,12,14,15])
Nd_10_2 = np.array([2,4,6,8,10,12,14,16,18,19])

#n2 = [0.125,0.25,0.5,10/16,12/16,14/16,15/16,14/24,22/24,23/24,31/32]
#t2 = [-0.8,-0.9,-0.6,-0.3,-0.3,-0.3,0.3,-0.3,0.2,0.3,0.3]
t1_4 = [-0.55,-0.85,-0.55,-0.35,-0.25,-0.45,0.25]
t1_6 = [-0.55,-1.25,-0.9,-0.95,-0.65,-0.55,-0.25,-0.25,-0.05,0.25]
t1_8 = [-0.55,-1.45,-1.25,-1.35,-1.35,-0.75,-0.5,-0.25,-0.15,0.15,0.25]
t1_6_2 = [-0.65,-0.85,-0.85,-0.65,0.1,0.25]
t1_8_2 = [-0.55,-0.85,-0.85,-0.75,-0.6,-0.2,0.25,0.25]
t1_10_2 = [-0.55,-0.85,-0.85,-0.75,-0.75,-0.65,-0.3,0.1,0.25,0.25]

t_all = [-1.5,-0.5]+t1[-1::-9]+t1_8_2[-5:]+[t1_8_2[-1],-1.5]
n_all = [0,0]+list(n1[-1::-9])+list(Nd_8_2[-5:]/16)+[1,1]

plt.fill(n_all,t_all,color='wheat')

plt.plot(n1,t1,'o',label="2 particles") #, OBC")
#plt.plot(n1p,t1p,'x',label="2 particles, PBC")

#plt.plot(Nd_4/16,t1_4,'.',color='C1',label=r"$4\times4$")
plt.plot(Nd_6/24,t1_6,'+',color='C1',label=r"$4\times6$")
plt.plot(Nd_8/32,t1_8,'x',color='C2',label=r"$4\times8$")
#plt.plot(Nd_6_2/12,t1_6_2,'d',color='C4',label=r"$2\times6$")
plt.plot(Nd_8_2/16,t1_8_2,'s',color='C3',label=r"$2\times8$")
plt.plot(Nd_10_2/20,t1_10_2,'^',color='C4',label=r"$2\times10$")
#plt.text(N2_i[2]/N2_i[0]/N2_i[1],t2_i+loc[i],
#         f"${N2_i[0]}\\times{N2_i[1]},{N2_i[2]}$")#,label="DMRG")
plt.legend(fontsize=16, loc="upper left", ncol=2)

plt.plot([1,1],[-1.6,1.6],'--',color='grey')
plt.plot([0,2],[0,0],'--',color='grey')
plt.plot([0,2],[1,1],'--',color='grey')
plt.plot([0,2],[-1,-1],'--',color='grey')
#plt.text(1,1.2,"half filling",color='grey')
#plt.text(1.5,0.05,"square",color='grey')
#plt.text(1.5,1.05,"triangle",color='grey')
#plt.text(1.5,-0.95,"triangle",color='grey')

plt.xlabel("$n$",fontsize=24)
plt.ylabel("$t\'$",fontsize=24)
plt.xlim([-0.1,1.1])
plt.ylim([-1.5,1.1])
plt.xticks([0,0.25,0.5,0.75,1],fontsize=20)
plt.yticks(fontsize=20)
plt.tight_layout()
plt.savefig("pd.pdf",format='pdf')

'''
plt.figure()
plt.fill([-0.5,1,1,-0.5],[0,0,1.5,1.5],fill=False,hatch='//',ec='C0')
plt.fill([1,2.5,2.5,1],[-1.5,-1.5,0,0],fill=False,hatch='//',ec='C0')
plt.fill([-0.5,1,1,-0.5],[-1.5,-1.5,0,0],fill=False,hatch=r'\\',ec='C1')
plt.fill([1,2.5,2.5,1],[0,0,1.5,1.5],fill=False,hatch=r'\\',ec='C1')

plt.plot([1,1],[-1.6,1.6],'--',lw=3,color='grey')
plt.plot([0,2],[0,0],'--',lw=3,color='grey')
plt.plot([0,2],[1,1],'--',lw=3,color='grey')
plt.plot([0,2],[-1,-1],'--',lw=3,color='grey')

plt.text(1.1,1.15,"half filling",color='grey',backgroundcolor='w',fontsize=20)
plt.text(1.6,-0.25,"square",color='grey',backgroundcolor='w',fontsize=20)
plt.text(1.6,0.75,"triangle",color='grey',backgroundcolor='w',fontsize=20)
plt.text(1.6,-1.25,"triangle",color='grey',backgroundcolor='w',fontsize=20)

plt.plot([1.3,0.7],[0.8,-0.8],'o',color='C2')
plt.arrow(1.2,0.55,-0.4,-1.1,width=0.01,head_width=0.05,head_length=0.16,color='k')
plt.arrow(0.8,-0.55,0.4,1.1,width=0.01,head_width=0.05,head_length=0.16,color='k')

plt.xlabel("$n$",fontsize=24)
plt.ylabel("$t\'$",fontsize=24)
plt.xlim([-0.1,2.1])
plt.ylim([-1.4,1.4])
plt.xticks([0,0.5,1,1.5,2],fontsize=20)
plt.yticks(fontsize=20)
plt.tight_layout()
#plt.savefig("pd0.pdf",format='pdf')

plt.figure()
plt.plot(n1,t1,'x',label="2 particles") #, OBC")
#plt.plot(n1p,t1p,'x',label="2 particles, PBC")

plt.plot([1,1],[-1.2,1.2],'--',color='grey')
plt.plot([0,1],[0,0],'--',color='grey')
plt.plot([0,1],[1,1],'--',color='grey')
plt.plot([0,1],[-1,-1],'--',color='grey')
plt.text(0.9,1.2,"half filling",color='grey')
plt.text(0.8,0.05,"square",color='grey')
plt.text(0.8,1.05,"triangle",color='grey')
plt.text(0.8,-0.95,"triangle",color='grey')
plt.legend()

plt.xlabel("$n$")
plt.ylabel("$t\'$")
#plt.xlim([0,1.25])
#plt.ylim([-1.5,0.5])
plt.savefig("pd_2p.pdf",format='pdf')
'''

plt.show()
