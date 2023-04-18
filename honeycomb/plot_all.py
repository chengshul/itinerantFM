import numpy as np
import matplotlib.pyplot as plt

n = 1/np.array([12,12,15,15,18,18,20,20])
l = [[3,4],[4,3],[3,5],[5,3],[3,6],[6,3],[4,5],[5,4]]
t1 = [-0.35,-0.38,-0.33,-0.36,-0.33,-0.35,-0.32,-0.36]

n_cyl = 1/np.array([9,12,15,18,15,20,25])
l_cyl = [[2,3],[2,4],[2,5],[2,6],[4,3],[4,4],[4,5]]
t1_cyl = [-0.37,-0.35,-0.34,-0.34,-0.34,-0.3,-0.29]

l_obc = [[2,2],[3,2],[4,2],[5,2],[6,2],[2,4],[3,4],[4,4]]
t1_obc = [-0.3,-0.28,-0.27,-0.26,-0.26,-0.32,-0.28,-0.26]
n_obc = []
for l in l_obc:
    lx,ly = l
    n_obc.append(1/(2*(ly+1)*lx+2*ly))

n_3_4_0 = np.array([2,4,12,28,29])/30
n_3_6_0 = np.array([2,4,8,12,16,20,40,41])/42
n_2_4_0 = np.array([2,4,8,12,27])/28
n_2_6_0 = np.array([2,4,8,12,16,20,36,39])/40

n_all = [n_2_4_0,n_2_6_0,n_3_4_0,n_3_6_0]

t1_3_4_0 = [[-0.35,-2],[-0.25,-0.45],[-0.25,-0.55],[0.05,-2],[0.05,-2]]
t1_3_6_0 = [[-0.35,-2],[-0.25,-0.35],[-0.25,-0.45],[-0.15,-0.45],[-0.1,-0.4],
            [-0.1,-0.4],[0.05,-2],[0.05,-2]]
t1_2_4_0 = [[-0.25,-2],[-0.35,-0.65],[-0.3,-0.9],[-0.05,-0.65],[0.05,-2]]
t1_2_6_0 = [[-0.25,-2],[-0.25,-0.65],[-0.35,-0.75],[-0.25,-0.9],[-0.25,-0.55],
            [0.05,-0.65],[-0.05,-2],[0.05,-2]]

t1_all = [t1_2_4_0,t1_2_6_0,t1_3_4_0,t1_3_6_0]

plt.plot(n,t1,'.',color='C0',label="$N=2$")# torus")
plt.plot(n_cyl,t1_cyl,'.',color='C0')#,label="2 particles cylinder")
plt.plot(n_obc,t1_obc,'.',color='C0')#,label="2 particles obc")

label = [r"$2\times 4$",r"$2\times 6$",
         r"$3\times 4$",r"$3\times 6$"]
marker = ['>-','<-','^-','v-']

for i in range(len(n_all)):
    for j in range(len(n_all[i])):
        n = n_all[i][j]
        t = t1_all[i][j]
        if j == 0:
            plt.plot([n,n],t,marker[i],color=f'C{i+1}',label=label[i])
        else:
            plt.plot([n,n],t,marker[i],color=f'C{i+1}')


plt.plot([0,1],[-0.25,-0.25],'--',color='grey')
#plt.text(0.8,-0.3,"Lifshitz")
plt.xlabel("$n$",fontsize=40)
plt.ylabel("$t\'$",fontsize=40)
plt.xticks(fontsize=30)
plt.yticks([-1,-0.6,-0.2,0.2],fontsize=30)
plt.ylim([-1,0.2])
plt.legend(fontsize=18,loc=[0.52,0.02])
plt.tight_layout()
plt.savefig("pd_honeycomb.pdf",format='pdf')

plt.show()
