import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as an


pos1 = (np.loadtxt("datos/TSI_1D_pos1.txt", delimiter='\t', dtype=float))
vel1 = (np.loadtxt("datos/TSI_1D_vel1.txt", delimiter='\t', dtype=float))

pos2 = (np.loadtxt("datos/TSI_1D_pos2.txt", delimiter='\t', dtype=float))
vel2 = (np.loadtxt("datos/TSI_1D_vel2.txt", delimiter='\t', dtype=float))

fig = plt.figure(figsize=(4, 3))
ax = fig.subplots()

ax.tick_params(direction='in')

plt.scatter(pos1[0], vel1[0], color='r', s=1)
plt.scatter(pos2[0], vel2[0], color='c', s=1)

plt.xlabel("$x/\lambda_D$")
plt.ylabel("$v/v_{th}$")

plt.ylim(-15,15)
plt.xlim(0,100)
plt.savefig('graficas/TSI_1D_phaseSpace_0.png',bbox_inches="tight")
def update(i):
    
    plt.cla()
    
    plt.ylim(-15,15)
    plt.xlim(0,100)
    
    plt.xlabel("$x/\lambda_D$")
    plt.ylabel("$v/v_{th}$")
    
    
    plt.scatter(pos1[i], vel1[i], color='r', s=1) 
    plt.scatter(pos2[i], vel2[i], color='c', s=1)
    if(i%10 == 0): plt.savefig('graficas/TSI_1D_phaseSpace_' + str(i) + '.png',bbox_inches="tight")



ani = an.FuncAnimation(fig, update, frames=100)

writergif = an.PillowWriter(fps=60) 
ani.save('graficas/TSI_1D_phaseSpace_ani.gif', writer=writergif)

##############################################################################

potencial  = (np.loadtxt("datos/TSI_1D_potencial.txt" ,delimiter='\t', dtype=float))
potencialK  = (np.loadtxt("datos/TSI_1D_pot_K.txt" ,delimiter='\t', dtype=float))
cinetica = (np.loadtxt("datos/TSI_1D_cinetica.txt" ,delimiter='\t', dtype=float))
x = list(range(0, 100))

plt.cla()


for i in [0,20,50]:
    plt.xlabel("$x/\lambda_D$")
    plt.ylabel("Potencial eléctrico normalizado")
    
    plt.ylim(-15,15)
    plt.xlim(0,100)
    potL = []
    for j in range(len(potencial[i])):
        if(j%8 == 0): 
            potL.append(potencial[i][j])
    potL.pop()
    plt.plot(x, potL, color='r')
    plt.savefig('graficas/TSI_1D_potencial_' + str(i) + '.png', bbox_inches="tight")
    plt.cla()
    



t = list(range(101))

Ue = [0.0]
for i in range(len(potencialK)):
    Ue.append(potencialK[i])
    
plt.xlabel("$t\omega_r$")
plt.ylabel("$U/K_0$")

plt.plot(t, Ue, color='r')
plt.savefig('graficas/TSI_1D_UvsT.png', bbox_inches="tight")
plt.cla()

plt.xlabel("$t\omega_r$")
plt.ylabel("$K/K_0$")

plt.plot(t, cinetica, color='c')
plt.savefig('graficas/TSI_1D_pKvsT.png', bbox_inches="tight")

plt.show()



