
import numpy as np
import matplotlib.pyplot as plt




pos1 = (np.loadtxt("datos/TSI_2D_10_pos1.txt", delimiter='\t', dtype=float))
vel1 = (np.loadtxt("datos/TSI_2D_10_vel1.txt", delimiter='\t', dtype=float))

pos2 = (np.loadtxt("datos/TSI_2D_10_pos2.txt", delimiter='\t', dtype=float))
vel2 = (np.loadtxt("datos/TSI_2D_10_vel2.txt", delimiter='\t', dtype=float))

E = np.loadtxt("datos/TSI_2D_10_E.txt", delimiter='\t', dtype=float)


fig = plt.figure(figsize=(4, 3))
ax = fig.subplots()

ax.tick_params(direction='in')


plt.scatter(pos1[0], vel1[0], color='c', s=1)
plt.scatter(pos2[0], vel2[0], color='r', s=1)

plt.ylim(-15,15)
plt.xlim(0,100)


for i in list(range(0,100,10)):
    
    plt.scatter(pos1[i], vel1[i], color='c', s=1) 
    plt.scatter(pos2[i], vel2[i], color='r', s=1) 
    plt.savefig('graficas/TSI_2D_phaseSpace10_' + str(i) + '.png',bbox_inches="tight")
    


plt.cla()

plt.xlabel("$x/\lambda_D$")
plt.ylabel("$y/\lambda_D$")
plt.ylim(0,10)
plt.xlim(0,100)
plt.pcolormesh(np.transpose(E))
plt.colorbar()
plt.savefig('graficas/TSI_2D_E_10.png',bbox_inches="tight")

