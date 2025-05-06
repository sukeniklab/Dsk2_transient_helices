import numpy as np 
import matplotlib.pyplot as plt

mods = ["",  "_bound"] 
colors = ['blue',  'green']
directory = [ "WT_Unbound",  "WT_Bound"]

labels = [ "Unbound Simulation",  "Bound Simulation"]


trials = list(range(1, 11))


bins = np.arange(0, 85, 0.75)
bin_centers = 0.5 * (bins[1:] + bins[:-1])
fig, ax = plt.subplots(figsize=(5, 4))

for j,th in enumerate(directory): 

    avg_rg = []
    histo_rg = []

    for trial in trials:
        rg = 10 * np.load(f"../../Dsk2/Data/WT{mods[j]}/CALVADOS3COM_2.0_MD_gpu_trial{trial}_Dsk2_full{mods[j]}/Dsk2_full{mods[j]}/0/Rg_traj.npy")


        h, bins_edge = np.histogram(rg, bins=bins)

        histo = h/ np.sum(h)

        # ax.stairs(histo, bins_edge, color=colors[j], alpha=0.2)

        histo_rg.append(histo)
        avg_rg.append(np.mean(rg))
    
    mean = np.mean(histo_rg, axis=0)
    std = np.std(histo_rg, axis=0)

    ##Plot Distribution 
    ax.plot(bin_centers, mean, color=colors[j],label=labels[j])
    ax.fill_between(bin_centers,mean-std, mean+std, color=colors[j], alpha = 0.2, zorder = 50, linewidth=0)

    ax.fill_between(bin_centers, mean, color=colors[j], alpha=0.1)

    ax.plot([np.mean(avg_rg),np.mean(avg_rg)], [0, 0.1], color=colors[j], linestyle='--')

    ax.fill_betweenx(y = [0, 0.1], x1=np.mean(avg_rg)-np.std(avg_rg) , x2=np.mean(avg_rg)+np.std(avg_rg), color=colors[j], alpha=0.3, zorder = 50)




    print(f"{labels[j]}", "avg mean:", np.mean(avg_rg), "avg mean:", np.std(avg_rg))

ax.plot([35.1, 35.1], [0, 0.1], color='fuchsia', linestyle='--')

ax.fill_betweenx(y = [0, 0.1], x1=35.1-0.2 , x2=35.1+0.2, color='fuchsia', alpha=0.3, zorder = 50)

ax.plot([35.1, 35.1], [0, -0.1], color='fuchsia', label='Experimental')

ax.plot([35.1, 35.1], [0, 0.1], color='black', linestyle='--', label=r'Avgerage $R_g$')



fig.legend(bbox_to_anchor=(0.9, 0.89))
ax.set_ylabel('Probability', fontsize=12)
ax.set_xlabel('$R_g$ ($\AA$)', fontsize=12)
ax.set_xlim(20,80)
ax.set_ylim(0,0.065)
ax.set_yticks([0.0, 0.02, 0.04, 0.06])
plt.savefig('Rg_unbound_v_bound.png', dpi=300, bbox_inches="tight")
plt.savefig("Rg_unbound_v_bound.svg", dpi=300, format="svg", bbox_inches="tight")
plt.show()


