import numpy as np 
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
plt.rcParams['svg.fonttype'] = 'none'
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import ScalarFormatter


numReplicas=20
sim_time = 70 #ns
print_interval = 0.01

snaps = numReplicas * int((sim_time/print_interval) - (0.05 * (sim_time/print_interval)))


trials = list(range(1, 11))
domains = [[75,145], [223,325]]

mods = ["",  "bound_"] 
colors = ['blue',  'green']
labels = [ "Unbound",  "Bound"]
directory = [ "WT_Unbound",  "WT_Bound"]

STI1_locs = [[146,222]]

### set up figure and insets ###
fig, axes = plt.subplots(1, 2, figsize=(10, 4))
ax_inset = inset_axes(axes[0], width="30%", height="30%", loc="lower right", borderpad=0)  
ax_inset.set_axes_locator(lambda ax, r: [0.335, 0.17, 0.125, 0.26]) 

ax_inset2 = inset_axes(axes[1], width="30%", height="30%", loc="lower left") 
ax_inset2.set_axes_locator(lambda ax, r: [0.616, 0.17, 0.125, 0.26]) 


for j,th in enumerate(directory): 
    for i, section in enumerate(domains): 

        all_histograms = np.zeros((len(trials), len(np.arange(section[0], section[1]+1, 1))))

        all_histograms_nolog = np.zeros((len(trials), len(np.arange(section[0], section[1]+1, 1))))
        for trial in trials:
            print(f"{th} {section} {trial}")
            try: 
                ## if file found, plot it 
                histogram = np.load(f"Histos/{th}/Dsk2_full_{mods[j]}trial{trial}_innerVol_idr{section[0]}_{section[1]}.npy")
            
                histogram /= snaps
 
                all_histograms[trial-1, :] = np.log10(histogram) 
                all_histograms_nolog[trial-1, :] = (histogram)

            except Exception as e:
                print(f'missing Dsk2_full_{mods[j]}trial{trial}_innerVol_idr{section[0]}_{section[1]}')
                continue
      


        all_histograms = all_histograms[~np.all(all_histograms == 0, axis=1)] ##ensure that no arrays of all 0s are still left --- check if exception!

        ### take avg and get std between independent runs ###
        mean_histogram = np.mean(all_histograms, axis=0)
        std_histogram = np.std(all_histograms, axis=0)
        std_nolog = np.std(all_histograms_nolog, axis=0)
        
        residues = np.arange(section[0], section[1]+1, 1)
        axes[i].set_ylim( -5, -1.5)


        
        if i ==0:  ### if first IDR region           

            diff_STI1 = section[0] - STI1_locs[0][0]
            dist_STI1 = np.arange(diff_STI1, 0, 1)

            print(dist_STI1[0], dist_STI1[-1])
            axes[i].set_xlim(dist_STI1[0], -5) 

            ## main plot
            axes[i].plot(dist_STI1, mean_histogram, color=colors[j], label=labels[j], zorder = 100)
            axes[i].fill_between(dist_STI1, mean_histogram - std_histogram, mean_histogram + std_histogram,
                     color=colors[j], alpha=0.15, zorder=50,linewidth=0)

            axes[i].text(120-STI1_locs[0][0], -1.75, "TH1", fontsize=12, color='black')
            axes[i].grid(color='grey', linestyle='-', linewidth=0.25, alpha=0.5)

          
            ## inset 
            ax_inset.plot(dist_STI1, 1e2 * 10**mean_histogram, color=colors[j])
            ax_inset.fill_between(dist_STI1,
                      1e2 * 10**mean_histogram - 1e2 * std_nolog,
                      1e2 * 10**mean_histogram + 1e2 * std_nolog,
                      color=colors[j], alpha=0.15, zorder=50,linewidth=0)

            ax_inset.set_xlim(-40, -10)
            ax_inset.set_ylim(0.0,  0.75)
            ax_inset.set_ylabel("$P_{\mathrm{g}}$",labelpad=-1.5)
            ax_inset.set_ylabel(r"$P_{\mathrm{g}}$ ($\times10^{-2}$)", labelpad=1.5)
            ax_inset.xaxis.set_minor_locator(MultipleLocator(5))
            ax_inset.yaxis.set_minor_locator(MultipleLocator(0.5))
            ax_inset.grid(color='grey', linestyle='-', linewidth=0.25, which='both', alpha=0.5)

            ##plot upper residue number
            secax0 = axes[i].secondary_xaxis('top', functions=(lambda x: 146 + x , lambda x: x - 146   ))
            x_coords = secax0.get_xticks()

        
           
        elif i == 1: ### if part of second idr
            
            diff_STI1 = section[1] - STI1_locs[0][1]
            dist_STI1 = np.arange(0,diff_STI1, 1)
            axes[i].set_xlim(5, dist_STI1[-1]) 

            ## main plot 
            axes[i].plot(dist_STI1, mean_histogram, color=colors[j], zorder = 100)
            axes[i].fill_between(dist_STI1, mean_histogram - std_histogram, mean_histogram + std_histogram,
                     color=colors[j], alpha=0.15, zorder=50,linewidth=0)
            axes[i].text(279- STI1_locs[0][1], -1.75, "TH2", fontsize=12, color='black')
            axes[i].text(302- STI1_locs[0][1], -1.75, "TH3", fontsize=12, color='black')
            axes[i].grid(color='grey', linestyle='-', linewidth=0.25, alpha=0.5)
            axes[i].set_xticks([10, 20, 30, 40, 50, 60, 70, 80, 90, 100]) # to match the inset

            ##inset
            ax_inset2.plot(dist_STI1, 1e3 * 10**mean_histogram, color=colors[j])
            ax_inset2.fill_between(dist_STI1,
                                1e3 * 10**mean_histogram - 1e3 * std_nolog,
                                1e3 * 10**mean_histogram + 1e3 * std_nolog,
                                color=colors[j], alpha=0.15, zorder=50,linewidth=0)

            ax_inset2.set_xlim(50, 100)
            ax_inset2.set_ylim(0.0, 0.75)
            ax_inset2.set_ylabel(r"$P_{\mathrm{g}}$ ($\times10^{-3}$)", labelpad=1.5)            
            ax_inset2.xaxis.set_minor_locator(MultipleLocator(5))
            ax_inset2.yaxis.set_minor_locator(MultipleLocator(0.5))
            ax_inset2.grid(color='grey', linestyle='-', linewidth=0.25, which='both', alpha=0.5)


            ##set upper residue number 
            secax1 = axes[i].secondary_xaxis('top', functions=(lambda x: x + 223, lambda x: x+325))
            


### mark TH regions 
axes[0].axvspan(-13, -33, color='orange', alpha=0.3)
ax_inset.axvspan(-13, -33, color='orange', alpha=0.3)

axes[1].axvspan(56, 68, color='orange', alpha=0.3)
axes[1].axvspan(80, 90, color='orange', alpha=0.3)   

ax_inset2.axvspan(56, 68, color='orange', alpha=0.3)
ax_inset2.axvspan(80, 90, color='orange', alpha=0.3)   

axes[0].get_position()


fig.legend( bbox_to_anchor=(0.25, 0.88))
fig.supylabel("log($P_{\mathrm{g}}$)", x=0.05) 
fig.supxlabel("Distance from STI1", y=-0.02)    
fig.text(0.5, 1.0, "Residue Number", ha='center', va='center', fontsize='large')  # Secondary x-axis label

# plt.savefig('Dsk2_full_bound_v_unbound.png', dpi=300, bbox_inches="tight")
plt.savefig("Dsk2_full_bound_v_unbound.svg", dpi=300, format="svg", bbox_inches="tight")
plt.show()
