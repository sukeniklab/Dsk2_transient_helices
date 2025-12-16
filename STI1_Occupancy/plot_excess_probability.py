import numpy as np 
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
plt.rcParams['svg.fonttype'] = 'none'
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import ScalarFormatter
from scipy.interpolate import interp1d, UnivariateSpline


def smooth_log_histogram(log_histogram, smooth):

    # Find valid (non-NaN) points
    valid_mask = ~np.isnan(log_histogram)
    
    if np.sum(valid_mask) < 4:  # Need minimum points for spline
        print("Warning: Too few valid points for spline interpolation")
        return log_histogram
    
    # Create x-coordinates for valid points
    x_coords = np.arange(len(log_histogram))
    valid_x = x_coords[valid_mask]
    valid_y = log_histogram[valid_mask]
    
    # Fit spline to valid points
    try:
        spline = UnivariateSpline(valid_x, valid_y, s=smooth, k=3)
        # spline = InterpolatedUnivariateSpline(valid_x, valid_y)
        # Evaluate spline on all points
        smoothed = np.full_like(log_histogram, np.nan)
        
        # Only interpolate within the range of valid data (don't extrapolate)
        min_valid, max_valid = valid_x.min(), valid_x.max()
        interp_mask = (x_coords >= min_valid) & (x_coords <= max_valid)
        smoothed[interp_mask] = spline(x_coords[interp_mask])
        
        return smoothed
    except:
        print("Spline fitting failed, returning original data")
        return log_histogram

numReplicas=20
sim_time = 70 #ns
print_interval = 0.01

snaps = numReplicas * int((sim_time/print_interval) - (0.05 * (sim_time/print_interval)))


trials = list(range(1, 11))
domains = [[75,145], [223,325]]

mods = ["",  "bound_"]
colors = ['blue',  'green']
labels = [ "Unbound",  "Bound"]
directory = [ "WT_unbound",  "WT_bound"]

STI1_locs = [[146,222]]

### set up figure and insets ###
fig, axes = plt.subplots(1, 2, figsize=(10, 4))


for j,th in enumerate(directory): 


    for i, section in enumerate(domains): 

        all_histograms = np.zeros((len(trials), len(np.arange(section[0], section[1]+1, 1))-4))
        all_histograms_ex = np.zeros((len(trials), len(np.arange(section[0], section[1]+1, 1))))
        

        all_histograms_nolog = np.zeros((len(trials), len(np.arange(section[0], section[1]+1, 1))))
        for trial in trials:
            print(f"{th} {section} {trial}")
            try: 
                ## if file found, plot it 
                histogram = np.load(f"Histos/{th}/Dsk2_full_{mods[j]}trial{trial}_innerVol_idr{section[0]}_{section[1]}.npy")

                print(f"{th}/Dsk2_full_{mods[j]}trial{trial}_innerVol_idr{section[0]}_{section[1]}")
                print('att', histogram)
            
                histogram /= snaps

                histogram_ex = np.load(f"Histos/WT_unbound_excluded_volume/Dsk2_full_excluded_volume_trial{trial}_innerVol_idr{section[0]}_{section[1]}.npy")
                
                log_ex = np.log10(histogram_ex / (10* snaps))
                log_ex[np.isinf(log_ex)] = np.nan
                

                if i ==0: 
                    smoothed_log = smooth_log_histogram(log_ex[:-4], 1)  # Adjust smoothing_factor as needed
                    smooth_ev = 10**smoothed_log
                    ratio = histogram[:-4] / smooth_ev


                elif i ==1:
                    smoothed_log = smooth_log_histogram(log_ex[4::],2.4)
                    smooth_ev = 10**smoothed_log
                    ratio = histogram[4::] / smooth_ev

                ratio[np.isinf(ratio)] = np.nan

                all_histograms[trial-1, :] = ratio

                

                
            except FileNotFoundError:
                print(f'missing Dsk2_full_{mods[j]}trial{trial}_innerVol_idr{section[0]}_{section[1]}')
                continue
      


        all_histograms = all_histograms[~np.all(all_histograms == 0, axis=1)] ##ensure that no arrays of all 0s are still left --- check if exception!
        
        print(all_histograms)
        mean_histogram = np.nanmean(all_histograms, axis=0)



        
        ### take avg and get std between independent runs ###
        std_histogram = np.nanstd(all_histograms, axis=0)
        
        residues = np.arange(section[0], section[1]+1, 1)
        axes[i].set_ylim(0, 250)


        
        if i ==0:  ### if first IDR region           

            diff_STI1 = section[0] - STI1_locs[0][0]
            dist_STI1 = np.arange(diff_STI1, 0, 1)

            print(dist_STI1[0], dist_STI1[-1])
            axes[i].set_xlim(dist_STI1[0], -5) 

            ## main plot
            axes[i].plot(dist_STI1[:-4], mean_histogram, color=colors[j],  zorder = 100)
            axes[i].scatter(dist_STI1[:-4], mean_histogram, color=colors[j], label=labels[j],  s=10, zorder = 100)


            axes[i].fill_between(dist_STI1[:-4], mean_histogram - std_histogram, mean_histogram + std_histogram,
                     color=colors[j], alpha=0.15, zorder=50,linewidth=0)

            axes[i].text(120-STI1_locs[0][0], 237, "TH1", fontsize=12, color='black')
            axes[i].grid(color='grey', linestyle='-', linewidth=0.25, alpha=0.5)

        
            ##plot upper residue number
            secax0 = axes[i].secondary_xaxis('top', functions=(lambda x: 146 + x , lambda x: x - 146   ))
            x_coords = secax0.get_xticks()

        
           
        elif i == 1: ### if part of second idr
            
            diff_STI1 = section[1] - STI1_locs[0][1]
            dist_STI1 = np.arange(0,diff_STI1, 1)
            axes[i].set_xlim(5, dist_STI1[-1]) 

            ## main plot 
            axes[i].plot(dist_STI1[4::], mean_histogram, color=colors[j], zorder = 100)
            axes[i].scatter(dist_STI1[4::], mean_histogram, color=colors[j],s=10, zorder = 100)

            axes[i].fill_between(dist_STI1[4::], mean_histogram - std_histogram, mean_histogram + std_histogram,
                     color=colors[j], alpha=0.15, zorder=50,linewidth=0)
            axes[i].text(279- STI1_locs[0][1], 237, "TH2", fontsize=12, color='black')
            axes[i].text(302- STI1_locs[0][1], 237, "TH3", fontsize=12, color='black')
            axes[i].grid(color='grey', linestyle='-', linewidth=0.25, alpha=0.5)
            axes[i].set_xticks([10, 20, 30, 40, 50, 60, 70, 80, 90, 100]) # to match the inset


            ##set upper residue number 
            secax1 = axes[i].secondary_xaxis('top', functions=(lambda x: x + 223, lambda x: x+325))
            


### mark TH regions 
axes[0].axvspan(-13, -33, color='orange', alpha=0.3)
axes[1].axvspan(56, 68, color='orange', alpha=0.3)
axes[1].axvspan(80, 90, color='orange', alpha=0.3)   



axes[0].get_position()


fig.legend( bbox_to_anchor=(1.05, 0.55))
fig.supylabel("Excess Probability ($P_{\mathrm{g}}$/$P_{\mathrm{g,EV}}$)", x=0.05) 
fig.supxlabel("Distance from STI1", y=-0.02)    
fig.text(0.5, 1.0, "Residue Number", ha='center', va='center', fontsize='large')  # Secondary x-axis label

plt.savefig('Dsk2_full_bound_v_unbound_withEV_ratio_tounbound_v2.png', dpi=300, bbox_inches="tight")
plt.savefig("Dsk2_full_bound_v_unbound_ratio_v2.svg", dpi=300, format="svg", bbox_inches="tight")
plt.show()
