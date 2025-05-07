import numpy as np
import MDAnalysis as mda
from sklearn.neighbors import NearestNeighbors
from dask import delayed, compute
from dask.distributed import Client

def compute_covariance_matrix(points):
    centered_points = points - np.mean(points, axis=0)
    return np.dot(centered_points.T, centered_points)

def curvature_filter(inside_mask, idr_atoms, ref_positions, center ,max_ref_distance):
    filtered_mask = inside_mask.copy()
    inside_indices = np.where(inside_mask)[0]
    
    
    ref_cov_matrix = compute_covariance_matrix(ref_positions)
    eigenvalues, eigenvectors = np.linalg.eig(ref_cov_matrix)
    
    # Primary curvature direction
    primary_curvature_vector = eigenvectors[:, np.argmax(eigenvalues)]
    
   
    for idx in inside_indices:
        res_pos = idr_atoms[idx].position
        
        # Vector from residue to center of ref.
        res_vector = res_pos - center
        
        # Project atom vector onto primary curvature direction
        projection = np.dot(res_vector, primary_curvature_vector)
        
        # Check if atom is too far from the primary curvature line
        distance_from_curve = np.linalg.norm(
            res_vector - projection * primary_curvature_vector
        )
        
        # If atom is too far from the primary curvature, remove it
        if distance_from_curve > max_ref_distance * 0.8:
            filtered_mask[idx] = False
    
    return filtered_mask


def check_termini(inside_mask, idr_atoms, ref_positions, center):
    rescued_mask = inside_mask.copy()
    
    ref_distances = np.linalg.norm(ref_positions - center, axis=1) #TODO: update this for efficency? 
    tight_cutoff = np.max(ref_distances) * 0.5 # Tighter than original alpha filter

    # check that termini are within tighter cutoff
    terminiAtoms = [144, 145, 223, 224, 225, 226] 
    for idx, atom in enumerate(idr_atoms):
        if not inside_mask[idx] and (atom.resid in terminiAtoms):
    
            dist_to_center = np.linalg.norm(atom.position - center)

            # Re-include if within tighter boundary
            if dist_to_center < tight_cutoff:
                rescued_mask[idx] = True

    resids = idr_atoms.resids
    resid_to_idx = {resid: idx for idx, resid in enumerate(resids)}

    # verify that termini atoms are marked as in 
    # if residue near termini is in 
    if 143 in resids and rescued_mask[resid_to_idx[143]]:
        for r in [144, 145]:
            if r in resid_to_idx:
                rescued_mask[resid_to_idx[r]] = True

    if 226 in resids and rescued_mask[resid_to_idx[226]]:
        for r in [223, 224, 225]:
            if r in resid_to_idx:
                rescued_mask[resid_to_idx[r]] = True

    return rescued_mask


def exterior_filter(inside_mask, exterior_ids, idr_atoms, max_ref_distance):

    filtered_mask = inside_mask.copy()
    indices_marked_inside = np.where(inside_mask)[0]

    # get exterior residues and their positions 
    exterior_atoms = u.select_atoms(f"resid {' '.join(map(str, exterior_ids))}")
    exterior_positions = exterior_atoms.positions

    # filter to get closest distance 
    ext_nbrs = NearestNeighbors(n_neighbors=1, algorithm='ball_tree').fit(exterior_positions)
    idr_positions = idr_atoms.positions
    ext_distances, _ = ext_nbrs.kneighbors(idr_positions)
    ext_distances = ext_distances.flatten()

    exterior_threshold = max_ref_distance * 0.7

    for i, idx in enumerate(indices_marked_inside):
        # If atom is too close to an exterior atom, mark it as outside
        if ext_distances[idx] < exterior_threshold:
            filtered_mask[idx] = False

    return filtered_mask


def atoms_in_volume(universe, frame_index, reference_atom_ids, exterior_ids, idr_indices, alpha):

    universe.trajectory[frame_index]

    # Select inner STI1 atoms (aka reference atoms)
    ref_atoms= universe.select_atoms(f"resid {' '.join(map(str, reference_atom_ids))}")
    ref_positions = ref_atoms.positions

    # Get IDR residues and their positions 
    idr_atoms = universe.select_atoms(f"resid {idr_indices[0]}:{idr_indices[1]}") 
    idr_positions = idr_atoms.positions
    
    # Find the center of reference atoms
    center = np.mean(ref_positions, axis=0)
    
    # Get distance from avg center position to each reference atom
    ref_distances = np.linalg.norm(ref_positions - center, axis=1)
    max_ref_distance = np.max(ref_distances)

    # Get distance between center and each IDR atom (in current section given!)
    idr_vectors = idr_positions - center
    center_distances_idr = np.linalg.norm(idr_vectors, axis=1)
    
    # Filter IDR residues using scikit learn nearest neighbors 
    # (makes calculation faster so we don't need to compare everything, just those we know are close)
    nbrs = NearestNeighbors(n_neighbors=1, algorithm='ball_tree').fit(ref_positions)
    distances, _ = nbrs.kneighbors(idr_positions)
    distances = distances.flatten()
    
    # Atoms are inside if:
    # 1. They are closer to center than the furthest reference atom
    # 2. They are not too far from any reference atom 
    inside_mask = (
    (center_distances_idr < max_ref_distance * alpha) & 
    (distances  < max_ref_distance * 0.6) 
    )   

    ## extra checks 
    inside_mask = exterior_filter(inside_mask, exterior_ids, idr_atoms, max_ref_distance)
    inside_mask = curvature_filter(inside_mask, idr_atoms, ref_positions, center, max_ref_distance)
    inside_mask = check_termini(inside_mask, idr_atoms, ref_positions, center)

    ### Return the atoms inside the volume for current snapshot
    atoms_inside = idr_atoms[inside_mask]
    count_trial = np.zeros(len(np.arange(idr1_range[0], idr1_range[1]+1, 1)))

    for i in atoms_inside: 
        index = i.resid - idr1_range[0]
        count_trial[index] += 1

    return count_trial

trials = list(range(1, 11))
ref_ids = [146, 147, 150, 153, 154, 159, 163, 166, 167, 175, 176, 179, 180, 181, 182, 183,186, 193, 194, 199, 202, 203, 209, 212,213, 216] ##Atoms that make up inner part of STI1 -- acounts for starting at 0
exterior_ids = [148, 149, 151, 152, 155, 156, 157, 158, 160, 161, 162, 164, 165, 168, 169, 170, 171, 172, 173, 174, 177, 178, 184, 185, 187, 188, 189, 190, 191, 192, 195, 196, 197, 198, 200, 201, 204, 205, 206, 207, 208, 210, 211, 214, 215, 217, 218, 219, 220, 221, 222]
idr1_range = [75,145]  # or [223,325] for second half
for trial in trials: 
  
    path = f"../Data/WT_bound/CALVADOS3COM_2.0_MD_gpu_trial{trial}_Dsk2_full_bound/Dsk2_full_bound/0/"
    u = mda.Universe(f"{path}/Dsk2_full_bound.pdb", f"{path}/Dsk2_full_bound.dcd")
    
    count_trial = np.zeros(len(np.arange(idr1_range[0], idr1_range[1]+1, 1)))

    numFrames = 0

    job_list = []
    results = np.zeros(len(np.arange(idr1_range[0], idr1_range[1]+1, 1))) #[]
    for frame_index in range(u.trajectory.n_frames):
        
        #job_list.append(delayed(atoms_in_volume)(u, frame_index, ref_ids, idr1_range, alpha=2.0)) # parallel
        job_list.append(delayed(atoms_in_volume(u, frame_index, ref_ids, exterior_ids, idr1_range, alpha=2.0))) # serial

        if (frame_index + 1) % 1000 == 0:
            batch_number = (frame_index + 1) / 1000
            print(f'processing batch {batch_number}...')
            result = compute(*job_list)

            job_list = []
            for array in result:
                results += np.array(array) # Add's each snapshot to a total array that represents a single run

    if len(job_list) != 0:
        print(f'Processing remaining tasks (n-tasks {len(job_list)})')
        result = compute(*job_list)
        for array in result:
            results += np.array(array) # Add's each snapshot to a total array that represents a single run
        
    np.save(f"Histos/Dsk2_full_bound_trial{trial}_innerVol_idr{idr1_range[0]}_{idr1_range[1]}.npy", results)

    
    
