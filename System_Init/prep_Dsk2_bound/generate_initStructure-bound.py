from modeller import *
from modeller import restraints
from modeller.automodel import *
import shutil


def extract_sequence(mdl, chain_id):
    
    # convert from three-letter to one-letter amino acid sequence 
    aa_dict = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G', 'HIS': 'H',
        'ILE': 'I', 'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q',
        'ARG': 'R', 'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
    }

    sequence = ""
    for chain in mdl.chains:
        if chain.name == chain_id:
            for residue in chain.residues:

                sequence += aa_dict.get(residue.name, 'X')  # Default to 'X' for unknown residues
    return sequence


pdb_file = "Dsk2_full_bound_relaxed_structure.pdb" 
chain_id = "A"           
numTrials = 10

class MyModel(allhmodel):
        def select_atoms(self):
            return selection(self)
env = Environ()
mdl = Model(env)
struc = mdl.read(file= pdb_file)


seq = extract_sequence(mdl, chain_id)
seq_length=str(len(seq))

with open("Dsk2_full_bound-seq.txt" , "w") as f:

    f.write("%s" %(seq))
    f.close()


    
pdb_code = pdb_file[:-4]

#### Create the alignment file ####
alignment_file = "alignment.ali"
first_residue = " A" 
chain_id = " 1"
with open(alignment_file, "w") as f:
        f.write(f">P1;{pdb_code}\n")
        f.write(f"structure:{pdb_file}:{chain_id}:{first_residue}:{seq_length}:::::\n")
        f.write(seq + "*\n\n")
        f.write(f">P1;target\n")
        f.write(f"sequence:target:{chain_id}:{first_residue}:{seq_length}:::::\n")
        f.write(seq + "*\n")

with open(alignment_file, "r") as f:
    print(f.read())


#### Generate initial structures ####
env.edat.nonbonded_sel_atoms = 1 

a = MyModel(
        env,
        alnfile="alignment.ali",
        knowns=pdb_code,
        sequence="target"
    )

a.starting_model = 1
a.ending_model = numTrials

a.make()







