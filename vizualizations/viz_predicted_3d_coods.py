__author__ = "Anowarul Kabir"
__updated__ = "2020-08-29 22:33:57"

import sys
sys.path.append("../csbw_20_project")

from Bio.PDB import *
from Bio.PDB.PDBIO import PDBIO, Select
from datasets.Selector import AllBackboneAtomSelector

import torch

predicted_coords_dir = "/media/akabir/New Volume1/2nd_time/DCNN_debug_2/DCNN_debug/"
filename = "result_0.pt"
predicted_coords_filename = predicted_coords_dir + filename

fragments_dir = "data/fragments/"
fragment_id = "2bwrA_0"
chain_id = "A"
fragment_filename = fragments_dir + fragment_id + ".pdb"

parser = PDBParser(QUIET=True)
pdbio = PDBIO()

predicted_3d_coords = torch.load(predicted_coords_filename, map_location=torch.device('cpu')).numpy()[0]

structure = parser.get_structure(fragment_id, fragment_filename)
for i, residue in enumerate(structure[0][chain_id]):
    j = i*4
    CA_xyz = predicted_3d_coords[j]
    residue["CA"].set_coord(CA_xyz)
    if "GLY" not in residue.get_resname():
        CB_xyz = predicted_3d_coords[j+1]
        residue["CB"].set_coord(CB_xyz)
    N_xyz = predicted_3d_coords[j+2]
    residue["N"].set_coord(N_xyz)
    O_xyz = predicted_3d_coords[j+3]
    residue["O"].set_coord(O_xyz)
    
pdbio.set_structure(structure)
fragment_filename = fragments_dir + fragment_id + "_predicted.pdb"
pdbio.save(fragment_filename, select=AllBackboneAtomSelector(chain_id))
