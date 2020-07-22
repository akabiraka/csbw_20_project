import sys
sys.path.append('../csbw_20_project')
import numpy as np
from Bio.PDB import *

import pyrosetta
import rosetta
pyrosetta.init()

import configs as CONFIGS
from utils.clean_slate import CleanSlate
from datasets.distance_matrix import DistanceMatrix
from datasets.coordinates import Coordinates
import utils.data_utils as DataUtils

def get_pdb_id_chain_id(line):
    """
    Given a line where first item is a pdb_id with chain_id like '1A62A',
    this method returns the pdb_id and chain_id. If ther are multiple chain_ids 
    like '1A62ACD", it only returns first chain id like 'A'.
    """
    line = line.split()
    pdb_id = line[0][:4].lower()
    chain_id = line[0][4]
    return pdb_id, chain_id

def get_bb_atoms_coordinates(pose):
    coordinates = []
    for i in range(1, pose.total_residue()+1):
        N_xyz = pose.residue(i).xyz("N")
        coordinates.append(np.array(N_xyz))
        CA_xyz = pose.residue(i).xyz("CA")
        coordinates.append(np.array(CA_xyz))
        O_xyz = pose.residue(i).xyz("O")
        coordinates.append(np.array(O_xyz))
        if "GLY" not in pose.residue(i).name():
            CB_xyz = pose.residue(i).xyz("CB")
            coordinates.append(np.array(CB_xyz))
    
    coordinates = np.array(coordinates)
    return coordinates

def check_RMSD_of_centroid_translation():
    ground_truth_pose = pyrosetta.pose_from_file("PyRosetta/data/fragments/4yycA_0.pdb")
    translated_pose = pyrosetta.pose_from_file("PyRosetta/data/translated_fragments/4yycA_0.pdb")
    ground_truth_pose_coordinates = get_bb_atoms_coordinates(ground_truth_pose)
    translated_pose_coordinates = get_bb_atoms_coordinates(translated_pose)
    from models.rmsd_loss_quaternion import RMSDLossQuaternion
    from models.rmsd_loss_biopython import RMSDLossBiopython
    from models.rmsd_loss_tpl import RMSDLossTPL
    import torch
    rmsdLoss = RMSDLossQuaternion()
    loss = rmsdLoss(torch.tensor(ground_truth_pose_coordinates).unsqueeze(dim=0), torch.tensor(translated_pose_coordinates, requires_grad=True).unsqueeze(dim=0))
    print(loss)
    rmsdLoss = RMSDLossBiopython()
    loss = rmsdLoss(torch.tensor(ground_truth_pose_coordinates).unsqueeze(dim=0), torch.tensor(translated_pose_coordinates, requires_grad=True).unsqueeze(dim=0))
    print(loss)
    rmsdLoss = RMSDLossTPL(device='cpu')
    loss = rmsdLoss(torch.tensor(ground_truth_pose_coordinates).unsqueeze(dim=0), torch.tensor(translated_pose_coordinates, requires_grad=True).unsqueeze(dim=0))
    print(loss)
# check_RMSD_of_centroid_translation()

def translate_to_centroid(pose):
    coordinates = get_bb_atoms_coordinates(pose)
    centroid = coordinates.mean(axis=0)
    pose.translate(-centroid)


# fragments_dir = "PyRosetta/data/test_set_fragments/fragments/"
# fragment_ids_file = "PyRosetta/data/test_set_good_fragment_ids.txt"
# fragments_dir = "PyRosetta/data/val_set_fragments/fragments/"
# fragment_ids_file = "PyRosetta/data/val_set_good_fragment_ids.txt"
# fragments_dir = "PyRosetta/data/train_set_fragments/"
# fragment_ids_file = "PyRosetta/data/train_set_good_fragment_ids.txt"
fragments_dir = "PyRosetta/data/fragments/"
fragment_ids_file = "PyRosetta/data/fragment_ids.txt"
translated_fragments_dir = "PyRosetta/data/translated_fragments/"
translated_fragments_cmap_coords_dir = "PyRosetta/data/translated_fragments_cmap_coord_pairs/"

distanceMatrix = DistanceMatrix(output_matrix_type="4N4N", input_dir=translated_fragments_dir, input_file_format=CONFIGS.DOT_PDB, parser=PDBParser(QUIET=True), save=False)
coordinates = Coordinates(input_dir=translated_fragments_dir, input_file_format=CONFIGS.DOT_PDB, parser=PDBParser(QUIET=True), save=False, atoms=CONFIGS.BACKBONE_ATOMS)

file_content = open(fragment_ids_file, "r")
for i, fragment_id in enumerate(file_content):
    print("Processing {}th fragment ... ...".format(i+1))
    fragment_id = fragment_id.rstrip()
    pdb_id, chain_id = get_pdb_id_chain_id(fragment_id)
    
    print("Translating {}".format(fragment_id))
    pose = pyrosetta.pose_from_file("{}{}.pdb".format(fragments_dir, fragment_id))
    translate_to_centroid(pose)
    pose.dump_pdb("{}{}.pdb".format(translated_fragments_dir, fragment_id))
    
    dist_matrix = distanceMatrix.generate(filename=fragment_id, pdb_id=pdb_id, chain_id=chain_id)
    d3_coords = coordinates.generate(filename=fragment_id, pdb_id=pdb_id, chain_id=chain_id)
    DataUtils.save_using_pickle([dist_matrix, d3_coords], "{}{}.pkl".format(translated_fragments_cmap_coords_dir, fragment_id))
    