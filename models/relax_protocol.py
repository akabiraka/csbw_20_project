import sys
sys.path.append("../csbw_20_project")

import time
import numpy as np

import pyrosetta
import rosetta
pyrosetta.init()


import utils.data_utils as DataUtils

# http://faculty.washington.edu/dimaio/files/HW_2_2020.pdf
#                          GLY         ASP/ASN     ILE/VAL     PRO         OTHERS
ANGLE_DISTRIBUTIONS = [[[87, 7], [-140, 165], [-132, 153], [-64, 145], [-136, 153]],
                        [[-66, -35], [-78, 141], [-86, 127], [-60, -29], [-76, 143]],
                        [[70, -144], [-108, 103], [-118, 125], [-60, -29], [-112, 119]],
                        [[105, 170], [-97, 5], [-91, -9], [-77, 161], [-91, -9]],
                        [[-171, 177], [-64, -39], [-63, -42], [-77, 161], [-63, -42]],
                        [[-87, 163], [57, 39], [57, 39], [-84, -2], [57, 39]]]
BACKBONE_ATOMS = ['CA', 'CB', 'N', 'O']

class RelaxProtocol(object):
    def __init__(self):
        super(RelaxProtocol, self).__init__()
        
    def init_ramachandran_angles(self, pose, angles):
        n_residues = pose.total_residue()
        for i in range(1, n_residues+1):
            ith_residue = pose.residue(i).name()
            if "GLY" in ith_residue:
                pose.set_phi(i, angles[0][0])
                pose.set_psi(i, angles[0][1])
            elif "ASP" in ith_residue or "ASN" in ith_residue:
                pose.set_phi(i, angles[1][0])
                pose.set_psi(i, angles[1][1])
            elif "ILE" in ith_residue or "VAL" in ith_residue:
                pose.set_phi(i, angles[2][0])
                pose.set_psi(i, angles[2][1])
            elif "PRO" in ith_residue:
                pose.set_phi(i, angles[3][0])
                pose.set_psi(i, angles[3][1])
            else:
                pose.set_phi(i, angles[4][0])
                pose.set_psi(i, angles[4][1])
                
        return pose
    
    def get_atom_pair_ids(self, pose, i=0, ith_atom_idx="CA", j=0, jth_atom_idx="CA"):
        res_name_1 = pose.residue(i).name()
        res_name_2 = pose.residue(j).name()
        
        if ith_atom_idx=="CB" and 'GLY' in res_name_1:
            ith_atom_idx = "CA"
        
        if jth_atom_idx=="CB" and "GLY" in res_name_2 :
            jth_atom_idx = "CA"
        
        id_i = pyrosetta.AtomID(pose.residue(i).atom_index(ith_atom_idx), i)
        id_j = pyrosetta.AtomID(pose.residue(j).atom_index(jth_atom_idx), j)
        # print(id_i, id_j)
        return id_i, id_j
        
    def add_constraint_from_4N4N_distance_map(self, pose, dist_mat):
        """[Add distance constraints to each residue]

        Args:
            pose (pyrosetta pose): 
            dist_mat (2d matrix): (4n, 4n) distance matrix where n is the length of the sequence in pose
        """
        # dist_mat: (4n, 4n)
        n_residues = pose.total_residue()
        for i in range(1, n_residues+1):
            for j in range(1, n_residues+1):
                for k, atom_1 in enumerate(BACKBONE_ATOMS):
                    for l, atom_2 in enumerate(BACKBONE_ATOMS):
                        distance_value = dist_mat[4*(i-1)+k][4*(j-1)+l]
                        if np.isnan(distance_value):
                            continue
                        distance_lower_bound = distance_value-1 if distance_value-1 > 0 else 0
                        distance_upper_bound = distance_value+1
                        id_i, id_j = self.get_atom_pair_ids(pose, i, atom_1, j, atom_2)
                        ij_func = rosetta.core.scoring.constraints.BoundFunc(distance_lower_bound, distance_upper_bound, 1.0, "cst1")
                        ij_cnstrnt = rosetta.core.scoring.constraints.AtomPairConstraint(id_i, id_j, ij_func)
                        pose.add_constraint(ij_cnstrnt)
        return pose
    
    def add_constraint_from_NN_distance_map(self, pose, dist_mat, atom_1="CA", atom_2="CA"):
        n_residues = pose.total_residue()
        for i in range(1, n_residues+1):
            for j in range(1, n_residues+1):
                distance_value = dist_mat[i][j]
                if np.isnan(distance_value):
                    continue
                distance_lower_bound = distance_value-1 if distance_value-1 > 0 else 0
                distance_upper_bound = distance_value+1
                id_i, id_j = self.get_atom_pair_ids(pose, i, atom_1, j, atom_2)
                ij_func = rosetta.core.scoring.constraints.BoundFunc(distance_lower_bound, distance_upper_bound, 1.0, "cst1")
                ij_cnstrnt = rosetta.core.scoring.constraints.AtomPairConstraint(id_i, id_j, ij_func)
                pose.add_constraint(ij_cnstrnt)
        return pose
            
    def run(self, primary_seq, dist_mat, dist_map_type="4N4N", atom_1="CA", atom_2="CA"):
        """[summary]

        Args:
            primary_seq ([type]): [description]
            distance_matrix ([type]): [description]
            dist_map_type (str, optional): 4N4N, NN. Defaults to "4N4N".
                If NN, atom_1 and atom_2 will be used.

        Returns:
            [type]: [description]
        """
        scorefxn_scores = []
        best_pose = None
        best_energy_score = float('inf') #np.inf
        
        start_time = time.time()
        for i, angle in enumerate(ANGLE_DISTRIBUTIONS):
            print("Running relax protocol using {}th angle distribution".format(i+1))
            pose = pyrosetta.pose_from_sequence(primary_seq)
            print("dist-mat: {}, residues: {}".format(dist_mat.shape, pose.total_residue()))
            
            # initializing pose with ramachandran angle distribution
            pose = self.init_ramachandran_angles(pose, angle)

            # adding distance constrains for all pair of residues
            if dist_map_type=="4N4N":
                pose = self.add_constraint_from_4N4N_distance_map(pose, dist_mat)
            else:
                pose = self.add_constraint_from_NN_distance_map(pose, dist_mat, atom_1="CA", atom_2="CA")
            
            # defining scoring function        
            scorefxn = pyrosetta.get_fa_scorefxn()
            scorefxn.set_weight(rosetta.core.scoring.atom_pair_constraint, 1.0)

            # defining FastRelax protocol
            print("estimating pose using relax protocol ... ...")
            fastrelax = rosetta.protocols.relax.FastRelax(scorefxn, 10)
            fastrelax.apply(pose)
        
            energy_score = scorefxn(pose)
            scorefxn_scores.append(energy_score)
            if energy_score < best_energy_score:
                best_energy_score = energy_score
                best_pose = pose
                
            # if i==0: break
        
        run_time = (time.time() - start_time)/60
        return scorefxn_scores, best_pose, run_time
