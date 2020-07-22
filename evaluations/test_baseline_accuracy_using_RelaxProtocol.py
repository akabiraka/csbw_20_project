import sys
sys.path.append("../csbw_20_project")    

import pyrosetta
import rosetta
pyrosetta.init()

from models.relax_protocol import RelaxProtocol
import utils.data_utils as DataUtils

relaxProtocol = RelaxProtocol()

fragments_dir = "data/test_set_fragments/fragments/"
cmap_coord_pairs_dir = "data/test_set_cmap_coord_pairs/c_map_vs_coord_pairs/"
best_poses_by_energy_score_dir = "outputs/best_poses_RelaxProtocol_full_atom_distance/"
rmsd_scores_log_file = "outputs/results/rmsd_on_test_fragments_RelaxProtocol.csv"
file_content = open("data/test_set_good_fragment_ids.txt")

n_fragments_to_skip = 0
total_fragments_to_evaluate = n_fragments_to_skip + 20

for i, fragment_id in enumerate(file_content):
    if i < n_fragments_to_skip: continue
    
    fragment_id = fragment_id.rstrip()
    print("Processing {}th fragment:{} ... ...:".format(i+1, fragment_id))
    native_pose_pdb = pyrosetta.pose_from_file(fragments_dir+fragment_id+".pdb")
    sequence = native_pose_pdb.sequence()
    distance_matrix = DataUtils.load_using_pickle(cmap_coord_pairs_dir+fragment_id+".pkl")[0]
    # distance_matrix[distance_matrix <= 8.0] = 6.0
    # distance_matrix[distance_matrix > 8.0] = None
    print(distance_matrix.shape, len(sequence))
    
    energy_scores, best_pose, run_time = relaxProtocol.run(sequence, distance_matrix, dist_map_type="4N4N", atom_1="CA", atom_2="CA")
    best_pose.dump_pdb("{}{}.pdb".format(best_poses_by_energy_score_dir, fragment_id)) 
    rmsd_score = pyrosetta.rosetta.core.scoring.all_atom_rmsd(native_pose_pdb, best_pose) 
    
    print("{:.3f} {:.3f}\n".format(rmsd_score, run_time))
    with open(rmsd_scores_log_file, "a") as rmsd_scores_log_file_handle:
        rmsd_scores_log_file_handle.write("{:.3f} {:.3f}\n".format(rmsd_score, run_time))
        
    if i+1 == total_fragments_to_evaluate: break