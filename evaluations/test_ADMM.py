import sys
sys.path.append("../csbw_20_project")

import torch
import numpy as np
import cvxpy as cp
import time

from datasets.protein_dataset import ProteinDataset
from models.ADMM import ADMM
from models.RMSD import RMSD

data_dir="data/cmap_coord_pairs/"
input_file="data/good_fragment_ids.txt"
predicted_fragments_dir = "outputs/predicted_fragments/"
result_log_file="outputs/results/result_on_test_fragments_ADMM.csv"

pd = ProteinDataset(data_dir=data_dir, file=input_file)
data_loader = torch.utils.data.DataLoader(pd, batch_size=1, shuffle=False)
ADMM = ADMM()
RMSD = RMSD()

for i, data  in enumerate(data_loader):
    print("Evaluating {}th fragment ... ...".format(i+1))
    
    full_dist_map = data[0].squeeze(dim=0).numpy()
    full_native_3d_coords = data[1].squeeze(dim=0).numpy()
    CA_dist_map = full_dist_map[0:256:4, 0:256:4]
    CA_coords = full_native_3d_coords[0:256:4]
    print(CA_dist_map.shape, CA_coords.shape)
    
    start_time = time.time()
    predicted_3d_coords = ADMM.solve_sdp(CA_dist_map)
    bio_rmsd = RMSD.compute_by_SVD(CA_coords, predicted_3d_coords)
    run_time = (time.time()-start_time)/60
    
    print("{:.3f} {:.3f}\n".format(bio_rmsd, run_time))
    
    with open(result_log_file, "a") as result_log_handle:
        result_log_handle.write("{:.3f} {:.3f}\n".format(bio_rmsd, run_time))
        
    np.savetxt(predicted_fragments_dir+pd.get_record_id(i)+".txt", predicted_3d_coords, fmt="%.3f", delimiter=" ", newline=" ")
