import sys
sys.path.append("../csbw_20_project")

import torch
from sklearn.decomposition import TruncatedSVD
from datasets.protein_dataset import ProteinDataset
from Bio.SVDSuperimposer import SVDSuperimposer
import numpy as np
import cvxpy as cp
import time

class TestADMM(object):
    def __init__(self, data_dir, input_file, result_log_file, n_fragments_to_evaluate=None):
        super(TestADMM, self).__init__()
        self.superimposer = SVDSuperimposer()
        pd = ProteinDataset(data_dir=data_dir, file=input_file)
        self.data_loader = torch.utils.data.DataLoader(pd, batch_size=1, shuffle=False)
        self.result_log_file = result_log_file
        self.n_fragments_to_evaluate = n_fragments_to_evaluate
        
    def solve_sdp(self, dist_map):
        mat_len = dist_map.shape[0]
        X = cp.Variable((mat_len,mat_len))
        constraints = [X >> 0]
        # prob = cp.Problem(cp.Minimize(np.sum(([(X[i][i]+X[j][j]-2*X[i][j]-dist_map[i][j]**2)**2 for i in range(mat_len) for j in range(mat_len)]))),constraints)
        # #prob = cp.Problem(cp.Minimize(cp.trace(X)), constraints)
        # prob.solve()
        
        objective = []
        for i in range(mat_len):
            for j in range(mat_len):
                if np.isnan(dist_map[i][j]):
                    continue
                objective.append((X[i][i]+X[j][j]-2*X[i][j]-dist_map[i][j]**2)**2)
        print(len(constraints), len(objective))
        prob = cp.Problem(cp.Minimize(np.sum(objective)), constraints)     
        prob.solve()   
        print(X.value.shape)
        
        svd = TruncatedSVD(n_components=3, n_iter=7, random_state=42)
        svd.fit(X.value)
        res = svd.transform(X.value)
        return res

    def compute_RMSD(self, native_3d_coords, predicted_3d_coords):
        self.superimposer.set(native_3d_coords, predicted_3d_coords)
        self.superimposer.run()
        return self.superimposer.get_rms()
        
    def compute_base_line_test_accuracy(self, th=None, is_contact_map=False):
        """[summary]

        Args:
            th (float): Threshold. If distance is higher than th, replace that with a 0
        """
        for i, data  in enumerate(self.data_loader):
            print("Evaluating {}th fragment ... ...".format(i+1))
            
            full_dist_map = data[0].squeeze(dim=0).numpy()
            full_native_3d_coords = data[1].squeeze(dim=0).numpy()
            
            CA_dist_map = full_dist_map[0:256:4, 0:256:4]
            if is_contact_map:
                CA_dist_map[CA_dist_map < th] = 0.
                CA_dist_map[CA_dist_map > th] = 1.
            if th is not None:
                CA_dist_map[CA_dist_map > th] = 0
            CA_coords = full_native_3d_coords[0:256:4]
            print(CA_dist_map.shape, CA_coords.shape)
            
            start_time = time.time()
            predicted_3d_coords = self.solve_sdp(CA_dist_map)
            # raise
            bio_rmsd = self.compute_RMSD(CA_coords, predicted_3d_coords)
            run_time = (time.time()-start_time)/60
            
            print("{:.3f} {:.3f}\n".format(bio_rmsd, run_time))
            with open(self.result_log_file, "a") as result_log_handle:
                result_log_handle.write("{:.3f} {:.3f}\n".format(bio_rmsd, run_time))

            if self.n_fragments_to_evaluate is not None and (i+1)==self.n_fragments_to_evaluate: break
            
        return predicted_3d_coords
        
testADMM = TestADMM(data_dir="data/test_set_cmap_coord_pairs/c_map_vs_coord_pairs/",
                    input_file="data/test_set_good_fragment_ids.txt",
                    result_log_file="outputs/results/result_on_test_fragments_ADMM.csv",
                    n_fragments_to_evaluate=20)
testADMM.compute_base_line_test_accuracy(th=8, is_contact_map=False)
