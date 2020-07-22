import sys
sys.path.append("../csbw_20_project")

import cvxpy as cp
import numpy as np
from sklearn.decomposition import TruncatedSVD

class ADMM(object):
    def __init__(self):
        super(ADMM, self).__init__()
        
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
        # print(len(constraints), len(objective))
        prob = cp.Problem(cp.Minimize(np.sum(objective)), constraints)     
        prob.solve()   
        # print(X.value.shape)
        
        svd = TruncatedSVD(n_components=3, n_iter=7, random_state=42)
        svd.fit(X.value)
        res = svd.transform(X.value)
        return res