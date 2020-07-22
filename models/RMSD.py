import sys
sys.path.append("../csbw_20_project")

from Bio.SVDSuperimposer import SVDSuperimposer
from Bio.PDB.QCPSuperimposer import QCPSuperimposer

class RMSD(object):
    def __init__(self):
        super(RMSD, self).__init__()
        self.mySVDSuperimposer = SVDSuperimposer()
        self.myQCPSuperimposer = QCPSuperimposer()
        
    def compute_by_SVD(self, native_3d_coords, predicted_3d_coords):
        """[summary]

        Args:
            native_3d_coords: an NxDIM array
            predicted_3d_coords: an NxDIM array

        Returns:
            RMSD: float
        """
        self.mySVDSuperimposer.set(native_3d_coords, predicted_3d_coords)
        self.mySVDSuperimposer.run()
        return self.mySVDSuperimposer.get_rms()
    
    def compute_by_QCP(self, native_3d_coords, predicted_3d_coords):
        """[summary]
        QCP: Quaternion Characteristic Polynomial
        Args:
            native_3d_coords: an NxDIM array
            predicted_3d_coords: an NxDIM array

        Returns:
            RMSD: float
        """
        self.myQCPSuperimposer.set(native_3d_coords, predicted_3d_coords)
        self.myQCPSuperimposer.run()
        return self.myQCPSuperimposer.get_rms()
    

# TESTING
# import numpy as np
# native_coords = np.random.uniform(low=-30, high=30, size=(5, 3))
# predicted_coords = np.random.uniform(low=-30, high=30, size=(5, 3))
# RMSD = RMSD()
# rmsd_1 = RMSD.compute_by_SVD(native_coords, predicted_coords)
# rmsd_2 = RMSD.compute_by_QCP(native_coords, predicted_coords)
# print(rmsd_1, rmsd_2)