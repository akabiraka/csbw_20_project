import sys
sys.path.append("../csbw_20_project")

from Bio.SVDSuperimposer import SVDSuperimposer

class RMSD(object):
    def __init__(self):
        super(RMSD, self).__init__()
        self.superimposer = SVDSuperimposer()
        
    def compute_by_biopython(self, native_3d_coords, predicted_3d_coords):
        self.superimposer.set(native_3d_coords, predicted_3d_coords)
        self.superimposer.run()
        return self.superimposer.get_rms()
    
