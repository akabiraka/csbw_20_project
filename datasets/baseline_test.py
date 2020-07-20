import sys
sys.path.append('../csbw_20_project')
import traceback
from Bio.PDB import *

import configs as CONFIGS
from utils.clean_slate import CleanSlate
from datasets.my_pdb_data import MyPDBData, CAAtomSelector
from datasets.distance_matrix import DistanceMatrix
from datasets.coordinates import Coordinates
import utils.data_utils as DataUtils

def get_pdb_id(line):
        """
        Given a line where first item is a pdb_id with chain_id like '1A62A',
        this method returns the pdb_id and chain_id. If ther are multiple chain_ids 
        like '1A62ACD", it only returns first chain id like 'A'.
        """
        line = line.split()
        pdb_id = line[0][:4].lower()
        chain_id = line[0][4:]
        return pdb_id, chain_id
    
myPDBData = MyPDBData()
cln = CleanSlate()
distanceMatrix = DistanceMatrix(output_matrix_type="NN", input_dir=CONFIGS.CLEAN_PDB_DIR, input_file_format=CONFIGS.DOT_PDB, parser=PDBParser(QUIET=True), save=False, atom_1="CA", atom_2="CA")
coordinates = Coordinates(input_dir=CONFIGS.CLEAN_PDB_DIR, input_file_format=CONFIGS.DOT_PDB, parser=PDBParser(QUIET=True), save=False, atoms=["CA"])

good_pdb_ids_handle = open("data/good_pdb_ids_handle.txt", "a")    
file_content = open("data/test.txt", "r")

n_proteins_to_skip = 0
total_proteins_to_evaluate = n_proteins_to_skip + 5

for i, line in enumerate(file_content):
    if i < n_proteins_to_skip: continue
    print("Processing {}th protein ... ...:".format(i+1))
    
    pdb_id, chain_id = get_pdb_id(line)
    if len(chain_id) > 1: continue
    
    myPDBData.download(pdb_id)
    cleaned_structure = myPDBData.clean(pdb_id, chain_id, select=CAAtomSelector(chain_id))
    chain = myPDBData.get_chain_from_structure(cleaned_structure, chain_id)
    
    ca_dist_matrix = distanceMatrix.generate(filename=pdb_id+chain_id, pdb_id=pdb_id, chain_id=chain_id)
    ca_coordinates = coordinates.generate(filename=pdb_id+chain_id, pdb_id=pdb_id, chain_id=chain_id)
    DataUtils.save_using_pickle([ca_dist_matrix, ca_coordinates], CONFIGS.CONTACT_MAP_VS_COORDINATES_DIR+pdb_id+chain_id+CONFIGS.DOT_PKL)
    print(ca_dist_matrix.shape, ca_coordinates.shape)
    
    good_pdb_ids_handle.write(pdb_id+chain_id+"\n")
    if i+1 == total_proteins_to_evaluate: break
        
    cln.clean_all_files(mydir=CONFIGS.PDB_DIR, ext=CONFIGS.DOT_CIF)
    cln.clean_all_files(mydir=CONFIGS.CLEAN_PDB_DIR, ext=CONFIGS.DOT_PDB)    
    print()