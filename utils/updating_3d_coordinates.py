# this file contains 2 example codes to update 3d coordinates using Biopython and pyrosetta

# using Biopython
# parser = PDBParser(QUIET=True)
# pdbio = PDBIO()

# predicted_3d_coords = torch.load(predicted_data_path+"result_0.pt", map_location=torch.device('cpu')).numpy()[0]
# structure = parser.get_structure(fragment_id, fragment_filename)
# for i, residue in enumerate(structure[0][chain_id]):
#     j = i*4
#     CA_xyz = predicted_3d_coords[j]
#     residue["CA"].set_coord(CA_xyz)
#     if "GLY" not in residue.get_resname():
#         CB_xyz = predicted_3d_coords[j+1]
#         residue["CB"].set_coord(CB_xyz)
#     N_xyz = predicted_3d_coords[j+2]
#     residue["N"].set_coord(N_xyz)
#     O_xyz = predicted_3d_coords[j+3]
#     residue["O"].set_coord(O_xyz)
    
# pdbio.set_structure(structure)
# fragment_filename = fragments_dir + fragment_id + "_onlybackboneatoms.pdb"
# pdbio.save(fragment_filename, select=AllBackboneAtomSelector(chain_id))

# using pyrosetta
# pose = pyrosetta.pose_from_file(fragment_filename)

# for i in range(1, pose.total_residue()+1):
#     j = (i-1)*4
#     xyz = predicted_3d_coords[j]
#     CA_xyz = pyrosetta.rosetta.numeric.xyzVector_double_t(xyz[0], xyz[1], xyz[2])
#     xyz = predicted_3d_coords[j+1]
#     CB_xyz = pyrosetta.rosetta.numeric.xyzVector_double_t(xyz[0], xyz[1], xyz[2])
#     xyz = predicted_3d_coords[j+2]
#     N_xyz = pyrosetta.rosetta.numeric.xyzVector_double_t(xyz[0], xyz[1], xyz[2])
#     xyz = predicted_3d_coords[j+3]
#     O_xyz = pyrosetta.rosetta.numeric.xyzVector_double_t(xyz[0], xyz[1], xyz[2])
    
#     pose.residue(i).set_xyz("CA", CA_xyz)
#     if "GLY" not in pose.residue(i).name():
#         pose.residue(i).set_xyz("CA", CB_xyz)
#     pose.residue(i).set_xyz("N", N_xyz)
#     pose.residue(i).set_xyz("O", O_xyz)
    
# pose.dump_pdb("{}{}_updated.pdb".format(fragments_dir, fragment_id))