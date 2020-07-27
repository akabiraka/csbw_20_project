import sys
sys.path.append("../csbw_20_project")

import torch
from datasets.protein_dataset import ProteinDataset

data_dir = "/media/akabir/New Volume1/2nd_time/train_set_translated_fragments_cmap_coord_pairs/"
fragment_ids_file = "/media/akabir/New Volume1/2nd_time/train_set_good_fragment_ids.txt"
problematic_fragment_ids_file = "outputs/results/problematic_fragment_ids.txt"

proteinDataset = ProteinDataset(data_dir=data_dir, file=fragment_ids_file)
data_loader = torch.utils.data.DataLoader(proteinDataset, batch_size=1, shuffle=False)

for i, data  in enumerate(data_loader):
    fragment_id = proteinDataset.get_record_id(i)
    print("Processing {}:{}".format(i+1, fragment_id))
    full_dist_map = data[0].squeeze(dim=0).numpy()
    full_native_3d_coords = data[1].squeeze(dim=0).numpy()
    # print(full_dist_map.shape, full_native_3d_coords.shape)
    rows, cols = full_dist_map.shape
    if rows!=256 or cols!=256:
        print("Culprit found: {}:{}:{}".format(i+1, fragment_id, full_dist_map.shape))
        with open(problematic_fragment_ids_file, "a") as problematic_fragment_ids_file_handle:
            problematic_fragment_ids_file_handle.write("{} {}\n".format(fragment_id, full_dist_map.shape))
    
    # if i+1==5:
    #     break