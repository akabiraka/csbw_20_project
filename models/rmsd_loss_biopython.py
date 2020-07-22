import sys
sys.path.append('../FDFAPBG')

from Bio.SVDSuperimposer import SVDSuperimposer
import torch
import torch.nn as nn

class RMSDLossBiopython(nn.Module):
    def __init__(self):
        super(RMSDLossBiopython, self).__init__()
        self.superImposer = SVDSuperimposer()
        
    def forward(self, y_prime, y):
        """[summary]

        Args:
            y_prime ((batch_size, n_atoms, 3) dimensional tensor): Predicted
            y ((batch_size, n_atoms, 3) dimensional tensor): Ground truth
        """
        batch_size, n_atoms, d = y_prime.shape
        total_loss = 0.0
        for batch in range(batch_size):
            _y_prime = y_prime[batch].detach().numpy()
            _y = y[batch].detach().numpy()
            self.superImposer.set(_y, _y_prime)
            self.superImposer.run()
            total_loss += self.superImposer.get_rms()
        
        # print(total_loss, total_loss/batch_size)
        return total_loss/batch_size
    
# y_prime = torch.randn(size=(30, 256, 3), dtype=torch.float32, requires_grad=True)
# y = torch.randn(size=(30, 256, 3), dtype=torch.float32)
# criterion = RMSDLoss()
# criterion(y_prime, y)
