#convertrxmaps.py
#loads reaction maps from .mat file

import scipy.io

mat = scipy.io.loadmat('reactionmaps.mat')
rmaps = mat['reactionmaps']

Fmap = rmaps[0][0][0]
Cmap = rmaps[0][0][1]
Gmap = rmaps[0][0][2]

Fmap2 = rmaps[0][0][3]
Cmap2 = rmaps[0][0][4]
Gmap2 = rmaps[0][0][5]

