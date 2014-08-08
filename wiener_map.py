import numpy as np
import h5py
import healpy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import ndimage
from scipy import signal


with h5py.File('yamldriver/timestream/map_full.hdf5', 'r') as f:
    map_data = f['map'][...]
# inmap = map_data[0, 0, :].reshape(64, -1)
inmap = map_data[0, 0, :][:, np.newaxis]
print inmap.shape
wiener_map = signal.wiener(inmap).reshape(-1)
fig = plt.figure(1, figsize=(13,5))
healpy.mollview(wiener_map, fig=1, title='')
healpy.graticule()
fig.savefig('wiener_map_full.png')
fig.clf()

