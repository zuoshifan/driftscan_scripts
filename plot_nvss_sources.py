import numpy as np
import aipy as a
import matplotlib
matplotlib.use('Agg')
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt


flux = 5.0 # Jy
# flux = 10.0 # Jy
frequency = 750 # MHz
catalog = 'nvss'
# catalog = 'wenss'
src = '%f/%f' % (flux, frequency / 1.0e3)
srclist, cutoff, catalogs = a.scripting.parse_srcs(src, catalog)
cat = a.src.get_catalog(srclist, cutoff, catalogs)
nsrc = len(cat) # number of sources in cat
ras = [ np.degrees(cat.values()[i]._ra) for i in xrange(nsrc) ]
decs = [ np.degrees(cat.values()[i]._dec) for i in xrange(nsrc) ]
jys = [ cat.values()[i].get_jys() for i in xrange(nsrc) ]

print len(ras)


# lon_0 is central longitude of projection.
m = Basemap(projection='moll', lon_0=0, celestial=True)

# draw parallels and meridians.
m.drawparallels(np.arange(-90.,120.,30.))
m.drawmeridians(np.arange(0.,360.,30.))

# plt.title("Mollweide Projection")
# convert to map projection coords.
# Note that ra, dec can be scalars, lists or numpy arrays.
xpt, ypt = m(ras, decs)
# convert back to ra/dec
rapt, decpt = m(xpt, ypt, inverse=True)
# m.plot(xpt, ypt, 'bo')  # plot a blue dot there
# m.scatter(xpt, ypt, s=jys, color='r', alpha=.9)
m.scatter(xpt, ypt, s=jys, c='#9a0200', edgecolors='k', alpha=1) # deep red
# put some text next to the dot, offset a little bit
# (the offset is in map projection coordinates)
# plt.text(xpt+100000, ypt+100000, 'Boulder (%5.1fW,%3.1fN)' % (rapt, decpt))
plt.savefig('sources.png')
plt.close()