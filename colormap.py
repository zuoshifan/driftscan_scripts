"""
User defined color-maps.
"""

import matplotlib as mpl

# rjet comes from some modification of 'jet' colormap
rjet_dict = {'red': ((0., 0, 0),
                     (0.35, 0, 0),
                     (0.73, 0.7, 0.7),
                     (0.8,1, 1),
                     (1, 1, 1)),
           'green': ((0., 0, 0),
                     (0.125,0, 0),
                     (0.45,1, 1),
                     (0.84,1, 1),
                     (0.91,0,0),
                     (1, 0, 0)),
            'blue': ((0., 0.5, 0.5),
                     (0.11, 1, 1),
                     (0.34, 1, 1),
                     (0.65,0, 0),
                     (1, 0, 0))}

# jet09: the scale of 'jet' cmap dived by 0.89, so the most red is at the top of the colorbar
jet09_dict = {'red': ((0., 0, 0),
                     (0.393, 0, 0),
                     (0.7416, 1, 1),
                     (1,1, 1)),
                     # (1, 0.5, 0.5)),
           'green': ((0., 0, 0),
                     (0.14,0, 0),
                     (0.421,1, 1),
                     (0.719,1, 1),
                     (1,0,0)),
                     # (1, 0, 0)),
            'blue': ((0., 0.5, 0.5),
                     (0.1236, 1, 1),
                     (0.382, 1, 1),
                     (0.73,0, 0),
                     (1, 0, 0))}


rjet = mpl.colors.LinearSegmentedColormap('rjet',rjet_dict,256)
jet09 = mpl.colors.LinearSegmentedColormap('jet09',jet09_dict,256)

if __name__ == '__main__':
    import numpy as np
    import matplotlib.pyplot as plt
    
    plt.pcolor(np.random.rand(10,10),cmap=jet09)
    plt.colorbar()
    plt.show()