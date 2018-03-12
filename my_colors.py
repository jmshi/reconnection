import matplotlib  # in python
import matplotlib.pyplot as plt
import numpy as np

def get_jhcolors(debug=False):
  ct = np.loadtxt('jhcolors.tab')
  cm = matplotlib.colors.ListedColormap(ct/255.0)
  if debug:
    plt.imshow([[-1,1],[1,-1]], cmap=cm) # for example
    plt.colorbar(pad=0.01,fraction=0.045,orientation='vertical')
  return cm
