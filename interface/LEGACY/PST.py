import numpy as np
import argparse as ap
import matplotlib.pyplot as plt
import matplotlib.tri as tri

from interface.Global import load_mesh, read_soln

def GetCentres(nodes, eles):
    centres = np.zeros((eles.NC,2))
    for i in range(eles.NC):
        centres[i,:] = np.mean( nodes[ eles.cn_con[ eles.cn_off[i]:eles.cn_off[i+1] ],:] , axis=0)
    return centres

def process_PST(args):
    fmap = {'r':0, 'u':1, 'v':2, 'p':3, 'U':4, 'M':5}

    _, nodes, eles = load_mesh(args.mesh)
    centres = GetCentres(nodes, eles)
    soln = read_soln(args.soln)
    soln = np.column_stack((soln, np.sqrt(soln[:,1]*soln[:,1]+soln[:,2]*soln[:,2])))
    soln = np.column_stack((soln, soln[:,4]/np.sqrt(1.4*soln[:,3]/soln[:,0])))

    triang = tri.Triangulation(centres[:,0], centres[:,1])
    plt.tripcolor(triang, soln[:,fmap[args.field]], cmap='magma')
    plt.colorbar()
    ax = plt.gca()
    ax.set_aspect('equal', adjustable='box')
    plt.savefig(args.visout)
