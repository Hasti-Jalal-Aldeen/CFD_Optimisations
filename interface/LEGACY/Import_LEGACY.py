import numpy as np
from interface.Global import Elements, eletype_names, eletype_nnodes, eletype_nsides
import scipy as sp
import struct

import matplotlib.pyplot as plt

def parsePhysNames(PhsyNamesSec):
    n_physnames = int(PhsyNamesSec[0])
    physnames = []
    for i in range(1,n_physnames):
        ln = PhsyNamesSec[i].split(' ')
        tag = int(ln[1])
        name = ln[2]
        physnames.append( name.strip('"') )
    return np.array(physnames)

def parseEntities(EntitiesSec):
    n_entities = EntitiesSec[0].split(' ')
    n_entities = [int(entity) for entity in n_entities]
    npts, ncurves, nsurfs, nvols = n_entities

    if nvols == 0: dim = 2
    else: dim = 3
    if nsurfs == 0: raise Exception("No surface. ")
    
    curves = [None]
    for i in range(npts+1,npts+1+ncurves):
        line = EntitiesSec[i].split(' ')
        if int(line[7]) > 1: raise Exception("Invalid number of physical tags on curve. ")
        curves.append( int(line[8]) )
    
    return curves, dim

def parseNodes(NodesSec):
    edata = NodesSec[0].split(' ')
    edata = [int(num) for num in edata]
    _,n_nodes,*_ = edata

    nodes = []
    head = 1
    while head < len(NodesSec):
        eblock = NodesSec[head].split(' ')
        eblock = [int(num) for num in eblock]
        *_,para,n_nodesinblock = eblock
        if para == 1: raise Exception("Parametric nodes not supported. ")

        ncoords = NodesSec[head+1+n_nodesinblock : head+1+2*n_nodesinblock]

        ncoords = [ [float(ncoord.split(' ')[0]), float(ncoord.split(' ')[1])] for ncoord in ncoords]

        nodes = nodes+ncoords
        head = head+1+2*n_nodesinblock
    return np.array(nodes)

def parseElements2D(ElesSec, curves):
    idxs = [None, 1, 0]

    eldata = ElesSec[0].split(' ')
    eldata = [int(num) for num in eldata]
    _,n_eles,*_ = eldata
    
    eletypes = [[],[]]
    elennodes = [[],[]]
    nodes = [[],[]]
    phystags = [[],[]]
    phystags_n = [[],[]]
    head = 1
    while head < len(ElesSec):
        eblock = ElesSec[head].split(' ')
        eblock = [int(num) for num in eblock]
        endim,entag,eletype,n_elesinblock = eblock

        eletypestr = eletype_names[eletype]
        nnodes = eletype_nnodes[eletypestr]
        phystag = curves[entag]
        
        eles_block = [ [int(val) for val in ele.split(' ')[:-1]] for ele in ElesSec[head+1 : head+1+n_elesinblock] ]
        elnodes = np.array([[node-1 for node in ele[1:]] for ele in eles_block]).flatten()

        idx = idxs[endim]
        eletypes[idx] += n_elesinblock*[eletype]
        elennodes[idx] += n_elesinblock*[nnodes]
        nodes[idx] += [*elnodes]
        phystags[idx].append(phystag)
        phystags_n[idx].append(n_elesinblock)

        head = head+1+n_elesinblock

    eles = Elements()
    eles.dim = 2

    eles.N = n_eles
    eles.NC = len(elennodes[idxs[2]])
    eles.eletypes = np.array([*eletypes[0], *eletypes[1]], dtype=np.int32 )

    eles.cn_con = np.array( [*nodes[0] , *nodes[1]], dtype=np.int32 )
    nnodes = np.array( [*elennodes[0] , *elennodes[1]], dtype=np.int32 )
    eles.cn_off = np.zeros((nnodes.shape[0]+1,), dtype=np.int32)
    eles.cn_off[1:] = np.cumsum(nnodes)
    eles.ptags = np.array( phystags[idxs[1]] , dtype=np.int32 )
    eles.ptagsN = np.array( phystags_n[idxs[1]] , dtype=np.int32 )

    return eles

def getFaces2D(eles, nodes):
    nc_con = []
    for i in range(nodes.shape[0]):
        nc_con.append([])
    for i in range(eles.N):
        nds = eles.cn_con[eles.cn_off[i]:eles.cn_off[i+1]]
        for nd in nds:
            nc_con[nd].append(i)
    
    fn_con = []
    fc_con = []
    
    for i in range(eles.NC):
        nds = eles.cn_con[ eles.cn_off[i]:eles.cn_off[i+1] ]
        nnds = nds.shape[0]
        for j in range(nnds-1):
            ndcon1 = nc_con[nds[j]]
            for k in range(j+1, nnds):
                ndcon2 = nc_con[nds[k]]
                comm = list(set(ndcon1) & set(ndcon2))
                if len(comm) == 2:
                    if min(comm) == i: continue
                    fn_con.append(nds[j])
                    fn_con.append(nds[k])
                    fc_con.append(comm[0])
                    fc_con.append(comm[1])

    for i in range(eles.NC, eles.N):
        nds = eles.cn_con[ eles.cn_off[i]:eles.cn_off[i+1] ]
        fn_con.append(nds[0])
        fn_con.append(nds[1])

        ndcon1 = nc_con[nds[0]]
        ndcon2 = nc_con[nds[1]]

        comm = list(set(ndcon1) & set(ndcon2))
        comm.sort()
        
        fc_con.append(comm[0])
        fc_con.append(comm[1])
    
    eles.Nf = int(len(fc_con)/2)
    eles.fn_con = np.array(fn_con, dtype=np.int32)
    eles.fc_con = np.array(fc_con, dtype=np.int32)

def gen_cf_con(eles):
    cf_con = []
    for i in range(eles.N):
        cf_con.append([])
    for i in range(eles.Nf):
        cf_con[eles.fc_con[2*i+0]].append(i)
        cf_con[eles.fc_con[2*i+1]].append(i)

    c_nf = np.array( [len(cf_con_i) for cf_con_i in cf_con], dtype=np.int32 )
    cf_off = np.zeros((eles.N+1,), dtype=np.int32)
    cf_off[1:] = np.cumsum(c_nf)

    cf_con = np.array([con for cell_con in cf_con for con in cell_con], dtype=np.int32)

    eles.cf_con = cf_con
    eles.cf_off = cf_off

def sortfcols(eles):
    NC  = eles.NC
    Nf  = eles.Nf
    NBC = eles.N-eles.NC

    fn_con = eles.fn_con[:2*(Nf-NBC)]
    fc_con = eles.fc_con[:2*(Nf-NBC)]

    cf_con = []
    for i in range(NC):
        cf_con.append([])
    for i in range(Nf-NBC):
        cf_con[fc_con[2*i+0]].append(i)
        cf_con[fc_con[2*i+1]].append(i)
    
    Ncols = 6
    colopts = np.array(range(Ncols), dtype=np.int32)
    colcount = np.zeros((Ncols,), dtype=np.int32)
    cols = -np.ones((Nf-NBC,), dtype=np.int32)
    for i in range(NC):
        faces = np.array(cf_con[i])
        fcols = list(cols[faces])
        cleanfaces = faces[fcols == np.int64(-1)]
        N_cleanfaces = len(cleanfaces) 
        freecols = [col for col in colopts if col not in fcols]
        cols[cleanfaces] = freecols[:N_cleanfaces]
        colcount[freecols[:N_cleanfaces]] += 1
    colcount = [col for col in colcount if col != 0]

    sortIDX = np.argsort(cols)
    fn_con = fn_con.reshape((-1,2))[sortIDX,:]
    fc_con = fc_con.reshape((-1,2))[sortIDX,:]

    colIDX = np.array([0,*np.cumsum(colcount)], dtype=np.int32)
    for i in range(len(colIDX)-1):
        subarr = np.sum(fc_con[colIDX[i]:colIDX[i+1],:], axis=1)
        subsortIDX = np.argsort(subarr)
        fn_con[colIDX[i]:colIDX[i+1],:] = fn_con[colIDX[i]+subsortIDX,:]
        fc_con[colIDX[i]:colIDX[i+1],:] = fc_con[colIDX[i]+subsortIDX,:]

    eles.fn_con[:2*(Nf-NBC)] = fn_con.ravel()
    eles.fc_con[:2*(Nf-NBC)] = fc_con.ravel()
    eles.colcount = colcount

def gen_regions(eles):
    Ncols = len(eles.colcount)
    tags = eles.ptags
    BCcount = eles.ptagsN
    NBCs = eles.ptags.shape[0]
    N_reg = Ncols + NBCs

    regns = np.zeros((N_reg,), dtype=np.int32)
    regntags = np.zeros((N_reg,), dtype=np.int32)
    regntags[Ncols:] = tags
  
    regns[:Ncols] = eles.colcount
    regns[Ncols:] = BCcount

    eles.regcount = regns
    eles.regtags = regntags 

def plotmesh(eles, nodes):
    centres = np.zeros((eles.N,2))
    for i in range(eles.N):
        centres[i,:] = np.mean( nodes[ eles.cn_con[ eles.cn_off[i]:eles.cn_off[i+1] ],:] , axis=0)
    
    # plt.plot(centres[:,0], centres[:,1], 'ko', ms=1)  

    fn_con = np.reshape(eles.fn_con, (-1,2))
    fc_con = np.reshape(eles.fc_con, (-1,2))
    x = 0.5*( nodes[fn_con[:,0],:] + nodes[fn_con[:,1],:] )
    fc = np.zeros((2*eles.Nf,2))
    fc[::2]  = nodes[fn_con[:,0],:]
    fc[1::2] = nodes[fn_con[:,1],:]
    rL = np.zeros((2*eles.Nf,2))
    rL[::2]  = centres[fc_con[:,0],:]
    rL[1::2] = x
    rR = np.zeros((2*eles.Nf,2))
    rR[::2]  = centres[fc_con[:,1],:]
    rR[1::2] = x
    
    for i in range(eles.Nf):
        I1 = 2*i+0; I2 = 2*i+1
        plt.plot([fc[I1,0], fc[I2,0]], [fc[I1,1], fc[I2,1]], 'k', linewidth=0.8)
        plt.plot([rL[I1,0], rL[I2,0]], [rL[I1,1], rL[I2,1]], 'r', linewidth=0.8)
        plt.plot([rR[I1,0], rR[I2,0]], [rR[I1,1], rR[I2,1]], 'g', linewidth=0.8)

    plt.savefig("temp.png", dpi=500)

def plotfmat(eles, nm):
    NBC = eles.N-eles.NC
    data = np.ones(eles.fc_con.shape, dtype=np.int32)
    fc_off = np.array([0, *np.cumsum( 2*np.ones((eles.Nf,)) )], dtype=np.int32)
    mat = sp.sparse.csr_matrix((data, eles.fc_con, fc_off), shape=(eles.Nf,eles.N)).tocoo()
    plt.plot(mat.col[:], mat.row[:], 's', color='red', ms=1)
    # plt.hlines(np.cumsum(eles.colcount), 0, eles.N)
    plt.hlines(eles.Nf-NBC, 0, eles.N, 'k')
    plt.vlines(eles.NC, 0, eles.Nf, 'k')
    plt.gca().invert_yaxis()
    plt.savefig(nm)

def writeHMSH(fname, physnames, nodes, eles):
    Nphysnames = physnames.shape[0]
    fmt = '<' + Nphysnames*'i'
    nameNs = np.array([len(physname) for physname in physnames], dtype=np.int32) + 1
    NphysnamesData = struct.pack(fmt, *nameNs)
    physnames_bytes = b''.join([s.encode('utf-8') + b'\0' for s in physnames])

    Ndata = struct.pack(
        '<iiiiiiiii', eles.dim, nodes.shape[0], 
        eles.N, eles.NC, eles.cn_con.shape[0], eles.Nf, 
        eles.cf_con.shape[0], Nphysnames, eles.regcount.shape[0]
    )

    with open(fname, "wb") as file:
        file.write(Ndata)
        file.write(NphysnamesData)
        file.write(physnames_bytes)
        file.write(eles.regcount.tobytes())
        file.write(eles.regtags.tobytes())
        file.write(nodes.tobytes()) 
        file.write(eles.cn_off.tobytes())
        file.write(eles.cn_con.tobytes())
        file.write(eles.fn_con.tobytes())
        file.write(eles.fc_con.tobytes())
        file.write(eles.cf_off.tobytes())
        file.write(eles.cf_con.tobytes())
        file.write(eles.eletypes.tobytes())

def process_import(args):
    fin = args.meshin

    secs = ["Entities", "PhysicalNames", "Nodes", "Elements"]
    stags = [f"${sec}" for sec in secs]
    etags = [f"$End{sec}" for sec in secs]

    secsraw = {}

    with open(fin, 'r') as rfile:
        file = rfile.read()
        
        for sec,stag,etag in zip(secs,stags,etags):
            secsraw[sec] = list( filter(None, file[file.index(stag)+len(stag) : file.index(etag)].split('\n') ) )
    
    physnames = parsePhysNames(secsraw["PhysicalNames"])
    curves, dim = parseEntities(secsraw["Entities"])
    nodes = parseNodes(secsraw["Nodes"])
    if dim == 2: 
        eles = parseElements2D(secsraw["Elements"], curves)
        getFaces2D(eles, nodes)
        sortfcols(eles)
        gen_cf_con(eles)
        gen_regions(eles)
        writeHMSH(args.meshout, physnames, nodes, eles)
    else: 
        raise Exception("3D not yet supported")