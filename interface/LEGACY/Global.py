import numpy as np
import struct

ctype_names = [None, 'line', 'tri', 'quad', 'tet', 'hex', 'pri', 'pyr']
ctype_nnodes = {'line':2, 'tri':3, 'quad':4, 'tet':4, 'hex':8, 'pri':6, 'pyr':5}
ctype_nsides = {'line':1, 'tri':3, 'quad':4, 'tet':4, 'hex':6, 'pri':5, 'pyr':5}

class Struct:
    pass

class Elements:
    def __init__(self):
        self.dim: int
        
        self.N: int
        self.NC: int
        self.Nf: int

        self.eletypes: list[str]

        self.cn_con: np.ndarray
        self.cn_off: np.ndarray

        self.fn_con: np.ndarray
        self.fn_off: np.ndarray
        self.fc_con: np.ndarray

        self.ptags: np.ndarray
        self.colcount: np.ndarray

def read_soln(fname):
    with open(fname, "rb") as file:
        bin = file.read()
    N = struct.unpack_from("<i", bin)[0]
    bindata = np.frombuffer(bin[4:], dtype=np.float64)
    bindata = np.reshape(bindata, (4,N))
    bindata = np.transpose(bindata)

    fdata = np.zeros((N,4), dtype=np.float64)
    fdata[:,0] = bindata[:,0]
    fdata[:,1] = bindata[:,1]/fdata[:,0]
    fdata[:,2] = bindata[:,2]/fdata[:,0]
    fdata[:,3] = 0.4 * (bindata[:,3] - 0.5*fdata[:,0]*(fdata[:,1]*fdata[:,1] + fdata[:,2]*fdata[:,2]))
    return fdata

def load_mesh(fname):
    sz_int = 4
    sz_dbl = 8
    
    with open(fname, "rb") as file:
        bin = file.read()
        head = bin

    Nr = np.frombuffer(head[:sz_int*9], dtype=np.int32); head = head[sz_int*9:]
    dim   = Nr[0]; Nnds  = Nr[1]
    N     = Nr[2]; NC    = Nr[3]
    Ncnc  = Nr[4]; Nf    = Nr[5]
    Ncfc  = Nr[6]; Npnms = Nr[7]
    Nregs = Nr[8]

    Npname = np.frombuffer(head[:sz_int*Npnms], dtype=np.int32); head = head[sz_int*Npnms:]
    Nchars = np.sum(Npname)

    eles = Elements()
    eles.N = N; eles.NC = NC; eles.Nf = Nf
    physnamesb = head[:Nchars]; head = head[Nchars:]
    eles.regcount = np.frombuffer(head[:Nregs*sz_int]   , dtype=np.int32  ); head = head[Nregs*sz_int:]
    eles.regtags  = np.frombuffer(head[:Nregs*sz_int]   , dtype=np.int32  ); head = head[Nregs*sz_int:]
    nodes         = np.frombuffer(head[:dim*Nnds*sz_dbl], dtype=np.float64); head = head[dim*Nnds*sz_dbl:]
    eles.cn_off   = np.frombuffer(head[:(N+1)*sz_int]   , dtype=np.int32  ); head = head[(N+1)*sz_int:]
    eles.cn_con   = np.frombuffer(head[:Ncnc*sz_int]    , dtype=np.int32  ); head = head[Ncnc*sz_int:]
    eles.fn_con   = np.frombuffer(head[:2*Nf*sz_int]    , dtype=np.int32  ); head = head[2*Nf*sz_int:]
    eles.fc_con   = np.frombuffer(head[:2*Nf*sz_int]    , dtype=np.int32  ); head = head[2*Nf*sz_int:]
    eles.cf_off   = np.frombuffer(head[:(N+1)*sz_int]   , dtype=np.int32  ); head = head[(N+1)*sz_int:]
    eles.cf_con   = np.frombuffer(head[:Ncfc*sz_int]    , dtype=np.int32  ); head = head[Ncfc*sz_int:]
    eles.eletypes = np.frombuffer(head[:N*sz_int]       , dtype=np.int32  ); head = head[N*sz_int:]

    physnames = np.array([s for s in physnamesb.decode().split('\x00') if s])
    nodes = np.reshape(nodes, (-1,2))
    return physnames, nodes, eles

