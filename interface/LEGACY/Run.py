import numpy as np
import configparser as cfgp
from ctypes import byref, POINTER, c_void_p, c_double, c_int
from interface.Global import Struct
from interface.PST import load_mesh

def get_ICs(cfg, nodes, eles):
    gma = float(cfg['consts']['gamma'])
    gmo = gma-1

    r = lambda x,y: eval(cfg['ICs']['r'])
    u = lambda x,y: eval(cfg['ICs']['u'])
    v = lambda x,y: eval(cfg['ICs']['v'])
    p = lambda x,y: eval(cfg['ICs']['p'])

    centres = np.zeros((eles.NC,2))
    for i in range(eles.NC):
        centres[i,:] = np.mean( nodes[ eles.cn_con[ eles.cn_off[i]:eles.cn_off[i+1] ],:] , axis=0)
    
    P = np.zeros((eles.NC, 4))
    P[:,0] = r(centres[:,0], centres[:,1])
    P[:,1] = u(centres[:,0], centres[:,1])
    P[:,2] = v(centres[:,0], centres[:,1])
    P[:,3] = p(centres[:,0], centres[:,1])

    U = np.zeros((eles.NC, 4))
    U[:,0] = P[:,0]
    U[:,1] = P[:,0]*P[:,1]
    U[:,2] = P[:,0]*P[:,2]
    U[:,3] = P[:,3]/gmo + 0.5*P[:,0]*(P[:,1]*P[:,1] + P[:,2]*P[:,2])
    U = U.ravel()

    U = U.ctypes.data_as(POINTER(c_double))
    return U

def get_consts(cfg, eles):
    consts = Struct()
    consts.nvars = eles.dim+2
    consts.ndims = eles.dim
    consts.gma = float(cfg["consts"]["gamma"])
    if cfg["solver"]["system"].lower() == "euler":
        consts.mu = 0
        consts.Pr = 0
    else:
        consts.mu = float(cfg["consts"]["mu"])
        consts.Pr = float(cfg["consts"]["Pr"])
    return consts

def get_BCs(cfg, physnames, eles):
    sections = cfg.sections()
    BCnames = [sec.replace('BC-', '') for sec in sections if 'BC-' in sec]

    Nf  = int(eles.Nf)
    NBC = int(eles.N-eles.NC)

    physIDXs = {}
    tag = int(eles.ptags[0])
    physIDXs[tag] = [(Nf-NBC), None]

    for f in range(NBC):
        if eles.ptags[f] != tag:
            physIDXs[tag][1] = f + (Nf-NBC)
            tag = int(eles.ptags[f])
            physIDXs[tag] = [f + (Nf-NBC), None]
    physIDXs[tag][1] = Nf

    BCs = []
    for BCname in BCnames:
        phystag = physnames.index(BCname)+1
        F = physIDXs[phystag]

    return BCs 

def prepareMesh(eles):
    eles.fc_con = eles.fc_con.ctypes.data_as(POINTER(c_int))
    eles.fn_con = eles.fn_con.ctypes.data_as(POINTER(c_int))

    eles.cf_con = eles.cf_con.ctypes.data_as(POINTER(c_int))
    eles.cf_off = eles.cf_off.ctypes.data_as(POINTER(c_int))

    eles.cn_con = eles.cn_con.ctypes.data_as(POINTER(c_int))
    eles.cn_off = eles.cn_off.ctypes.data_as(POINTER(c_int))

def init_regs(kernels, nodes, eles, ICs):
    regs = {
        "faces" : c_void_p(),
        "U"     : c_void_p(),
        "F"     : c_void_p(),
        "G"     : c_void_p(),
        "phi"   : c_void_p(),
        "ivols" : c_void_p(),
        "UfL"   : c_void_p(),
        "UfR"   : c_void_p()
    }

    nodes_in = nodes.ravel().ctypes.data_as(POINTER(c_double))

    kernels.Initialise(
        eles.N, eles.NC, eles.Nf, 4,
        ICs, nodes_in, 
        eles.fc_con, eles.fn_con,
        eles.cn_con, eles.cn_off,
        eles.cf_con, eles.cf_off,
        byref(regs["faces"]), byref(regs["U"   ]), 
        byref(regs["G"    ]), byref(regs["phi" ]), 
        byref(regs["ivols"])
    )

    kernels.CreateReg(eles.Nf, eles.dim+2, eles.dim, byref(regs["UfL"]))
    kernels.CreateReg(eles.Nf, eles.dim+2, eles.dim, byref(regs["UfR"]))
    kernels.CreateReg( eles.N, eles.dim+2, eles.dim, byref(regs["F"])  )

    return regs

def getTimeIntegrator(cfg):
    scheme = cfg["time-integrator"]["scheme"].lower()
    T = float( cfg["time-integrator"]["T"] )
    dt = float( cfg["time-integrator"]["dt"] )

    return None, None, None

def finaliseTEMP(regs, kernels):
    for reg in regs.values():
        kernels.Free(reg)

def process_run(args):
    cfg = cfgp.ConfigParser()
    cfg.read(args.config)

    physnames, nodes, eles = load_mesh(args.mesh)
    plugins, kernels = None, None
    ICs = get_ICs(cfg, nodes, eles)
    consts = get_consts(cfg, eles)
    BCs = get_BCs(cfg, physnames, eles)

    prepareMesh(eles)
    Intg, T, dt = getTimeIntegrator(cfg)
    regs = init_regs(kernels, nodes, eles, ICs)

    intg = Intg(consts, eles, kernels, plugins, BCs, regs)
    intg.run(dt, T)
