import numpy as np
from scipy import sparse
import struct
from dataclasses import dataclass
import matplotlib.pyplot as plt

gmsh_eltype_nnds = np.array([0, 2, 3, 4, 4, 8, 6, 5], dtype=np.int32)
gmsh_eltype_nsds = np.array([0, 1, 3, 4, 4, 6, 5, 5], dtype=np.int32)

@dataclass
class CellBlock:
    ptag: int
    celltype: int
    ncells: int
    cn_con: np.ndarray

    def __repr__(self):
        return f"ptag: {self.ptag}, celltype: {self.celltype}, ncells: {self.ncells}, cn_con shape: {self.cn_con.shape}\n"


class mshParserASCII:    
    def load_mesh(self, fname):
        self._readfile(fname)

        self._parsePhysNames()
        self._parseEntities()
        self._parseNodes()
        self._parseElementBlocks()
        
        return self.dim, self.n_nodes, self.n_blocks, self.n_cells, self.nodes, self.blocks, self.physnames
    
    def _readfile(self, fname):
        secs = np.array(["Entities", "PhysicalNames", "Nodes", "Elements"])
        stags = np.array([f"${sec}" for sec in secs])
        etags = np.array([f"$End{sec}" for sec in secs])

        self.secsraw = {}

        with open(fname, 'r') as file:
            file = file.read()
            
            for sec,stag,etag in zip(secs,stags,etags):
                I1 = np.char.find(file, stag) + len(stag) + 1
                I2 = np.char.find(file, etag) - 1
                self.secsraw[sec] = file[I1:I2]

    def _parsePhysNames(self):
        sec = self.secsraw["PhysicalNames"].split('\n')
        n_physnames = int(sec[0])

        physnames = []
        for i in range(1,n_physnames):
            ln = sec[i].split(' ')
            name = ln[2][1:-1]
            physnames.append( name ) 
        self.physnames = physnames

    def _parseEntities(self):
        sec = self.secsraw["Entities"].split('\n')
        npts, ncurves, nsurfs, nvols = [int(entity) for entity in sec[0].split(' ')]

        if nvols == 0: self.dim = 2
        else: self.dim = 3
        if nsurfs == 0: raise Exception("No surface. ")

        self.curves = [None]
        for i in range(npts+1,npts+1+ncurves):
            line = sec[i].split(' ')
            if int(line[7]) > 1: raise Exception("Invalid number of physical tags on curve. ")
            if int(line[7]) == 0: self.curves.append(None)
            else: self.curves.append( int(line[8]) )

    def _parseNodes(self):
        sec = self.secsraw["Nodes"]
        sec = np.char.replace(sec, '\n', ' ')
        sec = np.char.split(sec, ' ').item()
        sec = list(filter(None, sec))

        n_nodes = int(sec[1])

        nodes = np.zeros((n_nodes,3))
        ndhead = 0
        head = 4
        while head < len(sec):
            nblock = sec[head:]
            para = int(nblock[2])
            if para == 1: raise Exception("Parametric nodes not supported. ")
            n_nodesinblock = int(nblock[3])

            nds = np.array(nblock[4+n_nodesinblock:4+4*n_nodesinblock], dtype=np.float64)
            nds = np.reshape(nds, (-1,3))

            nodes[ndhead:ndhead+n_nodesinblock,:] = nds[:,:]
            
            ndhead += n_nodesinblock
            head = head + 4+4*n_nodesinblock

        self.n_nodes = n_nodes
        self.nodes = nodes[:,:self.dim]

    def _parseElementBlocks(self):
        sec = self.secsraw["Elements"]
        sec = np.char.replace(sec, '\n', ' ')
        sec = np.char.split(sec, ' ').item()
        sec = list(filter(None, sec))
        sec = np.array(sec, dtype=np.int32)
        
        self.n_blocks = sec[0]
        self.n_cells = sec[1]

        self.blocks = []
        head = 4
        block = 0
        while head < sec.shape[0]:
            eblock = sec[head:]
            endim,entag,eltype,n_cells = eblock[:4]
            nnds = gmsh_eltype_nnds[eltype]

            nds = eblock[4:4+n_cells*(1+nnds)].reshape((-1,1+nnds))
            nds = nds[:,1:].ravel()

            if endim == 1: ptag = self.curves[entag]
            else: ptag = 0

            self.blocks.append(CellBlock(ptag, eltype, n_cells, nds-1))
            head = head+4+n_cells*(1+nnds)
            block += 1

"""
class mshParserBIN:    
    def load_mesh(self, fname):
        self._readfile(fname)

        print(self.secsraw["Entities"])

        # self._parsePhysNames()
        # self._parseEntities()
        # self._parseNodes()
        # self._parseElementBlocks()

        # return self.dim, self.n_nodes, self.n_blocks, self.n_cells, self.nodes, self.blocks, self.physnames
    
    def _readfile(self, fname):
        secs = np.array(["Entities", "PhysicalNames", "Nodes", "Elements"])
        stags = np.array([f"${sec}".encode() for sec in secs])
        etags = np.array([f"$End{sec}".encode() for sec in secs])

        self.secsraw = {}

        with open(fname, 'rb') as file:
            file = file.read()
            
            for sec,stag,etag in zip(secs,stags,etags):
                I1 = np.char.find(file, stag) + len(stag) + 1
                I2 = np.char.find(file, etag) - 1
                self.secsraw[sec] = file[I1:I2]
"""

class Mesh:
    def generate_from_msh(self, fname):
        parser = mshParserASCII()
        self.dim, self.n_nodes, self.n_blocks, self.n_cells, self.nodes, self.blocks, self.physnames = parser.load_mesh(fname)

        self._sortBlocks()
        self._combineBlocks()
        self._expandBlocks()

        self._gen_cc_con()
        self._rcm_sort_mesh()
        self._gen_faces()

    def load_hmsh(self, fname):
        sz_chr = 1
        sz_int = 4
        sz_dbl = 8
        
        with open(fname, "rb") as file:
            bin = file.read()
            head = bin
        
        Ndat = np.frombuffer(head[:sz_int*8], dtype=np.int32)
        head = head[sz_int*8:]
        dim         = Ndat[0]
        n_nodes     = Ndat[1]
        n_blocks    = Ndat[2]
        n_cells     = Ndat[4]
        n_faces     = Ndat[5]
        n_chars     = Ndat[7]

        physnames_b   = head[:n_chars*sz_chr ]                                     ;  head = head[n_chars*sz_chr     :]
        self.fcoltag  = np.frombuffer(head[:n_blocks*sz_int    ], dtype=np.int32  );  head = head[n_blocks*sz_int    :]
        self.fcolcnt  = np.frombuffer(head[:(n_blocks+1)*sz_int], dtype=np.int32  );  head = head[(n_blocks+1)*sz_int:]
        self.nodes    = np.frombuffer(head[:n_nodes*dim*sz_dbl ], dtype=np.float64);  head = head[n_nodes*dim*sz_dbl :]
        self.fn_con   = np.frombuffer(head[:2*n_faces*sz_int   ], dtype=np.int32  );  head = head[2*n_faces*sz_int   :]
        self.fc_con   = np.frombuffer(head[:2*n_faces*sz_int   ], dtype=np.int32  );  head = head[2*n_faces*sz_int   :]
        self.cn_off   = np.frombuffer(head[:(n_cells+1)*sz_int ], dtype=np.int32  );  head = head[(n_cells+1)*sz_int :]
        n_cn_con = self.cn_off[-1]
        self.cn_con   = np.frombuffer(head[:n_cn_con*sz_int    ], dtype=np.int32  );  head = head[n_cn_con*sz_int    :]
        self.cc_off   = np.frombuffer(head[:(n_cells+1)*sz_int ], dtype=np.int32  );  head = head[(n_cells+1)*sz_int :]
        n_cc_con = self.cc_off[-1]
        self.cc_con   = np.frombuffer(head[:n_cc_con*sz_int    ], dtype=np.int32  );  head = head[n_cc_con*sz_int    :]
        self.c_ctypes = np.frombuffer(head[:n_cells*sz_int     ], dtype=np.int32  );  head = head[n_cells*sz_int     :]

        self.nodes      = self.nodes.reshape((-1,dim))
        self.physnames  = np.array([s for s in physnames_b.decode().split('\x00') if s])
        self.dim        = dim
        self.n_nodes    = n_nodes
        self.n_blocks   = n_blocks
        self.n_cells    = n_cells
        self.n_faces    = n_faces
    

    def _sortBlocks(self):
        ptags = []
        for i in range(self.n_blocks):
            ptags.append(self.blocks[i].ptag)

        sortIDX = np.argsort(ptags)
        blocks_sorted = []
        for i in range(self.n_blocks):
            blocks_sorted.append(self.blocks[sortIDX[i]])

        self.blocks = blocks_sorted

    def _combineBlocks(self): # combine like blocks
        ptags = []
        for i in range(self.n_blocks):
            ptags.append(self.blocks[i].ptag)

        blocks_c = [self.blocks[0]]
        acc = self.n_blocks*[0]
        acc[0] = 1
        block = 0
        while True:
            for i in range(1,self.n_blocks):
                curr_block = self.blocks[i]
                same_ptag = (curr_block.ptag == blocks_c[-1].ptag)
                same_celltype = (curr_block.celltype == blocks_c[-1].celltype)
                not_same_block = (i != block)

                if same_ptag and same_celltype and not_same_block:
                    blocks_c[-1].ncells += curr_block.ncells
                    blocks_c[-1].cn_con = np.concat((blocks_c[-1].cn_con, curr_block.cn_con))
                    acc[i] = 1

            if 0 not in acc: break

            block = acc.index(0)
            blocks_c.append(self.blocks[block])
            acc[block] = 1

        self.n_blocks  = len(blocks_c)
        self.n_bc_fcol = self.n_blocks-1
        self.blocks    = blocks_c

    def _expandBlocks(self): # blocks to arrays
        self.c_ptags  = np.zeros((self.n_cells,), dtype=np.int32)
        self.c_ctypes = np.zeros((self.n_cells,), dtype=np.int32)
        c_nnds   = np.zeros((self.n_cells,), dtype=np.int32)

        head = 0
        len_cn_con = 0
        for block in self.blocks:
            len_cn_con += block.ncells * gmsh_eltype_nnds[block.celltype]

            tail = head + block.ncells
            self.c_ptags[head:tail] = block.ptag
            self.c_ctypes[head:tail] = block.celltype
            c_nnds[head:tail] = gmsh_eltype_nnds[block.celltype]
            head = tail

        self.cn_off = np.zeros((self.n_cells+1,), dtype=np.int32)
        self.cn_off[1:] = np.cumsum(c_nnds)
        self.cn_con = np.zeros((len_cn_con,), dtype=np.int32)

        head = 0
        for block in self.blocks:
            tail = head + block.ncells * gmsh_eltype_nnds[block.celltype]
            self.cn_con[head:tail] = block.cn_con
            head = tail


    def _gen_cc_con(self):
        data_cn = np.ones_like(self.cn_con)
        cn = sparse.csr_matrix((data_cn, self.cn_con, self.cn_off))
        cc = cn @ cn.T
        cc.setdiag(0)

        mask = cc.data >= self.dim
        cc.data = mask
        cc.eliminate_zeros()

        self.cc_con = cc.indices
        self.cc_off = cc.indptr
    
    def _rcm_sort_mesh(self):
        data_cc = np.ones_like(self.cc_con)
        cc = sparse.csr_matrix((data_cc, self.cc_con, self.cc_off))
        c_perm = sparse.csgraph.reverse_cuthill_mckee(cc)
        cc = cc[c_perm,:][:,c_perm]

        data_cn = np.ones_like(self.cn_con)
        cn = sparse.csr_matrix((data_cn, self.cn_con, self.cn_off))
        cn = cn[c_perm,:]

        nn = cn.T @ cn
        nn.setdiag(0)
        n_perm = sparse.csgraph.reverse_cuthill_mckee(nn)[::-1]
        nn = nn[n_perm,:][:,n_perm]
        cn = cn[:,n_perm]
        self.nodes = self.nodes[n_perm,:]

        self.cc_con   = cc.indices
        self.cc_off   = cc.indptr
        self.c_ptags  = self.c_ptags[c_perm]
        self.c_ctypes = self.c_ctypes[c_perm]

        self.cn_con = cn.indices
        self.cn_off = cn.indptr

    def _gen_faces(self):
        # find fc_con
        self.n_faces = int(self.cc_con.shape[0]/2)
        data_cc = np.ones_like(self.cc_con)
        cc = sparse.csr_matrix((data_cc, self.cc_con, self.cc_off)).tocoo()
        mask = cc.col > cc.row
        cl = cc.row[mask]
        cr = cc.col[mask]
        fc_con = np.array([cl,cr], dtype=np.int32).transpose()

        # find ptags
        fl_ptags = self.c_ptags[cl]
        fr_ptags = self.c_ptags[cr]
        f_ptags = np.max([fl_ptags, fr_ptags], axis=0)

        # swap cells so right cell is always BC cell
        idx = np.where(fl_ptags > fr_ptags)[0]
        fc_con[idx] = fc_con[idx][:, [1, 0]]
        cl = fc_con[:,0]
        cr = fc_con[:,1]

        # find fn_con
        data = np.ones_like(self.cn_con, dtype=np.int32)
        cn = sparse.csr_matrix((data, self.cn_con, self.cn_off)).tocsc() # , (self.n_cells, self.n_nodes)
        off = np.arange(self.n_faces+1, dtype=np.int32)
        data = np.ones_like(cl, dtype=np.int32)
        fcl = sparse.csr_matrix((data, cl, off), (self.n_faces, self.n_cells))
        fcr = sparse.csr_matrix((data, cr, off), (self.n_faces, self.n_cells))
        fnl = fcl @ cn
        fnr = fcr @ cn
        fn = fnl.multiply(fnr)
        fn.eliminate_zeros()
        fn_con = fn.indices.reshape((-1,2))

        # sort by ptag
        self.fc_con = np.zeros_like(fc_con, dtype=np.int32)
        self.fn_con = np.zeros_like(fn_con, dtype=np.int32)
        self.fcoltag = np.zeros((self.n_blocks,), dtype=np.int32)
        self.fcolcnt = np.zeros((self.n_blocks+1,), dtype=np.int32)

        head = 0
        for block in range(self.n_blocks):
            idx = np.where(f_ptags == block)[0]
            tail = head+idx.shape[0]

            self.fc_con[head:tail,:] = fc_con[idx,:]
            self.fn_con[head:tail,:] = fn_con[idx,:]
            self.fcoltag[block]   = block
            self.fcolcnt[block+1] = tail
            head = tail

        # print(self.fcolcnt[1])
        # plt.plot(self.fc_con[:self.fcolcnt[1],0], 'k')
        # plt.plot(self.fc_con[:self.fcolcnt[1],1], 'r')
        # plt.savefig("t.png")

        self.fc_con = self.fc_con.ravel()
        self.fn_con = self.fn_con.ravel()        

    def plot_mesh(self, fname):
        plt.clf()
        nodes = self.nodes

        ndcount = self.cn_off[1:] - self.cn_off[:-1]
        mat = sparse.csr_matrix((np.ones_like(self.cn_con), self.cn_con, self.cn_off))
        xc = mat @ self.nodes
        xc[:,0] = xc[:,0]/ndcount
        xc[:,1] = xc[:,1]/ndcount
    
        cl = self.fc_con[0::2]
        cr = self.fc_con[1::2]

        nl = nodes[self.fn_con[0::2],:]
        nr = nodes[self.fn_con[1::2],:]
        xf = 0.5*(nl+nr)

        rl = xf[:,:] - xc[cl,:]
        xc[cr[2550:],:] = xf[2550:,:] + rl[2550:,:]
        rr = xf[:,:] - xc[cr,:]

        Af = np.zeros((self.n_faces,2))
        Af[:,0] =  (nr[:,1]-nl[:,1])
        Af[:,1] = -(nr[:,0]-nl[:,0])

        Af_dot_rl = Af[:,0]*rl[:,0] + Af[:,1]*rl[:,1]
        cond = Af_dot_rl < 0
        Af[cond,:] = - Af[cond,:] 

        plt.plot([nl[ :,0], nr[:,0]], [nl[ :,1], nr[:,1]], 'k', linewidth=0.5)
        plt.plot([xc[cl,0], xf[:,0]], [xc[cl,1], xf[:,1]], 'r', linewidth=0.3)
        plt.plot([xc[cr,0], xf[:,0]], [xc[cr,1], xf[:,1]], 'b', linewidth=0.3)
        # plt.plot([xf[:,0], xf[:,0]+Af[:,0]], [xf[:,1], xf[:,1]+Af[:,1]], 'g', linewidth=0.3)
        plt.plot(nodes[:,0], nodes[:,1], 'kx', ms=0.8)
        plt.plot(xc[:,0], xc[:,1], 'r+', ms=0.8)
        # plt.gca().set_aspect('equal')
        plt.savefig(fname, dpi=1000)

    def plot_BC_faces(self, fname):
        plt.clf()
        nodes = self.nodes
        int_fcs = np.where(self.fcoltag == 0)[0][-1]

        F1 = self.fcolcnt[int_fcs+1]
        F2 = self.fcolcnt[-1]

        fn_con = self.fn_con[2*F1:2*F2]
        nl = nodes[fn_con[0::2],:]
        nr = nodes[fn_con[1::2],:]

        plt.plot([nl[ :,0], nr[:,0]], [nl[ :,1], nr[:,1]], 'kx-', linewidth=1)
        plt.gca().set_aspect('equal')
        plt.savefig(fname, dpi=1000)

    def plot_mat(self, con, off, fname):
        plt.clf()
        data = np.ones_like(con)
        mat = sparse.csr_matrix((data, con, off))
        plt.spy(mat, precision=0, marker='s', color='red', markersize=1)
        plt.savefig(fname, dpi=800)

    def write_hmsh(self, fname):
        n_physnames = len(self.physnames)
        physnames_b = b''.join([s.encode('utf-8') + b'\0' for s in self.physnames])
        n_chars = len(physnames_b)
        
        Ndata = struct.pack(
            '<iiiiiiii', self.dim, self.n_nodes, self.n_blocks, 
            self.n_bc_fcol, self.n_cells, self.n_faces, n_physnames, n_chars
        )

        with open(fname, "wb") as file:
            file.write(Ndata)
            file.write(physnames_b)
            file.write(self.fcoltag.tobytes())
            file.write(self.fcolcnt.tobytes())
            file.write(self.nodes.tobytes()) 
            file.write(self.fn_con.tobytes())
            file.write(self.fc_con.tobytes())
            file.write(self.cn_off.tobytes())
            file.write(self.cn_con.tobytes())
            file.write(self.cc_off.tobytes())
            file.write(self.cc_con.tobytes())
            file.write(self.c_ctypes.tobytes())


def process_import(args):
    fin  = args.meshin
    fout = args.meshout

    mesh = Mesh()
    mesh.generate_from_msh(fin)
    mesh.write_hmsh(fout)

def process_checkcell(args):
    fin  = args.meshin
    cell = int(args.cell)

    mesh = Mesh()
    mesh.load_hmsh(fin)

    ctype = mesh.c_ctypes[cell]
    nodes = mesh.cn_con[mesh.cn_off[cell]:mesh.cn_off[cell+1]]
    cells = mesh.cc_con[mesh.cc_off[cell]:mesh.cc_off[cell+1]]
    xc    = np.mean(mesh.nodes[nodes,:], axis=0)

    print(f"Cell: {cell}:")
    print(f"\ttype:       {ctype}")
    print(f"\tnodes:      {nodes}")
    print(f"\tneighbours: {cells}")
    print(f"\tcentre:     {xc}")

def process_plot_BC_faces(args):
    fin  = args.meshin
    fout = args.imgout

    mesh = Mesh()
    mesh.load_hmsh(fin)
    mesh.plot_BC_faces(fout)
