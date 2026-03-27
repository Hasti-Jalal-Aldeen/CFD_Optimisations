import glob
import re
import numpy as np
import vtk
from vtk.util import numpy_support
import struct

from interface.Mesh import Mesh

ctype_names = [None, 'line', 'tri', 'quad', 'tet', 'hex', 'pri', 'pyr']
ctype_nnodes = {'line':2, 'tri':3, 'quad':4, 'tet':4, 'hex':8, 'pri':6, 'pyr':5}
ctype_nsides = {'line':1, 'tri':3, 'quad':4, 'tet':4, 'hex':6, 'pri':5, 'pyr':5}

def read_soln(fname):
    sz_chr = 1
    sz_int = 4
    sz_dbl = 8

    with open(fname, "rb") as file:
        bin = file.read()
    N = struct.unpack_from("<iii", bin)
    Nc       = N[0]
    nvars    = N[1]
    nchars   = N[2]
    head     = bin[3*sz_int:]
    varnames = head[:nchars*sz_chr]
    head     = head[nchars*sz_chr:]
    data     = np.frombuffer(head[:nvars*Nc*sz_dbl], dtype=np.float64)
    data     = np.reshape(data, (nvars,Nc))

    varnames = varnames.split(b'\x00')[:-1]
    return varnames, data

def write_vtk(varnames: list[str], filename: str, mesh: Mesh, cell_data: np.ndarray):
    nodes3 = np.hstack([mesh.nodes, np.zeros((mesh.nodes.shape[0], 1))])
    vtk_points = vtk.vtkPoints()
    vtk_points.SetData(numpy_support.numpy_to_vtk(nodes3, deep=True))

    con = mesh.cn_con
    off = mesh.cn_off

    vtk_cell_type_map = np.array([0, vtk.VTK_LINE, vtk.VTK_TRIANGLE, vtk.VTK_QUAD], dtype=np.uint8)

    cell_types = vtk_cell_type_map[mesh.c_ctypes]
    vtk_types = numpy_support.numpy_to_vtk(cell_types, deep=True)

    vtk_cell_array = vtk.vtkCellArray()
    vtk_cell_array.SetData(
        numpy_support.numpy_to_vtk(off, deep=True),
        numpy_support.numpy_to_vtk(con, deep=True)
    )
    
    grid = vtk.vtkUnstructuredGrid()
    grid.SetPoints(vtk_points)
    grid.SetCells(vtk_types, vtk_cell_array)

    vel_mask = np.zeros((len(varnames),), dtype=bool)
    for var,varname in enumerate(varnames):
        if b"velocity" in varname.lower():
            vel_mask[var] = True
            continue
        vtk_fld = numpy_support.numpy_to_vtk(cell_data[var,:], deep=True)
        vtk_fld.SetName(varname)
        grid.GetCellData().AddArray(vtk_fld)

    if np.count_nonzero(vel_mask) > 0:
        if mesh.dim == 2:
            vel2d = cell_data[vel_mask,:]
            vel= np.vstack([vel2d, np.zeros((1,mesh.n_cells))])
        elif mesh.dim == 3:
            vel = cell_data[vel_mask,:]
        vel = vel.transpose()
        vtk_velocity = numpy_support.numpy_to_vtk(vel, deep=True)
        vtk_velocity.SetName("Velocity")
        grid.GetCellData().SetVectors(vtk_velocity)

    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetDataModeToBinary()
    writer.SetFileName(filename)
    writer.SetInputData(grid)
    writer.Write()

def process_export(args):
    meshnm = args.mesh
    solnin = args.solnin
    vtkout = args.solnout

    mesh = Mesh()
    mesh.load_hmsh(meshnm)
    varnames, fdata = read_soln(solnin)
    write_vtk(varnames, vtkout, mesh, fdata)

def process_batchexport(args):
    meshnm = args.mesh
    solnin_pattern = args.solnin_pattern
    vtkout_pattern = args.solnout_pattern
    Wnone = args.Wnone

    glob_pattern = solnin_pattern.replace('%d', '*')
    solnsin = glob.glob(glob_pattern)

    solnin_repattern_str = solnin_pattern.replace('%d', r'(\d+)')
    solnin_repattern = re.compile(solnin_repattern_str)

    mesh = Mesh()
    mesh.load_hmsh(meshnm)

    for solnin in solnsin:
        solnMatch = solnin_repattern.search(solnin)

        if (not solnMatch) and (not Wnone):
            print(f"Warning: Skipping '{solnin}' as it doesn't match the expected pattern.")
            continue
        
        varnames, fdata = read_soln(solnin)

        solnnum = solnMatch.group(1)
        vtkout = re.sub("%d", solnnum, vtkout_pattern)

        write_vtk(varnames, vtkout, mesh, fdata)

