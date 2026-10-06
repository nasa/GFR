from   typing import List
import math
import re

def read_obj_file(filename_obj: str):
    nodes = []
    faces = []
    try:
        file = open(filename_obj, 'r', encoding="ascii")
    except:
        print(f"Error: Unable to read file {filename_obj}.")
    with file:
        print(f"Reading file {filename_obj}.")
        line_num = 0
        for line in file:
            line_num += 1
            line = line.rstrip()  # strip trailing whitespace, eol
            if line.startswith('v '):
                # Process 'v x_flt y_flt z_flt'
                # split() returns List[str]. Convert to List[float].
                xyz_lst_str = line[2:].split()
                xyz_lst_flt = [float(x) for x in xyz_lst_str]
                nodes.append(xyz_lst_flt)
            elif line.startswith('f '):
                # Process 'f n1/i1/j1 n2/i2/j2 n3/i3/j3'
                tmp_lst_str = re.split(r'[ /]', line[2:])
                tmp         = [int(x) for x in tmp_lst_str]
                faces.append([tmp[0],tmp[3],tmp[6]])
    return [nodes,faces]

def write_tec_file(filename_tec: str,
                   nodes: List[List[float]],
                   faces: List[List[int]]) -> None :
    try:
        file = open(filename_tec, 'w', encoding="ascii")
    except:
        print(f"Error: Unable to write file {filename_tec}.")
    with file:
        print(f"Writing file {filename_tec}.")
        nnodes = len(nodes)
        nfaces = len(faces)
        file.write(f"Variables=\"x\",\"y\",\"z\"\n")
        file.write(f"Zone T=\"Surface\""
                   f", DATAPACKING=POINT"
                   f", NODES={nnodes}"
                   f", ELEMENTS={nfaces}"
                   f", ZONETYPE=FETRIANGLE\n")
        file.write(f"# Node data\n")
        for node in (nodes):
            x,y,z = node
            file.write(f"  {x:12.4f}{y:12.4f}{z:12.4f}\n")
        file.write(f"# Elem data\n")
        for face in (faces):
            n1,n2,n3 = face
            file.write(f"  {n1} {n2} {n3}\n")

def center_obj(nodes: List[List[float]]) -> List[List[float]] :
    # nodes enter as (x,y,z) = (span,up,streamwise)
    # nodes leave as (x,y,z) = (streamwise,up,span)
    #
    # get object bounds
    xmin = +1.0E+20
    xmax = -1.0E+20
    ymin = +1.0E+20
    ymax = -1.0E+20
    zmin = +1.0E+20
    zmax = -1.0E+20
    for n in range(len(nodes)):
        node = nodes[n]
        x,y,z = node
        xmin = min(xmin,x)
        xmax = max(xmax,x)
        ymin = min(ymin,y)
        ymax = max(ymax,y)
        zmin = min(zmin,z)
        zmax = max(zmax,z)
    obj_length = zmax-zmin
    obj_width  = xmax-xmin
    obj_height = ymax-ymin
    #
    # print object dimensions
    if True:
        print(f"length = {obj_length}")
        print(f"width  = {obj_width }")
        print(f"height = {obj_height}")
    #
    # set scaling factor
    fac = 1. / obj_length
    #
    # center the object in x=span
    x0 = ( xmax + xmin ) / 2.
    for n in range(len(nodes)):
        nodes[n][0] -= x0
    #
    # center the object in z=streamwise
    # mid chord estimated at z=-550.
    z0 = -550.
    for n in range(len(nodes)):
        nodes[n][2] -= z0
    #
    # now that object is centered in x,z:
    # rotate about y-axis so x=streamwise
    theta =  90. * math.acos(-1.) / 180.  # radians
    c = math.cos(theta)
    s = math.sin(theta)
    rot = [ [  c, 0, s ],
            [  0, 1, 0 ],
            [ -s, 0, c ] ]
    for n in range(len(nodes)):
        node = nodes[n]
        xyz  = node
        xyz_rot = [0,0,0]
        for jrow in range(3):
            for icol in range(3):
                xyz_rot[jrow] += fac * ( rot[jrow][icol] * xyz[icol] )
        nodes[n] = xyz_rot
    #
    # return modified nodes
    return nodes

if __name__ == "__main__":
    filename_obj= 'bell_x1_mesh.obj'
    filename_tec= 'bell_x1_mesh.tec'
    #
    # read the surface
    nodes, faces = read_obj_file(filename_obj)
    #
    # center the surface
    nodes = center_obj(nodes)
    #
    # write the surface
    write_tec_file(filename_tec,nodes,faces)




