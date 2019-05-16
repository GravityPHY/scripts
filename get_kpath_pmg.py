import sys
import os
##This script is written by Liz##
from pymatgen.core.structure import Structure
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


poscar_path = (sys.argv[1])

struc = Structure.from_file(poscar_path)
kpath = HighSymmKpath(struc)

file = open('kpath-mp','w')

spg = SpacegroupAnalyzer(struc).get_space_group_symbol()

file.write(spg + ' KPATH' + '\n')
file.write(str(16) + '\n')
file.write('Line-mode' + '\n')
file.write('reciprocal' + '\n')

kcoords = kpath.kpath['kpoints']
path = kpath.kpath['path']

for i in range(len(path)):
    last = len(path[i]) - 1
    for j in range(len(path[i])):
        coords = kcoords[path[i][j]]
        if j == 0:
            file.write(' ' + str(coords[0]) + ' ' + str(coords[1]) + ' ' + str(coords[2]) + ' ! ' + path[i][j] + '\n')
        elif (j > 0 and j < last):
            file.write(' ' + str(coords[0]) + ' ' + str(coords[1]) + ' ' + str(coords[2]) + ' ! ' + path[i][j] + '\n' + '\n')
            file.write(' ' + str(coords[0]) + ' ' + str(coords[1]) + ' ' + str(coords[2]) + ' ! ' + path[i][j] + '\n')
        elif j == last:
            file.write(' ' + str(coords[0]) + ' ' + str(coords[1]) + ' ' + str(coords[2]) + ' ! ' + path[i][j] + '\n' + '\n')

file.close()
