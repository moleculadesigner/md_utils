"""
This module allows to calculate residue fluctations in molecular dynamics trajectory.
"""
import re
import sys
from operator import itemgetter
#from collections import namedtuple

import rmsd # pip install rmsd
import numpy as np # pip install numpy



def get_trajectory(*pdbs):
    """
    Goes over all .pdb files in *pdbs*, collect coordinates of atom per residue and returns dictionary with these coordinates.

    ```
    xyz = { 
        res_number : [[frame1_coord], [frame2_coord], ... ] 
        ...
        }
    ```
    Example `ATOM` entry:
    ```
                          23 26    32                      56
                           v v      v                      v
    ATOM   1899  HD23LEU A 121      40.604  62.341   8.996        0.10000
    ```
    """
    atom_re = re.compile(r"^ATOM.+$", re.M)
    float_re = re.compile(r"[-]?\d+[.]?\d*")
    total_coords = {}
    for path in pdbs:
        """
        print(path)
        """
        with open(path, 'r') as pdb_frame:
            content = pdb_frame.read()
        atoms = atom_re.findall(content)
        frame_coord = {}
        for entry in atoms:
            residue = int(entry[23:26])
            if not (residue in frame_coord.keys()):
                frame_coord[residue] = list()
            frame_coord[residue].append(list(map(float, float_re.findall(entry[32:56]))))
            
        for key in frame_coord.keys():
            if not (key in total_coords.keys()):
                total_coords[key] = list()
            res_array = np.array(frame_coord[key])
            total_coords[key].append(res_array)
                          
    return total_coords

def main():
    #print(__doc__)
    trajectory = get_trajectory(*sys.argv[1:])
    rmsds = {}
    for key in trajectory.keys():
        rmsds[key] = list()
        for i in range(len(trajectory[key])):
            rmsds[key].append(rmsd.kabsch_rmsd(trajectory[key][0], trajectory[key][i]))
    
    # avreages
    avrgs = {key : np.average(rmsds[key]) for key in rmsds.keys()}

    sorted_keys = list(map(itemgetter(0) ,sorted([item for item in avrgs.items()], reverse=True, key=itemgetter(1))))
    print('"","' + '","'.join(list(map(str, sorted_keys))) + '"')
    print('"average",{}'.format(','.join(list(map(str, [avrgs[key] for key in sorted_keys])))))

    for i in range(len(sys.argv[1:])):
        row = [str(i + 1)]
        for key in sorted_keys:
            row.append(str(rmsds[key][i]))
        print(','.join(row))
    """    
    """

    

if __name__ == '__main__':
    main()