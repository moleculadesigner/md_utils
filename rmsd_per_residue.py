"""
This module allows to calculate residue fluctations in molecular dynamics trajectory.
"""
import re
import sys
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
ATOM   1899  HD23LEU A 121      40.604  62.341   8.996        0.10000
    ```
    """
    atom_re = re.compile(r"^ATOM.+$", re.M)
    for path in pdbs:
        """
        print(path)
        """
        with open(path, 'r') as pdb_frame:
            print(path)
            content = pdb_frame.read()
            atoms = atom_re.findall(content)
            for entry in atoms:
                print(int(entry[23:26]))

def main():
    print(__doc__)
    get_trajectory(*sys.argv[1:])

    

if __name__ == '__main__':
    main()