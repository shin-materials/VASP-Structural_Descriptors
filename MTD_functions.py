# -*- coding: utf-8 -*-
"""
Created on Thu Sep  6 19:51:52 2018

@author: Yongjin
"""

from pymatgen.io.vasp.outputs import Vasprun, Poscar
#from pymatgen.core.periodic_table import Element
#from pymatgen.electronic_structure.core import Spin
#from pymatgen.core.structure import Structure
from pymatgen import Structure
from pymatgen.util.coord import pbc_shortest_vectors, get_angle
from math import log10, floor
import sys
import re




#########################################################################

def convert_site_index(label2site_index,str_atom):
    """
    input: label2site_index -- dictionary defined in script
        str_atom -- str, atom label like Fe1, or could be user index.
    output: pmg_site_index -- int, pymatgen index, starting from 0
    User index starts from 1, while python site index starts from 0.
    Based on the format, this function returns the pymatgen index
    """
    if str_atom.isnumeric():
        # if is numeric, python index is -1 from user index.
        # user index starts from 1, python index from 0
        pmg_site_index=int(str_atom)-1
    else:
        pmg_site_index=label2site_index[str_atom]
    return pmg_site_index


def Bond_angle(pmg_struct, label2site_index, atom1_label,atom2_label,atom3_label):
    '''
    input: atom1_label, atom2_label, atom3_label -- str,  Ex) Fe5, O1, etc. 
            label2site_index -- dictionary defined in script
            pmg_struct -- pymatgen structure
    output: angle -- float, Angle of (atom1-atom2-atom3) will be calculated.
        Angle is calculated by using shortest vector of (atom1-atom2) and 
        (atom3-atom2) within periodic boundary condition.
    '''
    atom1_pmg_index=convert_site_index(label2site_index, atom1_label)
    atom2_pmg_index=convert_site_index(label2site_index, atom2_label)
    atom3_pmg_index=convert_site_index(label2site_index, atom3_label)
        
    # Get fractional coordinates for each atoms
    atom1_fcoords=pmg_struct.sites[atom1_pmg_index].frac_coords
    atom2_fcoords=pmg_struct.sites[atom2_pmg_index].frac_coords
    atom3_fcoords=pmg_struct.sites[atom3_pmg_index].frac_coords
    # pbc_shortest_vector from atom2 to atom1
    vector1=pbc_shortest_vectors(pmg_struct.lattice,atom2_fcoords,atom1_fcoords)
    # pbc_shortest_vector from atom2 to atom3
    vector2=pbc_shortest_vectors(pmg_struct.lattice,atom2_fcoords,atom3_fcoords)
    
    # angle betseen two vectors
    # pbc_shortest_vector can digest set of sites,
    #but here we only use one 1 site to get angle.
    angle=get_angle(vector1[0][0],vector2[0][0])
    return angle

def Bond_length(pmg_struct, label2site_index, atom1_label,atom2_label):
    '''
    input: atom1_label, atom2_label -- str,  Ex) Fe5, O1, 1, 8 etc. 
            label2site_index -- dictionary defined in script
            pmg_struct -- pymatgen structure
    output: bondlength -- float, bond length of (atom1-atom2) will be calculated.
    '''
    atom1_pmg_index=convert_site_index(label2site_index, atom1_label)
    atom2_pmg_index=convert_site_index(label2site_index, atom2_label)
    bondlength=pmg_struct.get_distance(atom1_pmg_index,atom2_pmg_index)
    return bondlength

def ave(numbers):
    return float(sum(numbers)) / max(len(numbers), 1)

def print_format(number):
    if floor(log10(abs(number)))>=4:
        str_form=format(number,'.3E')
    elif floor(log10(abs(number)))<=-1 and floor(log10(abs(number)))>=-3:
        str_form=format(number,'-5.5f')
    elif floor(log10(abs(number)))<-3:
        str_form=format(number,'.3E')
    else:
        str_form=format(number,'-.5g')
    return str_form