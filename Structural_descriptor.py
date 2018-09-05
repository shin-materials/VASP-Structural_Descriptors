# -*- coding: utf-8 -*-
"""
Created on Thu Jul 12 22:12:08 2018

@author: Yongjin
Features to implement:
    1. Make function to get angle (based on shortest vectors under periodic boundary condition)
    2. Digest system argument for bonding (regular expression?)
    6. Implementation to digest CONTCAR or POSCAR files
-----------------Done-----------------------------
    3. Cover variety (num: atom index from 1-n, atom_index=Fe1)
    4. Cover multiple statement of argument (ex. M-O bonds, M-O-M angles, searching method is needed)
    5. Cover multiple selection of target structure
    7. How to print. print labels first. Need to cover 4.

System arguments: XML file name, output filename, Type of elements
"""
#from pymatgen.electronic_structure import dos
from pymatgen.io.vasp.outputs import Vasprun, Poscar
#from pymatgen.core.periodic_table import Element
#from pymatgen.electronic_structure.core import Spin
#from pymatgen.core.structure import Structure
from pymatgen import Structure
from pymatgen.util.coord import pbc_shortest_vectors, get_angle
import sys
import re

#import os
#import numpy as np


#sys.argv=['script_name',"SFO_P30_AFM_G.xml","SFO_U0_G30",'Fe5-O18', 'Fe5-O1','Fe5-O18-Fe8','a']

###### System arguments: XML file name, output filename, Type of elements #################
xml_filename=sys.argv[1]
system_name=sys.argv[2]
#
#list_element=[]
#list_sites=[]
#for i in range(3,len(sys.argv)):
#    if sys.argv[i][-1].isdigit():
#        list_sites.append(sys.argv[i])
#    else:
#        list_element.append(sys.argv[i])

##imports vasprun.xml into dos object
#dos_vrun = Vasprun("SFO_P30_AFM_G.xml")


dos_vrun=Vasprun(xml_filename)
total_dos = dos_vrun.complete_dos

#get structure from vasprun
struct = dos_vrun.structures[-1]
#P=Poscar.from_file('SMO_U0_E0.vasp')
#struct = P.struct

### Make dictionary: from label to site index ############
#For easy use, Labe them in conventional manner. labeled as Element + site position in poscar 
#(i.e. BaTiO3 --> Ba1, Ti2, O1, O2, O3)
n_atom_count_dict=dict()
label2site_index=dict()
for i in range(0,struct.num_sites):
    # Update label for each element
    if struct[i].specie in list(n_atom_count_dict.keys()):
        n_atom_count_dict.update({struct[i].specie:n_atom_count_dict[struct[i].specie]+1})
    else:
        n_atom_count_dict.update({struct[i].specie:1})
    #Example: BaTiO3 --> Ba1:0 Ti1:1 O1:2 O2:3 O3:4
    label2site_index.update({'{0}{1}'.format(struct.species[i],n_atom_count_dict[struct[i].specie]):i})
##########################################################

def Bond_angle(pmg_struct, label2site_index, atom1_label,atom2_label,atom3_label):
    '''
    input: atom1_label, atom2_label, atom3_label -- str,  Ex) Fe5, O1, etc. 
            label2site_index -- dictionary defined in script
            pmg_struct -- pymatgen structure
    output: angle -- float, Angle of (atom1-atom2-atom3) will be calculated.
        Angle is calculated by using shortest vector of (atom1-atom2) and 
        (atom3-atom2) within periodic boundary condition.
    '''
    # Get fractional coordinates for each atoms
    atom1_fcoords=pmg_struct.sites[label2site_index[atom1_label]].frac_coords
    atom2_fcoords=pmg_struct.sites[label2site_index[atom2_label]].frac_coords
    atom3_fcoords=pmg_struct.sites[label2site_index[atom3_label]].frac_coords
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
    input: atom1_label, atom2_label -- str,  Ex) Fe5, O1, etc. 
            label2site_index -- dictionary defined in script
            pmg_struct -- pymatgen structure
    output: bondlength -- float, bond length of (atom1-atom2) will be calculated.
    '''
    # Get fractional coordinates for each atoms
    atom1_site_index=label2site_index[atom1_label]
    atom2_site_index=label2site_index[atom2_label]
    bondlength=pmg_struct.get_distance(atom1_site_index,atom2_site_index)
    return bondlength


### Construct entries and options ############################
list_entries=[]
val=dict()
# Make list of entries
for i in range(3,len(sys.argv)):
    #option arguments
    
    #option: Bond with one '-'
    if len(re.findall('-',sys.argv[i]))==1:
        current_entry=sys.argv[i]
        atom1_str,atom2_str=re.split('-',current_entry)
        current_value=Bond_length(struct,label2site_index,atom1_str,atom2_str)
    
    #option: Bond with two '-'
    elif len(re.findall('-',sys.argv[i]))==2:
        current_entry=sys.argv[i]
        atom1_str,atom2_str,atom3_str=re.split('-',current_entry)
        current_value=Bond_angle(struct,label2site_index,atom1_str,atom2_str,atom3_str)
    
    #option: lattice parameters
    elif sys.argv[i] == 'a':
        current_value=struct.lattice.a
    elif sys.argv[i] == 'b':
        current_value=struct.lattice.b
    elif sys.argv[i] == 'c':
        current_value=struct.lattice.c
    elif sys.argv[i] == 'alpha':
        current_value=struct.lattice.alpha
    elif sys.argv[i] == 'beta':
        current_value=struct.lattice.beta
    elif sys.argv[i] == 'gamma':
        current_value=struct.lattice.gamma
    elif sys.argv[i] == 'V':
        current_value=struct.lattice.volume
    
    #option: Band gap
    elif sys.argv[i] == 'Eg':
        current_value=total_dos.get_gap()
    
    
    # supporting simple operations
    # 
    else:
        current_entry=sys.argv[i]
        current_entry=re.sub(r'\[',"val['",current_entry)
        current_entry=re.sub(r'\]',"']",current_entry)
        current_value=eval(current_entry)
        
#    val.append(current_value)
    val.update({sys.argv[i]:current_value})
    # normal entry arguments
#    else:
#        list_entries.append(sys.argv[i])
##############################################################
