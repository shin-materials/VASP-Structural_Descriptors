# -*- coding: utf-8 -*-
"""
Created on Thu Sep  6 19:51:52 2018

@author: Yongjin
"""

from pymatgen.io.vasp.outputs import Vasprun, Poscar
from pymatgen.core.periodic_table import Element
#from pymatgen.electronic_structure.core import Spin
#from pymatgen.core.structure import Structure
from pymatgen import Structure
from pymatgen.util.coord import pbc_shortest_vectors, get_angle
from math import log10, floor
import os
import sys
import re
import numpy as np
import pandas as pd
from pymatgen.analysis.bond_valence import calculate_bv_sum, calculate_bv_sum_unordered, BVAnalyzer
from pymatgen.analysis.ewald import EwaldSummation

bv_analyzer=BVAnalyzer(max_radius=4)


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

# def GII_calculation(pmg_struct):
#     """
#     input: pmg_struct: pymatgen structure
#     output: GII calculation
#     """


#     return GII_val


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
    if abs(number)==0:
        str_form=format(number,'-.5g')
    elif floor(log10(abs(number)))>=4:
        str_form=format(number,'.3E')
    elif floor(log10(abs(number)))<=-1 and floor(log10(abs(number)))>=-3:
        str_form=format(number,'-5.5f')
    elif floor(log10(abs(number)))<-3:
        str_form=format(number,'.3E')
    else:
        str_form=format(number,'-.5g')
    return str_form

def Calc_Ewald(pmg_struct, formal_val=[]):
    """
    input: pmg_struct: pymatgen structure
           formal_val: list - list of valence for each atom

    TBD: use input formal valence list to decorate the structure
    """
    # GET valence list
    if len(formal_val)==0:
        bv_analyzer=BVAnalyzer(max_radius=4) #max_radius default is 4
        formal_val=bv_analyzer.get_valences(pmg_struct)

    # default for cutoff is set to None
    real_cutoff = None
    rec_cutoff = None

    # Oxidation states are decorated automatically.  
    decorated_struct = bv_analyzer.get_oxi_state_decorated_structure(pmg_struct)
    NUM_sites=pmg_struct.num_sites
    
    # Per atom
    ewald_per_atom=1/NUM_sites*EwaldSummation(decorated_struct,
        real_space_cut=real_cutoff, recip_space_cut=rec_cutoff).total_energy
    return ewald_per_atom


####################################################################
##### GII calculation functions: Modified from Nick's version ######
####################################################################

# Read Bond valence dataset
fileDir = os.path.dirname(os.path.realpath(sys.argv[0]))
bv_data = pd.read_csv(os.path.join(fileDir,"Bond_valences2016.csv"))
# Use element names and valences to lookup bond valence
def get_bv_params(cation, anion, cat_val, an_val):
    bond_val_list = bv_data[(bv_data['Atom1'] == cation) & (bv_data['Atom1_valence'] == cat_val)\
                  & (bv_data['Atom2'] == anion) & (bv_data['Atom2_valence'] == an_val)]
    return bond_val_list.iloc[0] # If multiple values exist, take first one

def gii_compute(pmg_struct, formal_val=[]):
    """
    input: formal_val: list - list of valence for each atom
    """
    # Define length cutoff
    cutoff = 3.5

    # GET valence list
    if len(formal_val)==0:
        # bv_analyzer=BVAnalyzer(max_radius=4) #max_radius default is 4
        formal_val=bv_analyzer.get_valences(pmg_struct)

    # loop for every site pair
    bv_diffs=[]
    for i,i_site in enumerate(pmg_struct.sites):
        bv_atom=0
        for j,j_site in enumerate(pmg_struct.sites):
            site_distance=pmg_struct.get_distance(i,j)
            #conditions: exclude itself, only pairs within cufoff range, and only cation-anion pairs 
            if site_distance > 0 \
            and site_distance <=cutoff \
            and formal_val[i]*formal_val[j] < 0:
                params_df=bv_data[(bv_data['Atom1'] == str(i_site.specie)) & (bv_data['Atom1_valence'] == formal_val[i])\
                 & (bv_data['Atom2'] == str(j_site.specie)) & (bv_data['Atom2_valence'] == formal_val[j])]
                # I ignore the case of not having matched item
                # Possible alternative can be made by adopting value from slightly different formal valence
                if params_df.empty:
                    bv_atom += 0
                else:
                    params=params_df.iloc[0]
                    bv_atom += np.exp((params['Ro']- site_distance)/params['B'])
        if bv_atom != 0:
            bv_diffs.append(np.power(formal_val[i] - bv_atom, 2))
    GII_val = np.sqrt(np.sum(bv_diffs)/pmg_struct.composition.num_atoms)
    return GII_val





# A function to calculate a generalized Goldschmidt tolerance factor for perovskites and RP phases
def calc_tol_factor(ion_list, valence_list, rp=0):
    if len(ion_list) > 4 or len(ion_list) < 3:
        print("Error: there should be three or four elements")
        return None
    if len(ion_list) < 4:
        for i in range(len(valence_list)): # If charge is 2-, make -2 to match tables
            if valence_list[i][-1] == '-':
                valence_list[i] = valence_list[i][-1] + valence_list[i][:-1]
        for i in range(len(valence_list)): # Similarly, change 2+ to 2
            valence_list[i] = int(valence_list[i].strip("+"))
        
    if len(ion_list) == 4:
#         print("RED ALERT: We are taking averages of bond valence parameters")
        AO_value1 = get_bv_params(ion_list[0], ion_list[-1], valence_list[0], valence_list[-1])
        AO_value2 = get_bv_params(ion_list[1], ion_list[-1], valence_list[1], valence_list[-1])
        AO_values = np.concatenate([AO_value1.values.reshape(1, len(AO_value1)), 
                                    AO_value2.values.reshape(1, len(AO_value2))])
        AO_B = np.average(AO_values[:, 4])
        AO_Ro = np.average(AO_values[:, 5])
        AO_valence = np.average(AO_values[:, 1]) # RED ALERT: We are taking averages of bond valence parameters
    else:
        AO_row = get_bv_params(ion_list[0], ion_list[-1], valence_list[0], valence_list[-1])
    
    BO_row = get_bv_params(ion_list[-2], ion_list[-1], valence_list[-2], valence_list[-1])
    
    
    if len(ion_list) != 4:
        if rp == 0:
            AO_bv = AO_row['Ro']-AO_row['B'] * np.log(AO_row['Atom1_valence']/12)
            BO_bv = BO_row['Ro']-BO_row['B'] * np.log(BO_row['Atom1_valence']/6)               
        else: # Currently for Ruddlesden-Popper phases a naive weighted sum is used between A-site coordination of 
              # 9 in the rocksalt layer and 12 in perovskite
            AO_bv = AO_row['Ro']-AO_row['B'] * np.log(AO_row['Atom1_valence']/((9+12*(rp-1))/rp))
            BO_bv = BO_row['Ro']-BO_row['B'] * np.log(BO_row['Atom1_valence']/6)
    else:
        if rp == 0:
            AO_bv = AO_Ro-AO_B * np.log(AO_valence/12)
            BO_bv = BO_row['Ro']-BO_row['B'] * np.log(BO_row['Atom1_valence']/6)               
        else: # Currently for Ruddlesden-Popper phases a naive weighted sum is used between A-site coordination of 
              # 9 in the rocksalt layer and 12 in perovskite
            AO_bv = AO_Ro-AO_B * np.log(AO_valence/((9+12*(rp-1))/rp))
            BO_bv = BO_row['Ro']-BO_row['B'] * np.log(BO_row['Atom1_valence']/6)
    
    tol_fact = AO_bv / (2**0.5 * BO_bv)
    
    return tol_fact


# def MObonds_greedy(structure,Msite, cutoff=3.0):
def MObonds_identifier(structure,Msite, cutoff=3.0):
    '''
    This function takes a pymatgen structure and perovskite Bsite and returns 
    a list of the bond lengths associated with the Bsite/oxygen bond lengths
    '''
    bond_lengths = []
    # determine Bsite and oxygen indexes
    for site in structure.sites:
        if Msite in str(site):
            neighbors = structure.get_neighbors(site, r = cutoff, include_index=True)
            for neighbor in neighbors:
                elems_on_neighsite = structure.species_and_occu[neighbor[2]].elements
                symbols = [elem.symbol for elem in elems_on_neighsite]
                if Msite in symbols:
                    continue
                else:
                    bond_lengths.append(neighbor[1])
            if not bond_lengths:
                neighbors = structure.get_neighbors(site, r = cutoff+0.6, include_index=True)
                for neighbor in neighbors:
                    elems_on_neighsite = structure.species_and_occu[neighbor[2]].elements
                    symbols = [elem.symbol for elem in elems_on_neighsite]
                    if Msite in symbols:
                        continue
                    else:
                        bond_lengths.append(neighbor[1])
                return bond_lengths
            else:
                return bond_lengths
    return bond_lengths








#     pymat_neighbors = struct.get_all_neighbors(cutoff, include_index=True)
#     if not formal_val:
#         print("Please specify formal valences of all species. Returning None")
#         return
    
#     # for loop to calculate the BV sum on each atom
#     bv = BVAnalyzer(max_radius=cutoff+0.1)
#     bv_diffs = []
#     for atom_indx, neighbors in enumerate(pymat_neighbors):
#         bv = 0
#         for pair in neighbors:
#             atom = struct.species[atom_indx].symbol
#             neighbor = struct.species[pair[2]].symbol
            
#             try:
#                 if iscation(atom) and isanion(neighbor):
#                     params = get_bv_params(cation=atom, anion=neighbor, 
#                                            cat_val=float(formal_val[atom]), an_val=float(formal_val[neighbor]))
#                     bv += np.exp((params['Ro']- pair[1])/params['B'])
#                 elif iscation(neighbor) and isanion(atom):
#                     params = get_bv_params(cation=neighbor, anion=atom, cat_val=float(formal_val[neighbor]), 
#                                            an_val=float(formal_val[atom]))
#                     bv += np.exp((params['Ro']- pair[1])/params['B'])
#             except:
#                 print("Trouble with atom: {} and neighbor: {} in {}".format(atom, neighbor, 
#                                                                             struct.formula))
#                 print('Looking for +/- 1 similar cation valence states in BV table')
#                 try:
#                     if iscation(atom) and isanion(neighbor):
#                         params = get_bv_params(cation=atom, anion=neighbor, 
#                                    cat_val=float(formal_val[atom])+1, an_val=float(formal_val[neighbor]))
#                         bv += np.exp((params['Ro']- pair[1])/params['B'])
#                     elif iscation(neighbor) and isanion(atom):
#                         params = get_bv_params(cation=neighbor, anion=atom, 
#                                    cat_val=float(formal_val[neighbor])+1, an_val=float(formal_val[atom]))
#                         bv += np.exp((params['Ro']- pair[1])/params['B'])
#                     formal_val[atom] = float(formal_val[atom])+1
#                 except:
#                     try:
#                         if iscation(atom) and isanion(neighbor):
#                             params = get_bv_params(cation=atom, anion=neighbor, 
#                                    cat_val=float(formal_val[atom])+1, an_val=float(formal_val[neighbor]))
#                             bv += np.exp((params['Ro']- pair[1])/params['B'])
#                         elif iscation(neighbor) and isanion(atom):
#                             params = get_bv_params(cation=neighbor, anion=atom, 
#                                    cat_val=float(formal_val[neighbor])+1, an_val=float(formal_val[atom]))
#                             bv += np.exp((params['Ro']- pair[1])/params['B'])
#                         formal_val[atom] = float(formal_val[atom])-1
#                     except:
#                         print("No similar valence states found. Returning None")
#                         return None
# #         print('Atom: {}, BV: {}'.format(struct.species[atom_indx].symbol, bv))

#         bv_diffs.append(np.power(abs(float(formal_val[struct.species[atom_indx].symbol])) - bv, 2))
# #         print('BV_diffs: {}'.format(bv_diffs))
        
#     GII_val = np.sqrt(np.sum(bv_diffs)/struct.composition.num_atoms)
#     return GII_val