# -*- coding: utf-8 -*-
"""
Created on Thu Jul 12 22:12:08 2018

@author: Yongjin
command
$ python Structural_descriptor.py [filenames] [entry or entries]

[filename]: name of vasprun.xml file or POSCAR file (wildcard (*,?) is usable!)
    if you use xml file, you can extract band gap despite longer execution time
[entry or entries]: Descripor of user choice
    - a, b, c, V: lattice parameters and volume
    - atom1-atom2: distance between two atoms are called with one dash (-)
        atom can be designated as VESTA format or index number starting from 1.
        Ex) In SrTiO3, Sr --> Sr1 or 1, Fe --> Fe1 or 2, O1 --> O1 or 3
    - Fe5-O18-Mn8: bond angle is callable with two dashes joining three atoms
        Note: This feature is to calcaulte bond angle so it automatically
        copies the atom so that 'bond angle' is calculated.
        Angle between 'far' apart atoms are not reliable with this script
    - SG, SG#: Space group with symbol or space group number
    - Eg: Band gap (xml file is needed for this)
    - Mathmatical operation of descriptors:
        If you format descriptor or column with bracket ('[' and ']'), 
        you can execute simple operations.
        Ex) [Fe5-O1]+[Fe5-O2]: sum of bond length of Fe5-O1 and Fe5-O2
        Ex) [11-O1]+[Fe5-O2]: sum of bond length of (11th atom)-O1 and Fe5-O2
        Ex) [Fe5-O1]/[Fe5-O2]: ratio of bond length between Fe5-O1 and Fe5-O2
        Ex) [1]/[2]: ratio between 1th column value and 2nd column value
            (Note that 0th column is for filename).
        
Features to implement:
----------------- Done-----------------------------
    1. Make function to get angle (based on shortest vectors under periodic boundary condition)
    2. Digest system argument for bonding (regular expression?)
    3. Important: How to recognize '-', and '--' with evaluation case.
    4. How to print. print labels first. Need to cover 4.
    5. Printing format: make function to implement .E, based on order
    6. Cover variety (num: atom index from 1-n, atom_index=Fe1)
    7. Cover multiple selection of target structure --> glob
    8. Implementation to digest CONTCAR or POSCAR files
    9. Space group is supported
-----------------To be Done-----------------------------
    19. Cover multiple statement of argument (ex. M-O bonds, M-O-M angles, searching method is needed) --> YDNL
    8. How to extract Magmom --> Using outcar, in separate script
"""
from pymatgen.io.vasp.outputs import Vasprun, Poscar
from pymatgen import Structure
from pymatgen.util.coord import pbc_shortest_vectors, get_angle
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.periodic_table import Element
from MTD_functions import *
from math import pi
import sys
import re
import glob
import numpy as np
from numpy import sqrt
import os

#start=timeit.default_timer()
if len(sys.argv) ==1:
    #sys.argv=['script_name',"*.vasp",'SG','GII','Mn5-O18', 'Mn5-O1', '13-17' ,'Mn5-O18-Mn8','a', '[Mn5-O18]+[Mn5-O19]', '([-1]-[-2])', 'V']
    sys.argv=['script_name',"*.vasp",'GII', 'Ewald']
    sys.argv=['script_name',"SMO*.vasp",'Mn1-O11', 'Mn1-O5-Mn4']


##### Part 1: Preparation. File lists and print label line  ###################
## Collect list of files with glob
#file_list=glob.glob(sys.argv[1])
#file_list.sort()
## Construct list_entries, which will be used for label of result
## sys.argv[0] is name of script,
## and sys.argv[1] is list of file names, so rename the labe as 'Filename'
#list_entries=sys.argv[1::]
#list_entries[0]='Filename'

file_list=list()
list_entries=list(['Filename'])
for entry in sys.argv[1::]:
    # if the entry is the file name, append it to file_list
    # otherwise, append it to list of entries.
    if os.path.isfile(entry):
        file_list.append(entry)
    elif '*' in entry or '?' in entry:
        temp_list=glob.glob(entry)
        temp_list.sort()
        file_list=file_list+temp_list
    else:
        list_entries.append(entry)

# Printing label
# formatted for better outlook.
# Position is determined by maximum length of file names
print(list_entries[0].ljust(len(max(file_list,key=len))+2)+\
      '  '.join('%-6s' % entry for entry in list_entries[1::]))

# list_errors is to collect error messages for user
list_errors=[]
###############################################################################


##### Part 2: Iteration for each file. ########################################
for filename in file_list:
    
    #### Part 2-1: determine the type of file. xml or POSCAR ####
    if filename[-4::]=='.xml':
        ## if file is xml
        dos_vrun=Vasprun(filename)
        total_dos = dos_vrun.complete_dos
        #get structure from vasprun
    #    struct = dos_vrun.structures[-1]
        struct = dos_vrun.final_structure
    else:
        # if file is POSCAR format
        P=Poscar.from_file(filename)
        struct = P.structure
    #############################################################
    
    #### Part 2-2: Make dictionary: from label to site index ####
    #For easy use, Labe them in conventional manner.
    #labeled as Element + site position in poscar 
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
        label2site_index.update({'{0}{1}'.format(struct.species[i], n_atom_count_dict[struct[i].specie]):i})
    #############################################################
    

    ### Part 2-3: Define function 'obtain_descriptor'. ##########
    # This part is edited when new descriptor is added 
    def obtain_descriptor(str_entry):
        """
        input:
            str_entry -- str, argument from command line
            ex) Fe1-O18, Fe1-O15-Fe8, a, V, Eg, etc
        output: current_value -- str or float,
        'N/A' if there is no match,
        str or float if there is.
        
        """
        
        try:
            #option: atomic features
#            if 
            
            #option: Bond with one '-'
            if len(re.findall('-',str_entry))==1:
                atom1_str,atom2_str=re.split('-',str_entry)
                current_value=Bond_length(struct, label2site_index,atom1_str,atom2_str)
            
            #option: Bond with two '-'
            elif len(re.findall('-',str_entry))==2:
                atom1_str,atom2_str,atom3_str=re.split('-',str_entry)
                current_value=Bond_angle(struct, label2site_index,atom1_str,atom2_str,atom3_str)
            
            #option: lattice parameters
            elif str_entry == 'a':
                current_value=struct.lattice.a
            elif str_entry == 'b':
                current_value=struct.lattice.b
            elif str_entry == 'c':
                current_value=struct.lattice.c
            elif str_entry == 'alpha':
                current_value=struct.lattice.alpha
            elif str_entry == 'beta':
                current_value=struct.lattice.beta
            elif str_entry == 'gamma':
                current_value=struct.lattice.gamma
            elif str_entry == 'V':
                current_value=struct.lattice.volume
            elif str_entry == 'GII':
                # Maybe insert a line to extract formal valcne list from CONTCAR
                current_value=gii_compute(struct)  ################ To be done
            elif str_entry == 'Ewald':
                # Maybe insert a line to extract formal valcne list from CONTCAR
                current_value=Calc_Ewald(struct)  ################ To be done
            #option: spacegroup
            elif str_entry == 'SG':
                current_value=SpacegroupAnalyzer(struct).get_space_group_symbol()
            elif str_entry == 'SG#':
                current_value=SpacegroupAnalyzer(struct).get_space_group_number()
            #option: Band gap
            elif str_entry == 'Eg':
                current_value=total_dos.get_gap()
                
            # the number of atoms in the lattice
            elif str_entry == 'natom':
                current_value=struct.composition.num_atoms
            
            # If there is no matching key
            else:
                current_value='N/A'
            
        except KeyError:
            print(str_entry)
            # Error message
            # This often happens when user tried incorrect atom label
            list_errors.append('Note: Check atom labels')
            current_value='N/A'
        
        return current_value
    #############################################################
    
    #### Part 2-4: Collect descriptors of current file ##########
    # list_values
    # This will be the list of data for current file
    list_values=[]
    list_values.append(filename)
    
    # Iteration for each argument(entry)
    for i in range(2,len(sys.argv)):
        
        ## Case 1: Single component in one entry.
        # This does not contain '[' or ']' in the argument
        if '[' not in sys.argv[i]:
            current_value=obtain_descriptor(sys.argv[i])
            
        ## Case 2: Simple operations base on eval function
        # Multiple component in one entry
        # When entry contains '[' or ']', eval function is used to calculate.
        # regular expression is used to seperate each component
        # re.findall('\[[\w\-\_]+\]',current_entry)
        # other characters are to be included when needed
        else:  
            current_entry=sys.argv[i]
            # component is a singly callable value as descriptor
            list_components=re.findall('\[[\w\-\_]+\]',current_entry)
            
            # iteration for each componenet
            for component in list_components:
                # component[1:-1] is to extract strings between brackets
                
                # Component case1: numeric and calls value from list_values
#                elif component[1:-1].isnumeric():
                if len(re.findall('[^0-9:-]',component[1:-1]))==0:
                    # This is when bracket calls label value
                    try:
                        comp_value=list_values[int(component[1:-1])]
                    except (TypeError, ValueError,IndexError) as e:
                        #Error message
                        list_errors.append('Note: Check the column indicies in the brackets')
#                        print ("There is a syntax problem within bracket")
                        break
                
                # Component case2: able to call descriptor function
                elif obtain_descriptor(component[1:-1]) != 'N/A':
                    comp_value=obtain_descriptor(component[1:-1])
                
                # Component case3: comp_value is N/A and it is not index
                else:
                    # Then the value of formula is not reliable, thus print 'N/A'
                    break
                current_entry=current_entry.replace(component,str(comp_value))
            
            try:
#                print(current_entry)
                current_value=eval(current_entry)
            except (ZeroDivisionError, ValueError, NameError, TypeError) as e:
                #Error message
                list_errors.append('Note: Check the value of called column')
                current_value='N/A'
            
        list_values.append(current_value)
    #############################################################
    
    
    #### Part 2-5: Printing values
    # Set starting position first, for better outlook
    str_value_line=(list_values[0].ljust(len(max(file_list,key=len))+2))
    # iterate for each value, depending on value type
    for i, value in enumerate(list_values[1::]):
        # value case 1: string
        if type(value) ==str:
            str_value_line=str_value_line+value.ljust(max(8,len(list_entries[i+1])+2)) 
        
        # value case 2: number
        # print_format is used for better outlook
        else:
            str_value_line=str_value_line+print_format(value).ljust(max(8,len(list_entries[i+1])+2))
    print(str_value_line)
    #############################################################
###############################################################################


##### Part 3: Print error messages  ###########################################
# Omit repeated error messages
list_errors=list(set(list_errors))
for error in list_errors:
    print(error)
###############################################################################
