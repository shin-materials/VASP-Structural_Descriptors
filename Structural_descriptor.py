# -*- coding: utf-8 -*-
"""
Created on Thu Jul 12 22:12:08 2018

@author: Yongjin
Features to implement:
    1. Make function to get angle (based on shortest vectors under periodic boundary condition)
    2. Digest system argument for bonding (regular expression?)
    3. Cover variety (num: atom index from 1-n, atom_index=Fe1)
    4. Cover multiple statement of argument (ex. M-O bonds, M-O-M angles, searching method is needed)
    5. Cover multiple selection of target structure
    6. Interface. Implementation to digest CONTCAR or POSCAR files
    7. How to print. print labels first. Need to cover 4.
    8. Will I make it print the pwd, and time?

System arguments: XML file name, output filename, Type of elements
"""
#from pymatgen.electronic_structure import dos
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.core.periodic_table import Element
from pymatgen.electronic_structure.core import Spin
from pymatgen.util.coord_utils import pbc_shortest_vectors, get_angle
import sys

#import os
#import numpy as np

###### System arguments: XML file name, output filename, Type of elements #################
#xml_filename=sys.argv[1]
#system_name=sys.argv[2]
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
xml_filename="SFO_P30_AFM_G.xml"
system_name="SFO_U0_G30"
list_element=['Sr','Fe','O']
list_sites=['Fe1']
list_bond=['Fe1-O18']

dos_vrun=Vasprun(xml_filename)
total_dos = dos_vrun.complete_dos

#get structure from vasprun
struct = dos_vrun.structures[-1]

#For easy use, Labe them in conventional manner. labeled as Element + site position in poscar (i.e. BaTiO3 --> Ba1, Ti2, O1, O2, O3)
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


#struct.sites[8].distance(struct.sites[2])
# Band gap
label1=total_dos.get_gap()
# l_ap
#label1=struct.sites[label2site_index['Fe5']].distance(struct.sites[label2site_index['O18']])
label2=struct.get_distance(label2site_index['Fe5'],label2site_index['O18'])
# l_bs_c
#label2=struct.sites[label2site_index['Fe5']].distance(struct.sites[label2site_index['O1']])
label3=struct.get_distance(label2site_index['Fe5'],label2site_index['O1'])
# l_bs_O19
#label3=struct.sites[label2site_index['Fe5']].distance(struct.sites[label2site_index['O19']])
label4=struct.get_distance(label2site_index['Fe5'],label2site_index['O19'])
# l_Bs_O15
#label4=struct.sites[label2site_index['Fe5']].distance(struct.sites[label2site_index['O15']])
label5=struct.get_distance(label2site_index['Fe5'],label2site_index['O15'])
# l_ave
label6=(label2+2*label3+label4+label5)/5
# l_bs_ave
label7=(2*label3+label4+label5)/4
# bond_ratio
label8=label2/label7
# a
label9=struct.lattice.a
# b
label10=struct.lattice.b
# c
label11=struct.lattice.c
#V
label12=struct.lattice.volume
print('{0}\t{1:.3f}\t{2:.4f}\t{3:.4f}\t{4:.4f}\t{5:.4f}\t{6:.4f}\t{7:.4f}\t{8:.4f}\t{9:.4f}\t{10:.4f}\t{11:.4f}\t{12:.4f}'.format(system_name,label1,label2,label3,label4,label5,label6,label7,label8,label9,label10,label11,label12))

vector1=pbc_shortest_vectors(struct.lattice, struct.sites[label2site['O18']].frac_coords,struct.sites[label2site['Fe5']].frac_coords)
latt=struct.lattice
frac1=struct.sites[label2site_index['O18']].frac_coords
frac2=struct.sites[label2site_index['Fe5']].frac_coords

vector2=pbc_shortest_vectors(struct.lattice, struct.sites[label2site['O18']].frac_coords,struct.sites[label2site['Fe8']].frac_coords)
print(get_angle(vector1[0][0],vector2[0][0]))

#### Define ISPIN #####
ispin=len(total_dos.densities)
spin_list=[Spin.up, Spin.down]


####### Get DOS of each element ############        
#element_dos=total_dos.get_element_dos()
# Usage: element_dos[Element["Fe"]].densities[Spin.up]

####### Print with p4vasp data format ###########

#### Print header, including band gap, and list of printed DOS ######
#out_file1=open(out_filename1,'w')
##out_file1.write('# Total(Up, Down), Sr(Up,Down), Fe(Up,Down), O(Up,Down), Fe1(Up,Down)\n')
#list_printed=list_element+list_sites
#if ispin==2:
#    out_file1.write('# BandGap: {0:.3f}eV, Label: Spin up/down for Total, '.format(total_dos.get_gap())+', '.join(list_printed)+'\n')
#else:
#    out_file1.write('# BandGap: {0:.3f}eV, Label: Spin up for Total, '.format(total_dos.get_gap())+', '.join(list_printed)+'\n')
#  
#                  
## Energy grid point is obtained from site_dos of the first atom
#n_E_grid=len(total_dos.get_site_dos(struct[0]).energies)
#
##Printing Total_up
#for i in range(n_E_grid):
#    out_file1.write('{0:.4f}\t{1:.4f}\n'.format(total_dos.energies[i]-total_dos.efermi,total_dos.densities[Spin.up][i]))
##p4vasp format: Spacer between different DOS
#out_file1.write('\n')
#
##Printing Total_down
#if ispin == 2:    
#    for i in range(n_E_grid):
#        out_file1.write('{0:.4f}\t{1:.4f}\n'.format(total_dos.energies[i]-total_dos.efermi,-1.0*total_dos.densities[Spin.down][i]))
#    #p4vasp format: Spacer between different DOS
#    out_file1.write('\n')
#
##Printing each Element
#Element_list=list_element
#for k_element in Element_list:
#    # Loop for each spin
#    for j in range(ispin):
#        #if ispin==2, this will go through Spin.down
#        #j_spin is either Spin.up or Spin.down
#        j_spin=spin_list[j]
#        # Loop for each energy grid
#        for i in range(n_E_grid):
#            out_file1.write('{0:.4f}\t{1:.4f}\n'.format(total_dos.energies[i]-total_dos.efermi,float(int(j_spin))*element_dos[Element[k_element]].densities[j_spin][i]))
#        #p4vasp format: Spacer between different array
#        out_file1.write('\n')
#
#

#out_file1.close()
#print ('Interpolated Band gap is: ')
#if ispin==1:
#    print('Only up-spin is printed for following components')
#else:
#    print('Spin up/down is printed for following components')
#print ('--Elements: '+str(list_element))
#print ('--Atoms: '+str(list_sites))

