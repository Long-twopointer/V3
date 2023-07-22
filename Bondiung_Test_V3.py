########################################################################################################################
#   Bonding_Test_V3.py ----- A code to build the .data file for a LAMMPS Simulation ; WIth new bonding Scheme
#
#   General Description - New and Improved Methods implemented and now\
#                           build new bonding scheme
#
#
########################################################################################################################

##      Import Extra Modules / Tools        ##

import numpy as np  # Numpy is used for Array Creation and Manipulation Mainly
from copy import deepcopy  # DeepCopy Is used to Copy Arrays


########################################################################################################################
##      Self Defined Functions      ##
def find_atom_id(x_find, y_find, z_find):
    for iV2, posV2 in enumerate(molecule_1):
        if (molecule_1[posV2].get('X') == x_find) \
                and (molecule_1[posV2].get('Y') == y_find) \
                and (molecule_1[posV2].get('Z') == z_find):
            return molecule_1[posV2].get('Atom_ID')


########################################################################################################################


###      Set Variables Used For Later       ###


#   Output File Options #

filename = input("Enter a name for the Output File (remember .data at the end!) :")  # String Name of File to be Created

#   Simulation Box Settings   #

sim_box_side_length = 100  # Each Side of the Simulation Box Will have sides of these lengths

SIM_SPACE = np.zeros((sim_box_side_length, sim_box_side_length, sim_box_side_length))
# Creates an Array Called SIM_SPACE (is a 3-D array (x,y,z of each simulation space position)

SIM_SPACE_PUREcopy = deepcopy(SIM_SPACE)  # Extra Copy of the Sim Space kept because its "pure"

x0 = int(np.floor(SIM_SPACE.shape[2] / 2))  # Finds the Center of the Sim Space ( x0, y0, z0 )
y0 = int(np.floor(SIM_SPACE.shape[1] / 2))
z0 = int(np.floor(SIM_SPACE.shape[0] / 2))

#   Molecule Geometry Creation Settings #

molecule_ids = 1  # LAMMPS Track molecules by their ID number ; this is the number of molecules created

molecule_1_thickness = 15  # The thickness is an extra range past the radius to include in the molecules

molecule_1_radius = 10  # Radius of each molecule

num_atoms_molecule_1 = 0  # Placeholder for total number of atoms created as part of the molecule

molecule_1_bond_length = 1  # Functions Similar to molecule radius but on the bond scale

molecule_1 = {}  # initializes the dictionary for molecule 1

molecule_1_COPY = {}  # Copy of Molecule dictionary

molecule_1_bonds = {}  # dictionary for bonds

num_molecule_1_bonds = 0

#   LAMMPS Data Inserts #

atom_types = 1  # LAMMPS Data Field ; Number of different types of Atoms

bond_types = 1  # LAMMPS Data Field ; Number of different types of Bonds

#   Program Performance Tracking

positions_checked = 0  # Keeps Count of How many points are iterated through ; easy way to check it's iterating correctly

########################################################################################################################

###     Creates 3-D Set of Points      ###

# Iterates Each Position of the sim (x,y,z

for x in range(0, sim_box_side_length):
    for y in range(0, sim_box_side_length):
        for z in range(0, sim_box_side_length):

            # This section measures the distance from the Center (x0,y0,z0) to each point and
            # deb: measures how far a coordinate in A is far from the center.
            # deb>=0: inside the sphere.
            # deb<0: outside the sphere

            deb = (molecule_1_radius ** 2) - ((x0 - x) ** 2) - ((y0 - y) ** 2) - ((z0 - z) ** 2)
            positions_checked = positions_checked + 1
            SIM_SPACE[x, y, z] = molecule_ids
            # if statement for deciding which points will be in MOLECULE 1

            if (deb >= 0) and (deb <= molecule_1_thickness):
                num_atoms_molecule_1 = num_atoms_molecule_1 + 1
                molecule_1[num_atoms_molecule_1] = {'Atom_ID': num_atoms_molecule_1,
                                                    'Molecule_ID': molecule_ids,
                                                    'Atom_Type': atom_types,
                                                    'X': x,
                                                    'Y': y,
                                                    'Z': z,
                                                    'Positions String': '{} {} {}'.format(x, y, z)}

                for x1 in range(x - molecule_1_bond_length, x + molecule_1_bond_length):
                    for y1 in range(y - molecule_1_bond_length, y + molecule_1_bond_length):
                        for z1 in range(z - molecule_1_bond_length, z + molecule_1_bond_length):
                            Atom_ID_2 = find_atom_id(x1, y1, z1)
                            if Atom_ID_2 is not None:
                                if Atom_ID_2 != num_atoms_molecule_1:
                                    num_molecule_1_bonds = num_molecule_1_bonds + 1
                                    molecule_1_bonds[num_molecule_1_bonds] = {'Bond_ID': num_molecule_1_bonds,
                                                                              'Bond_Type': bond_types,
                                                                              'atom1': num_atoms_molecule_1,
                                                                              'atom2': Atom_ID_2}

########################################################################################################################

###     Prints Sim Info To Make Sure Its Running Correct        ###

print("Number of Points in Simulation Space Iterated thorough: {} \n".format(positions_checked))

print("Number of Bonds Created: {} \n".format(molecule_1_bonds))

print("Sim Space Array: {} \n".format(SIM_SPACE))

########################################################################################################################

###     Write LAMMPS Data File      ###

with open(filename, 'w+') as fdata:  # opens a text file named a for the 'filename' variable
    fdata.write('{}\n\n'.format(filename))  # First line is a comment line

    ##     Header of Data File     ##

    fdata.write('{} atoms\n'.format(num_atoms_molecule_1))  # Specify number of atoms
    fdata.write('{} atom types\n'.format(atom_types))

    #   specify box dimensions      #

    fdata.write('{} {} xlo xhi\n'.format(0.0, sim_box_side_length))  # Writes X Position
    fdata.write('{} {} ylo yhi\n'.format(0.0, sim_box_side_length))  # Writes Y Position
    fdata.write('{} {} zlo zhi\n'.format(0.0, sim_box_side_length))  # Writes Z Position

    fdata.write('\n')  # Skips the Next line for Formatting

    ##      Data File Body      ##

    #       Atoms section       #
    fdata.write('Atoms\n\n')
    for i, pos in enumerate(molecule_1):
        fdata.write('{} {} {} {} \n'.format(molecule_1[pos].get('Atom_ID'),
                                            molecule_1[pos].get('Molecule_ID'),
                                            molecule_1[pos].get('Atom_Type'),
                                            molecule_1[pos].get('Positions String')))

    fdata.write('\n')

    #       Bonds Section       #
    fdata.write('Bonds\n\n')
    for iV2, posV2 in enumerate(molecule_1_bonds):
        fdata.write('{} {} {} {} \n'.format(molecule_1_bonds[posV2].get('Bond_ID'),
                                            molecule_1_bonds[posV2].get('Bond_Type'),
                                            molecule_1_bonds[posV2].get('atom1'),
                                            molecule_1_bonds[posV2].get('atom2')))

########################################################################################################################

###         Finished!!!     ###

print(' Data File Created Successfully!!! ; File name => {} '.format(filename))

########################################################################################################################
