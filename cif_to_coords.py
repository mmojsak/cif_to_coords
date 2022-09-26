#cif_path = input("Please give the path to your .cif: ")

#print(f"The path to your cif is: {cif_path}")

from gemmi import cif as cif_parser
import numpy as np
import math

def cert_float(string_with_uncertainty):
    """Converts a string representing a value with uncertainty to a 'certain' float."""
    return float(string_with_uncertainty.split('(')[0])

# read and parse a CIF file into a variable called cif
cif = cif_parser.read_file('test_AuCl4_L.cif').sole_block() # remove .sole_block() if file has more blocks
#print(cif.name)

# read unit cell dimensions into float variables
a = cert_float(cif.find_value('_cell_length_a'))
b = cert_float(cif.find_value('_cell_length_b'))
c = cert_float(cif.find_value('_cell_length_c'))

# read unit cell angles into float variables
alpha = cert_float(cif.find_value('_cell_angle_alpha'))
beta = cert_float(cif.find_value('_cell_angle_beta'))
gamma = cert_float(cif.find_value('_cell_angle_gamma'))
# read labels and types of all atoms
atom_label_list = list(cif.find_loop('_atom_site_label'))
atom_type_list = list(cif.find_loop('_atom_site_type_symbol'))

# read lists of xyz coordinates of all atoms
atom_fract_coord_x_list = list(cif.find_loop('_atom_site_fract_x'))
atom_fract_coord_y_list = list(cif.find_loop('_atom_site_fract_y'))
atom_fract_coord_z_list = list(cif.find_loop('_atom_site_fract_z'))

# convert all read coordinates into floats using the pre-defined cert_float function
for (index, coord) in enumerate(atom_fract_coord_x_list):
    coord_float = cert_float(coord)
    atom_fract_coord_x_list[index] = coord_float
    
for (index, coord) in enumerate(atom_fract_coord_y_list):
    coord_float = cert_float(coord)
    atom_fract_coord_y_list[index] = coord_float
    
for (index, coord) in enumerate(atom_fract_coord_z_list):
    coord_float = cert_float(coord)
    atom_fract_coord_z_list[index] = coord_float
    
# defining the orthogonalisation matrix
a_star = math.acos((math.cos(beta) * math.cos(gamma) - math.cos(alpha)) / (math.sin(beta)) * math.sin(gamma))
ax = a
ay = 0
az = 0
bx = b * math.cos(gamma)
by = b * math.sin(gamma)
bz = 0
cx = c * math.cos(beta)
#cy = c * n
#cz = c * math.sqrt(math.sin(beta)**2 - n**2)
#orth_matrix = np.array([np.array([ax,ay,az]),np.array([bx,by,bz]),np.array([cx,cy,cz])])
    
print(a_star)


class Atom(object):
    
    def __init__(self, label='Unlabelled', atom_type='Unknown', frac_xyz=np.array([0,0,0]), 
                 cart_xyz=np.array([0,0,0])):
        self.label = label
        self.atom_type = atom_type
        self.frac_xyz = frac_xyz
        self.cart_xyz = cart_xyz
        self.frac_x = frac_xyz[0]
        self.frac_y = frac_xyz[1]
        self.frac_z = frac_xyz[2]
        self.cart_x = cart_xyz[0]
        self.cart_y = cart_xyz[1]
        self.cart_z = cart_xyz[2]

    def __str__(self):
        return f"A CelestialBody object representing {self.name}."
 
atoms = list()    
for i in range(0,len(atom_label_list)):
    atom = Atom(label=atom_label_list[i], atom_type=atom_type_list[i],
                frac_xyz=np.array([atom_fract_coord_x_list[i], atom_fract_coord_y_list[i], atom_fract_coord_z_list[i]]))
    atoms.append(atom)
    
#for atom in atoms:
#    print(atom.frac_xyz)
    
    
    
    
#c = {"one": 1, "two": 2}
#for k, v in c.items():
#    exec(f"{k} = {v}") 

