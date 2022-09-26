#cif_path = input("Please give the path to your .cif: ")

#print(f"The path to your cif is: {cif_path}")

from gemmi import cif as cif_parser
import numpy as np
from numpy import sin, cos, radians

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
alpha = radians(cert_float(cif.find_value('_cell_angle_alpha')))
beta = radians(cert_float(cif.find_value('_cell_angle_beta')))
gamma = radians(cert_float(cif.find_value('_cell_angle_gamma')))

# read unit cell volume
vol = cert_float(cif.find_value('_cell_volume'))

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

#ax = a
#ay = 0
#az = 0
#bx = b * cos(gamma)
#by = b * sin(gamma)
#bz = 0
#cx = c * cos(beta)
#cy = (b * c - bx * cx) / by
#cz = sqrt(c**2 - cx**2 - cy**2)
#orth_matrix = np.array([np.array([ax,ay,az]),np.array([bx,by,bz]),np.array([cx,cy,cz])])
    

ax = a
ay = 0
az = 0
bx = b * cos(gamma)
by = b * sin(gamma)
bz = 0
cx = c * cos(beta)
cy = -c * sin(beta) * ((cos(beta)*cos(gamma)-cos(alpha))/(sin(beta)*sin(gamma)))
cz = vol / (a * b * sin(gamma))

matrix = np.array([np.array([ax,bx,cx]),np.array([ay,by,cy]),np.array([az,bz,cz])])
inv_matrix = np.linalg.inv(matrix)

class Atom(object):
    """An object representing an atom in the unit cell."""
    def __init__(self, label='Unlabelled', atype='Unknown', frac_xyz=np.array([0,0,0])):
        self.label = label
        self.atype = atype
        self.frac_xyz = frac_xyz
        #self.frac_x = frac_xyz[0]
        #self.frac_y = frac_xyz[1]
        #self.frac_z = frac_xyz[2]
        
        self.cart_xyz = frac_xyz.dot(matrix)
        
        #self.cart_x = self.cart_xyz[0]
        #self.cart_y = self.cart_xyz[1]
        #self.cart_z = self.cart_xyz[2]
        
    def cart_to_frac(self):
        """"Converts cartesian coordinates to fractional."""
        self.frac_xyz = self.cart_xyz.dot(inv_matrix)
 

    def __str__(self):
        return f"An Atom object representing {self.atype}, labelled {self.label}."
 
atoms = list()    
for i in range(0,len(atom_label_list)):
    atom = Atom(label=atom_label_list[i], atype=atom_type_list[i],
                frac_xyz=np.array([atom_fract_coord_x_list[i], atom_fract_coord_y_list[i], atom_fract_coord_z_list[i]]))
    atoms.append(atom)
    
#for atom in atoms:
#    print(atom.frac_xyz)
    
print(a, b, c)

gold = atoms[0]
print('Original fractional coords:')
print(gold.frac_xyz)
print('Calculated cartesian coords:')
print(gold.cart_xyz)

gold.cart_to_frac()
print('Recalculated fractional coords:')
print(gold.frac_xyz)
    
    
#c = {"one": 1, "two": 2}
#for k, v in c.items():
#    exec(f"{k} = {v}") 

