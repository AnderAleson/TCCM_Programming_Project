import argparse
import numpy as np
import os
import time

########################NOTe########################################

#THIS CODE CONVERGES WELL, THE THING IS STILL THE CONVERGENCE CRITERIA IS DIFFERENT FROM THE PROFESORRS.
#STILL I WILL UPLOAD TO DRIVE

#################### DEFINITION OF CLASSES  ################################

class Atom:
    def __init__(self,atm_num,x_coor,y_coor,z_coor,element):
        self.atm_num=int(atm_num)
        self.x_coor=float(x_coor)
        self.y_coor=float(y_coor)
        self.z_coor=float(z_coor)
        self.element= str(element)


class Bond:
    def __init__ (self, bond_number ,atom_1, atom_2, bond_type, bond_distance,bond_energy):
        self.bond_number=bond_number
        self.atom_1= int(atom_1)
        self.atom_2= int(atom_2)
        self.bond_type= bond_type
        self.bond_distance=bond_distance
        self.bond_energy=bond_energy

class Angle:
    def __init__ (self,atom_1, atom_2, atom_3, angle_value, angle_energy):
    
        self.atom_1=int(atom_1)
        self.atom_2=int(atom_2)
        self.atom_3=int(atom_3)
        self.angle_value=angle_value
        self.angle_energy=angle_energy

class Torsion:
    def __init__ (self, atom_1, atom_2, atom_3, atom_4, torsion_angle, torsion_energy):
        self.atom_1=atom_1
        self.atom_2=atom_2
        self.atom_3=atom_3
        self.atom_4=atom_4
        self.torsion_angle=torsion_angle
        self.torsion_energy=torsion_energy
    def __eq__(self, other):
        # Check for the same torsion, accounting for reversed atom order
        return (self.atom_1 == other.atom_1 and self.atom_2 == other.atom_2 and 
                self.atom_3 == other.atom_3 and self.atom_4 == other.atom_4) or \
               (self.atom_1 == other.atom_4 and self.atom_2 == other.atom_3 and 
                self.atom_3 == other.atom_2 and self.atom_4 == other.atom_1)


class VDW_interaction:
    def __init__ (self, atom_1, atom_2, distance, VDW_energy):
        self.atom_1=atom_1
        self.atom_2=atom_2
        self.distance=distance
        self.VDW_energy=VDW_energy






################  END DEFINITION OF CLASSES ################################

################# DEFINITION OF THE FUNCTIONS ##############################

def extract_file_data(file_path):
    with open(file_path, 'r') as file:
    
        lines=file.readlines()
        line_0=lines[0]
        words_in_line_0=line_0.split()
        number_of_atoms=int(words_in_line_0[0])
        number_of_bonds=int(words_in_line_0[1])
        number_of_different_bonds=words_in_line_0[2]
        number_of_type_of_bonds=words_in_line_0[3]
        lines_atom_type_and_coordinates=[]
        for i in range(number_of_atoms):
        
            lines_atom_type_and_coordinates.append(lines[i+1])
        lines_conectivity=lines[1+number_of_atoms:]
 
    #print("number of atoms: ",number_of_atoms)
    #print("number of bonds: ",number_of_bonds)
    #print("number of different bonds: ", number_of_different_bonds)
    #print("number of type of bonds: ",number_of_type_of_bonds)

    return lines_atom_type_and_coordinates, lines_conectivity, number_of_atoms, number_of_bonds


def get_atoms(lines_atom_data):
   
    atoms=[]
    for i in range(len(lines_atom_data)):
        atom_data = lines_atom_data[i].split()

        # Creating atom object and appending it to the list
        new_atom = Atom(atm_num=i+1, x_coor=atom_data[0], y_coor=atom_data[1], z_coor=atom_data[2], element=atom_data[3])
        atoms.append(new_atom)
    return atoms



def get_bonds(number_of_bonds,lines_conectivity):
    bonds=[]
    for i in range(number_of_bonds):
        bond_data=lines_conectivity[i].split()
        new_bond= Bond(bond_number=i+1, atom_1=bond_data[0], atom_2=bond_data[1], bond_type=bond_data[2], bond_distance=0, bond_energy=0)
        bonds.append(new_bond)


    return bonds

def compute_bond_energies(atoms,bond):
      #### here the -1 is fundamental since in pythone the index numbering in list starts from zerooooo !!!!!!!!!!!
    vec_atom_1_atom_2=np.array([atoms[bond.atom_1-1].x_coor-atoms[bond.atom_2-1].x_coor ,atoms[bond.atom_1-1].y_coor-atoms[bond.atom_2-1].y_coor , atoms[bond.atom_1-1].z_coor-atoms[bond.atom_2-1].z_coor])
    bond.bond_distance=np.linalg.norm(vec_atom_1_atom_2)
    if atoms[bond.atom_1-1].element == "C" and atoms[bond.atom_2-1].element == "C": #CALCULATING THE BOND ENERGIES FOR C-H BONDS 
        kb=300
        r0=1.53
        bond.bond_energy=kb*(bond.bond_distance-r0)**2
    if atoms[bond.atom_1-1].element == "C" and atoms[bond.atom_2-1].element == "H" or atoms[bond.atom_1-1].element == "H" and atoms[bond.atom_2-1].element == "C":  ##CALCULATING THE BOND ENERGY FOR C-H
        kb=350
        r0=1.11
        bond.bond_energy=kb*(bond.bond_distance-r0)**2
    return 

def get_angles(atoms,bonds):
    angles=[]
    for i in range(len(bonds)):            #Pufff... complex part better check further if every thing is working well...
        for j in range(i+1, len(bonds)):
            bond_1=bonds[i]
            bond_2=bonds[j]
            
            if bond_1.atom_1==bond_2.atom_1:
                new_angle=Angle(atom_1=bond_1.atom_2, atom_2=bond_1.atom_1, atom_3=bond_2.atom_2, angle_value=0, angle_energy=0)
                angles.append(new_angle)
            if bond_1.atom_2==bond_2.atom_1:
                new_angle=Angle(atom_1=bond_1.atom_1, atom_2=bond_1.atom_2, atom_3=bond_2.atom_2, angle_value=0, angle_energy=0)
                angles.append(new_angle)
            if bond_1.atom_1==bond_2.atom_2:
                new_angle=Angle(atom_1=bond_1.atom_2, atom_2=bond_1.atom_1, atom_3=bond_2.atom_1, angle_value=0, angle_energy=0)
                angles.append(new_angle)
            if bond_1.atom_2==bond_2.atom_2:
                new_angle=Angle(atom_1=bond_1.atom_1, atom_2=bond_1.atom_2, atom_3=bond_2.atom_1, angle_value=0, angle_energy=0)
                angles.append(new_angle)

    for angle in angles:

        #calculate the vector corresponding to atom2 from the angle to atom1. HERE THE -1 IS DUE TO THE POSITION IN THE LIST START FROM 0 !!!!!!!!!!!!!!!!!!!!!
        bond_1_vec=np.array([atoms[angle.atom_1-1].x_coor - atoms[angle.atom_2-1].x_coor,atoms[angle.atom_1-1].y_coor - atoms[angle.atom_2-1].y_coor, atoms[angle.atom_1-1].z_coor - atoms[angle.atom_2-1].z_coor])
        #Calculating the respecting vector modulus
        modulus_bond_1_vec=np.linalg.norm(bond_1_vec)

        #same for the other bond
        bond_2_vec=np.array([atoms[angle.atom_3-1].x_coor - atoms[angle.atom_2-1].x_coor,atoms[angle.atom_3-1].y_coor - atoms[angle.atom_2-1].y_coor, atoms[angle.atom_3-1].z_coor - atoms[angle.atom_2-1].z_coor])
        modulus_bond_2_vec=np.linalg.norm(bond_2_vec)


        #Calculating dot product between previous vectors.
        dot_product=np.dot(bond_1_vec,bond_2_vec)
        #Calculating angle energies THE POTENTIALS OR THE OPTIMAL ARE WRONG AND NEED TO BEEE CHANGEDDDD      
        angle_vec_1_vec_2=np.arccos(dot_product/(modulus_bond_1_vec*modulus_bond_2_vec))                     
        angle.angle_value=angle_vec_1_vec_2#*(180/np.pi) 


    return angles

def compute_angle_energies(atoms,angle):
    
    if atoms[angle.atom_1-1].element == "C" and atoms[angle.atom_2-1].element == "C" and atoms[angle.atom_3-1].element == "C":
        ka=60
        theta_0=109.5*(np.pi/180)
        angle.angle_energy=ka*(angle.angle_value-theta_0)**2
    if atoms[angle.atom_1-1].element == "H" and atoms[angle.atom_2-1].element == "C" and atoms[angle.atom_3-1].element == "C" or atoms[angle.atom_1-1].element == "C" and atoms[angle.atom_2-1].element == "H" and atoms[angle.atom_3-1].element == "C" or atoms[angle.atom_1-1].element == "C" and atoms[angle.atom_2-1].element == "C" and atoms[angle.atom_3-1].element == "H":
        ka=35
        theta_0=109.5*(np.pi/180)
        angle.angle_energy=ka*(angle.angle_value-theta_0)**2
    if atoms[angle.atom_1-1].element == "H" and atoms[angle.atom_2-1].element == "C" and atoms[angle.atom_3-1].element == "H":
        ka=35
        theta_0=109.5*(np.pi/180)
        angle.angle_energy=ka*(angle.angle_value-theta_0)**2
    return

def get_torsions(atoms,bonds,angles):
    torsions=[]
    for angle in angles:
    
        for bond in bonds:


            if angle.atom_1 == bond.atom_1 and angle.atom_2 != bond.atom_2:

                new_torsion = Torsion(atom_1=angle.atom_3, atom_2=angle.atom_2, atom_3=angle.atom_1, atom_4=bond.atom_2, torsion_angle=0, torsion_energy=0)

                if new_torsion not in torsions:
                
                    torsions.append(new_torsion)
            if angle.atom_1 == bond.atom_2 and angle.atom_2 != bond.atom_1:

                new_torsion = Torsion(atom_1=angle.atom_3, atom_2=angle.atom_2, atom_3=angle.atom_1, atom_4=bond.atom_1, torsion_angle=0, torsion_energy=0)

                if new_torsion not in torsions:
                
                    torsions.append(new_torsion)
            if angle.atom_3 == bond.atom_1 and angle.atom_2 != bond.atom_2:

                new_torsion = Torsion(atom_1=angle.atom_1, atom_2=angle.atom_2, atom_3=angle.atom_3, atom_4=bond.atom_2, torsion_angle=0, torsion_energy=0)

                if new_torsion not in torsions:
                
                    torsions.append(new_torsion)
            if angle.atom_3 == bond.atom_2 and angle.atom_2 != bond.atom_1:

                new_torsion = Torsion(atom_1=angle.atom_1, atom_2=angle.atom_2, atom_3=angle.atom_3, atom_4=bond.atom_1, torsion_angle=0, torsion_energy=0)

                if new_torsion not in torsions:
                
                    torsions.append(new_torsion)

    #print("number of torsions: ", len(torsions))


    for torsion in torsions:
        point_1 = np.array([atoms[torsion.atom_1-1].x_coor, atoms[torsion.atom_1-1].y_coor, atoms[torsion.atom_1-1].z_coor])
        point_2 = np.array([atoms[torsion.atom_2-1].x_coor, atoms[torsion.atom_2-1].y_coor, atoms[torsion.atom_2-1].z_coor])
        point_3 = np.array([atoms[torsion.atom_3-1].x_coor, atoms[torsion.atom_3-1].y_coor, atoms[torsion.atom_3-1].z_coor])
        point_4 = np.array([atoms[torsion.atom_4-1].x_coor, atoms[torsion.atom_4-1].y_coor, atoms[torsion.atom_4-1].z_coor])

        rAB = point_2 - point_1
        rBC = point_3 - point_2
        rCD = point_4 - point_3

        t = np.cross(rAB, rBC)
        u = np.cross(rBC, rCD)

        v = np.cross(t, u)

        rBC_unit = rBC / np.linalg.norm(rBC)

        cos_phi = np.dot(t, u) / (np.linalg.norm(t) * np.linalg.norm(u))
        sin_phi = np.dot(rBC_unit, v) / (np.linalg.norm(t) * np.linalg.norm(u))

        phi = np.arctan2(sin_phi, cos_phi)
        
        #VERY CAREFUL. FORUM: Yifan Jiang: value is larger than pi. If so, subtract or add it by 2pi
        


        torsion.torsion_angle = phi
    return torsions

def compute_torsion_energies(torsion):
    A=0.3
    torsion.torsion_energy=A*(1+np.cos(3*torsion.torsion_angle))
    return

def are_bonded(atom1, atom2, bonds,angles):
    for bond in bonds:         
                                   #When storing atom number in bonds, they were taken from the imput, Becareful because the numbers start from 1
                                                                                   #But when atoms are stored they are stored from 0 in the list
        if (bond.atom_1-1 == atom1 and bond.atom_2-1 == atom2) or (bond.atom_1-1 == atom2 and bond.atom_2-1 == atom1) :
                return True
        for angle in angles:
            if (angle.atom_1-1 == atom1 and angle.atom_3-1 == atom2) or (angle.atom_1-1 == atom2 and angle.atom_3-1 == atom1) :
                return True
    return False


def get_VDW_interactions(atoms,bonds,angles):
    VDW_interactions=[]
    for i in range(len(atoms)):
        for j in range(i+1, len(atoms)):

            if not are_bonded( i,j,bonds,angles):
                #print("VDW interaction detected between ",i+1,atoms[i].element,"and",j+1,atoms[j].element )
                new_VDW_interaction=VDW_interaction(atom_1=i, atom_2=j, distance=0, VDW_energy=0)
                VDW_interactions.append(new_VDW_interaction)

    return VDW_interactions


def compute_VDW_energy(VDW_interaction,atoms):
    
    distance_ij=np.sqrt((atoms[VDW_interaction.atom_1].x_coor-atoms[VDW_interaction.atom_2].x_coor)**2+(atoms[VDW_interaction.atom_1].y_coor-atoms[VDW_interaction.atom_2].y_coor)**2+(atoms[VDW_interaction.atom_1].z_coor-atoms[VDW_interaction.atom_2].z_coor)**2)

    if atoms[VDW_interaction.atom_1].element==atoms[VDW_interaction.atom_2].element=="H":
        A=4382.44
        B=22.932
    elif atoms[VDW_interaction.atom_1].element==atoms[VDW_interaction.atom_2].element=="C":
        A=946181.74
        B=514.714
    elif atoms[VDW_interaction.atom_1].element=="H" and atoms[VDW_interaction.atom_2].element=="C" or atoms[VDW_interaction.atom_1].element=="C" and atoms[VDW_interaction.atom_2].element=="H" :
        A=64393.99
        B=108.644

    Energy=(A/(distance_ij**(12)))-(B/(distance_ij**(6)))

    VDW_interaction.distance=distance_ij
    VDW_interaction.VDW_energy=Energy
    return

def compute_internal_coordinates(atoms,bonds,connectivity):
    #print("compute_internal_coordinates")
    internal_coordinates=[]

    bonds=get_bonds(len(bonds),connectivity)

    angles=get_angles(atoms,bonds)

    torsions=get_torsions(atoms,bonds,angles)


    for bond in bonds:
        compute_bond_energies(atoms,bond)
        internal_coordinates.append(bond.bond_distance)
  

  
    for angle in angles:
        compute_angle_energies(atoms,angle)
        internal_coordinates.append(angle.angle_value)



    for torsion in torsions:
        #print(torsion.torsion_angle)
        compute_torsion_energies(torsion)
        internal_coordinates.append(torsion.torsion_angle)

    return internal_coordinates



def Compute_tot_energy(atoms,bonds,angles,torsions,VDW_interactions,file_name):
    
    with open(file_name+".out", "a") as file:
        file.write("\nSYSTEM ENERGY DATA \n")
        total_bond_energy=0
        for bond in bonds:
            compute_bond_energies(atoms,bond)
            total_bond_energy=total_bond_energy+bond.bond_energy
        file.write(f"total bond enegy= {total_bond_energy} Kcal·mol⁻¹. \n")
        total_angle_energy=0
        for angle in angles:
            compute_angle_energies(atoms,angle)
            total_angle_energy=total_angle_energy+angle.angle_energy
        file.write(f"total_angle_energy {total_angle_energy} Kcal·mol⁻¹. \n")
        total_torsion_energy=0
        for torsion in torsions:
            compute_torsion_energies(torsion)
            total_torsion_energy=total_torsion_energy+torsion.torsion_energy
        file.write(f"total_torsion_energy {total_torsion_energy} Kcal·mol⁻¹. \n")
        total_VDW_Energy=0
        for VDW_interaction in VDW_interactions:
            compute_VDW_energy(VDW_interaction,atoms)
            total_VDW_Energy=total_VDW_Energy+VDW_interaction.VDW_energy
        file.write(f"total_VDW_Energy {total_VDW_Energy} Kcal·mol⁻¹. \n")

        total_energy=total_bond_energy+total_angle_energy+total_torsion_energy+total_VDW_Energy
        file.write(f"Total energy= {total_energy} Kcal·mol⁻¹. \n")
    file.close()

    return total_energy


def derib_bond(atoms,bond):
    derib=np.array([(atoms[bond.atom_1-1].x_coor-atoms[bond.atom_2-1].x_coor)/bond.bond_distance,(atoms[bond.atom_1-1].y_coor-atoms[bond.atom_2-1].y_coor)/bond.bond_distance,(atoms[bond.atom_1-1].z_coor-atoms[bond.atom_2-1].z_coor)/bond.bond_distance])
    return derib ,-derib


def derib_angle(atoms,angle):
    
    r_b_a=np.array([atoms[angle.atom_1-1].x_coor-atoms[angle.atom_2-1].x_coor,atoms[angle.atom_1-1].y_coor-atoms[angle.atom_2-1].y_coor,atoms[angle.atom_1-1].z_coor-atoms[angle.atom_2-1].z_coor])
    r_b_c=np.array([atoms[angle.atom_3-1].x_coor-atoms[angle.atom_2-1].x_coor,atoms[angle.atom_3-1].y_coor-atoms[angle.atom_2-1].y_coor,atoms[angle.atom_3-1].z_coor-atoms[angle.atom_2-1].z_coor])
    p=np.cross(r_b_a,r_b_c)
    deribA=((np.cross(-r_b_a,p))/(((np.linalg.norm(r_b_a))**2)*np.linalg.norm(p)))+((np.cross(r_b_c,p))/(((np.linalg.norm(r_b_c))**2)*np.linalg.norm(p)))
    deribB=((np.cross(r_b_a,p))/(((np.linalg.norm(r_b_a))**2)*np.linalg.norm(p)))
    deribC=((np.cross(-r_b_c,p))/(((np.linalg.norm(r_b_c))**2)*np.linalg.norm(p)))
    return deribA, deribB, deribC

def torsion_derib(atoms,torsion):
    point_1 = np.array([atoms[torsion.atom_1-1].x_coor, atoms[torsion.atom_1-1].y_coor, atoms[torsion.atom_1-1].z_coor])
    point_2 = np.array([atoms[torsion.atom_2-1].x_coor, atoms[torsion.atom_2-1].y_coor, atoms[torsion.atom_2-1].z_coor])
    point_3 = np.array([atoms[torsion.atom_3-1].x_coor, atoms[torsion.atom_3-1].y_coor, atoms[torsion.atom_3-1].z_coor])
    point_4 = np.array([atoms[torsion.atom_4-1].x_coor, atoms[torsion.atom_4-1].y_coor, atoms[torsion.atom_4-1].z_coor])

    
    rAB = point_2 - point_1
    #rAB= rAB/np.linalg.norm(rAB)
    rBC = point_3 - point_2
    #rBC= rBC/np.linalg.norm(rBC)
    rCD = point_4 - point_3
    #rBC= rCD/np.linalg.norm(rCD)
    rAC = point_3-point_1   #is correct not an error
    rBD=point_4-point_2
    
    t = np.cross(rAB, rBC)
    #t=t/np.linalg.norm(t)
    u = np.cross(rBC, rCD)
    #u=u/np.linalg.norm(u)


    deribA=np.cross(((np.cross(t,rBC))/((np.linalg.norm(t)**2)*np.linalg.norm(rBC))),rBC)
    deribB=np.cross(rAC,((np.cross(t,rBC))/((np.linalg.norm(t)**2)*np.linalg.norm(rBC))))+np.cross(((np.cross(-u,rBC))/((np.linalg.norm(u)**2)*np.linalg.norm(rBC))),rCD)
    deribC=np.cross(((np.cross(t,rBC))/((np.linalg.norm(t)**2)*np.linalg.norm(rBC))), rAB)+np.cross(rBD,((np.cross(-u,rBC))/((np.linalg.norm(u)**2)*np.linalg.norm(rBC))))
    deribD=np.cross(((np.cross(-u,rBC))/((np.linalg.norm(u)**2)*np.linalg.norm(rBC))),rBC)

    return deribA, deribB, deribC,deribD

def compute_wilson_B_mat (number_of_internal_coordinates, number_of_atoms,atoms,bonds,angles,torsions):
    wilson_B_mat=np.zeros((number_of_internal_coordinates, number_of_atoms * 3))
    i=0
    for bond in bonds:
        at_1=(bond.atom_1-1)*3
    
        at_2=3*(bond.atom_2-1)
        #print(at_1,at_2)
    
        at_1_array,at_2_array=derib_bond(atoms,bond)
        wilson_B_mat[i,at_1],wilson_B_mat[i,at_1+1],wilson_B_mat[i,at_1+2]=at_1_array[0],at_1_array[1],at_1_array[2]
        wilson_B_mat[i,at_2],wilson_B_mat[i,at_2+1],wilson_B_mat[i,at_2+2]=at_2_array[0],at_2_array[1],at_2_array[2]
        i=i+1

    for angle in angles:
        at_1=(angle.atom_2-1)*3
    
        at_2=3*(angle.atom_1-1)
        
        at_3=3*(angle.atom_3-1)
        #print(at_1,at_2,at_3)
        at_1_array,at_2_array,at_3_array=derib_angle(atoms,angle)
        wilson_B_mat[i,at_1],wilson_B_mat[i,at_1+1],wilson_B_mat[i,at_1+2]=at_1_array[0],at_1_array[1],at_1_array[2]
        wilson_B_mat[i,at_2],wilson_B_mat[i,at_2+1],wilson_B_mat[i,at_2+2]=at_2_array[0],at_2_array[1],at_2_array[2]
        wilson_B_mat[i,at_3],wilson_B_mat[i,at_3+1],wilson_B_mat[i,at_3+2]=at_3_array[0],at_3_array[1],at_3_array[2]
        i=i+1
    
    
    for torsion in torsions:
        at_1=(torsion.atom_2-1)*3
    
        at_2=3*(torsion.atom_1-1)
        
        at_3=3*(torsion.atom_3-1)
        at_4=3*(torsion.atom_4-1)
        #print(at_1,at_2,at_3,at_4)
        at_1_array,at_2_array,at_3_array,at_4_array=torsion_derib(atoms,torsion)
        wilson_B_mat[i,at_1],wilson_B_mat[i,at_1+1],wilson_B_mat[i,at_1+2]=at_2_array[0],at_2_array[1],at_2_array[2]
        wilson_B_mat[i,at_2],wilson_B_mat[i,at_2+1],wilson_B_mat[i,at_2+2]=at_1_array[0],at_1_array[1],at_1_array[2]
        wilson_B_mat[i,at_3],wilson_B_mat[i,at_3+1],wilson_B_mat[i,at_3+2]=at_3_array[0],at_3_array[1],at_3_array[2]
        wilson_B_mat[i,at_4],wilson_B_mat[i,at_4+1],wilson_B_mat[i,at_4+2]=at_4_array[0],at_4_array[1],at_4_array[2]
        i=i+1
    

    
    return wilson_B_mat

def compute_G_mat(wilson_B_mat):
    G=np.matmul(wilson_B_mat,np.transpose(wilson_B_mat))
    return G

def compute_gradient_cartesian(atoms, num_atoms, bonds,angles,torsions,VDW_interactions):
    gradient=np.zeros((num_atoms,3))
    gradient_stretch=np.zeros((num_atoms,3))
    gradient_bends=np.zeros((num_atoms,3))
    gradient_torsions=np.zeros((num_atoms,3))
    gradient_VDW=np.zeros((num_atoms,3))
    #print("PRINTING GRADIENT \n",gradient)
    for bond in bonds:
     #   print(bond.atom_1, bond.atom_2)
        if atoms[bond.atom_1-1].element== "H" or atoms[bond.atom_2-1].element== "H":
      #    print("1 bond distance diff",(bond.bond_distance-1.11))
          gradient_stretch[bond.atom_1-1,:]=(2*350*(bond.bond_distance-1.11)*derib_bond(atoms,bond)[0])+gradient_stretch[bond.atom_1-1,:]
          gradient_stretch[bond.atom_2-1,:]=-(2*350*(bond.bond_distance-1.11)*derib_bond(atoms,bond)[0])+gradient_stretch[bond.atom_2-1,:]
        else:
       #   print("2 bond distance diff",(bond.bond_distance-1.11))
          gradient_stretch[bond.atom_1-1,:]=(2*300*(bond.bond_distance-1.53)*derib_bond(atoms,bond)[0])+gradient_stretch[bond.atom_1-1,:]
          gradient_stretch[bond.atom_2-1,:]=-(2*300*(bond.bond_distance-1.53)*derib_bond(atoms,bond)[0])+gradient_stretch[bond.atom_2-1,:]
    #print("printing gradient stretchs")
    gradient=gradient+gradient_stretch
    
    for angle in angles:
            
            if atoms[angle.atom_1-1].element==atoms[angle.atom_2-1].element==atoms[angle.atom_3-1].element=="C":
                #print("1")
                theta=(109.5*(np.pi/180))
                gradient_bends[angle.atom_1-1]=(2*60*(angle.angle_value-theta)*derib_angle(atoms,angle)[1])+gradient_bends[angle.atom_1-1,:]
                gradient_bends[angle.atom_2-1]=(2*60*(angle.angle_value-theta)*derib_angle(atoms,angle)[0])+gradient_bends[angle.atom_2-1,:]
                gradient_bends[angle.atom_3-1]=(2*60*(angle.angle_value-theta)*derib_angle(atoms,angle)[2])+gradient_bends[angle.atom_3-1,:]
            if (atoms[angle.atom_1-1].element=="H" and atoms[angle.atom_2-1].element=="C" and atoms[angle.atom_3-1].element=="C") or (atoms[angle.atom_1-1].element=="C" and atoms[angle.atom_2-1].element=="C" and atoms[angle.atom_3-1].element=="H"):
                #print("2")
                theta=(109.5*(np.pi/180))
                gradient_bends[angle.atom_1-1]=(2*35*(angle.angle_value-theta)*derib_angle(atoms,angle)[1])+gradient_bends[angle.atom_1-1,:]
                gradient_bends[angle.atom_2-1]=(2*35*(angle.angle_value-theta)*derib_angle(atoms,angle)[0])+gradient_bends[angle.atom_2-1,:]
                gradient_bends[angle.atom_3-1]=(2*35*(angle.angle_value-theta)*derib_angle(atoms,angle)[2])+gradient_bends[angle.atom_3-1,:]
            if (atoms[angle.atom_1-1].element=="H" and atoms[angle.atom_2-1].element=="C" and atoms[angle.atom_3-1].element=="H"):
                #print("3")
                theta=(109.5*(np.pi/180))
                gradient_bends[angle.atom_1-1]=(2*35*(angle.angle_value-theta)*derib_angle(atoms,angle)[1])+gradient_bends[angle.atom_1-1,:]
                gradient_bends[angle.atom_2-1]=(2*35*(angle.angle_value-theta)*derib_angle(atoms,angle)[0])+gradient_bends[angle.atom_2-1,:]
                gradient_bends[angle.atom_3-1]=(2*35*(angle.angle_value-theta)*derib_angle(atoms,angle)[2])+gradient_bends[angle.atom_3-1,:]

    gradient=gradient+gradient_bends
    
    for torsion in torsions:
        #if atoms[torsion.atom_1-1].element==atoms[torsion.atom_4-1].element=="H":
            #print("1")
            gradient_torsions[torsion.atom_1-1]=-3*0.3*np.sin(3*torsion.torsion_angle)*torsion_derib(atoms,torsion)[0]+gradient_torsions[torsion.atom_1-1,:]
            gradient_torsions[torsion.atom_2-1]=-3*0.3*np.sin(3*torsion.torsion_angle)*torsion_derib(atoms,torsion)[1]+gradient_torsions[torsion.atom_2-1,:]
            gradient_torsions[torsion.atom_3-1]=-3*0.3*np.sin(3*torsion.torsion_angle)*torsion_derib(atoms,torsion)[2]+gradient_torsions[torsion.atom_3-1,:]
            gradient_torsions[torsion.atom_4-1]=-3*0.3*np.sin(3*torsion.torsion_angle)*torsion_derib(atoms,torsion)[3]+gradient_torsions[torsion.atom_4-1,:]
    
    gradient=gradient_torsions+gradient

    for VDW_interaction in VDW_interactions:
        distance_ij=np.sqrt((atoms[VDW_interaction.atom_1].x_coor-atoms[VDW_interaction.atom_2].x_coor)**2+(atoms[VDW_interaction.atom_1].y_coor-atoms[VDW_interaction.atom_2].y_coor)**2+(atoms[VDW_interaction.atom_1].z_coor-atoms[VDW_interaction.atom_2].z_coor)**2)

        if atoms[VDW_interaction.atom_1].element==atoms[VDW_interaction.atom_2].element=="H":
            A=4382.44
            B=22.932
        elif atoms[VDW_interaction.atom_1].element==atoms[VDW_interaction.atom_2].element=="C":
            A=946181.74
            B=514.714
        elif atoms[VDW_interaction.atom_1].element=="H" and atoms[VDW_interaction.atom_2].element=="C" or atoms[VDW_interaction.atom_1].element=="C" and atoms[VDW_interaction.atom_2].element=="H" :
            A=64393.99
            B=108.644

        gradient_VDW[VDW_interaction.atom_1]=np.array([(atoms[VDW_interaction.atom_1].x_coor-atoms[VDW_interaction.atom_2].x_coor)*(((-12*A)/(distance_ij**14))+((6*B)/(distance_ij**8))),(atoms[VDW_interaction.atom_1].y_coor-atoms[VDW_interaction.atom_2].y_coor)*(((-12*A)/(distance_ij**14))+((6*B)/(distance_ij**8))) ,(atoms[VDW_interaction.atom_1].z_coor-atoms[VDW_interaction.atom_2].z_coor)*(((-12*A)/(distance_ij**14))+((6*B)/(distance_ij**8))) ])+gradient_VDW[VDW_interaction.atom_1,:]
        #gradient_VDW[VDW_interaction.atom_2]=np.array([(atoms[VDW_interact.atom_1].x_coor-atoms[VDW_interact.atom_2].x_coor)*(((-12*A)/(distance_ij**14))+((6*B)/(distance_ij**8))),(atoms[VDW_interact.atom_1].y_coor-atoms[VDW_interact.atom_2].y_coor)*(((-12*A)/(distance_ij**14))+((6*B)/(distance_ij**8))) ,(atoms[VDW_interact.atom_1].z_coor-atoms[VDW_interact.atom_2].z_coor)*(((-12*A)/(distance_ij**14))+((6*B)/(distance_ij**8))) ])+gradient_VDW[VDW_interaction.atom_2]
        gradient_VDW[VDW_interaction.atom_2]=np.array([(atoms[VDW_interaction.atom_2].x_coor-atoms[VDW_interaction.atom_1].x_coor)*(((-12*A)/(distance_ij**14))+((6*B)/(distance_ij**8))),(atoms[VDW_interaction.atom_2].y_coor-atoms[VDW_interaction.atom_1].y_coor)*(((-12*A)/(distance_ij**14))+((6*B)/(distance_ij**8))) ,(atoms[VDW_interaction.atom_2].z_coor-atoms[VDW_interaction.atom_1].z_coor)*(((-12*A)/(distance_ij**14))+((6*B)/(distance_ij**8))) ])+gradient_VDW[VDW_interaction.atom_2,:]
    
    gradient=gradient+gradient_VDW

    return gradient

def compute_gradient_internal(cartesian_gradient,B_matrix,G_inverse):
    cartesian_gradient_array=np.zeros((np.shape(cartesian_gradient)[0]*np.shape(cartesian_gradient)[1],1))
    c=0
    for i in range(np.shape(cartesian_gradient)[0]):
        
        for j in range(np.shape(cartesian_gradient)[1]):
            
            
            cartesian_gradient_array[c]=cartesian_gradient[i,j]
            
            c=c+1
    
    internal_gradient=G_inverse @ B_matrix @ cartesian_gradient_array 
    return internal_gradient

###############################################################################################################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
###############################################################################################################################################################

def internal_to_cartesians(atoms,bonds,angles,torsions,connectivity_data,wilson_B_mat, G_inverse,s_q_k, file_name):
    #print("  Iterative determination of optimal Cartesian coordinates:  ")
    #print("s_q_k \n", np.shape(np.transpose(wilson_B_mat)),np.shape(G_inverse),np.shape(s_q_k))
    dx_max=0.000010   #set the maximun DX we will allow
    with open(file_name+".out", "a") as file:
        file.write(f"\n Iterative determination of optimal Cartesian coordinates: \n")
        file.write(f"   tolerated diff for Cartsian coordinates:    {dx_max} \n")

    ###############Lets set the introduced cartesian coordinates to an array to work

    cartesians_0=np.zeros((len(atoms)*3,1))
    i=0
    for atom in atoms:
        cartesians_0[i,0],cartesians_0[i+1,0],cartesians_0[i+2,0]=atom.x_coor,atom.y_coor,atom.z_coor
        i=i+3

    internal_0=np.array(compute_internal_coordinates(atoms,bonds,connectivity_data))
    ##############lets make the jump in cartesian coordinates ##########

    cartesians_1=cartesians_0+np.transpose(wilson_B_mat) @ G_inverse @ s_q_k

    #print(np.shape(cartesians_1))

    ############## lets compute the desired internals  ###############

    desired_internals=np.array(compute_internal_coordinates(atoms,bonds,connectivity_data))
    
    for i in range(len(internal_0)):
        desired_internals[i]=internal_0[i]+s_q_k[i,0]
    
    #print(desired_internals)
    #time.sleep(5)
    for i in range(len(bonds)+len(angles),len(bonds)+len(angles)+len(torsions)):
        print(desired_internals[i])
        if abs(desired_internals[i]) > np.pi and desired_internals[i] > 0:
            desired_internals[i]=desired_internals[i]-2*np.pi
            desired_internals[i]
        elif abs(desired_internals[i]) > np.pi and desired_internals[i] < 0:
            desired_internals[i]=desired_internals[i]+2*np.pi
            desired_internals[i]

    #print("desired_internals\n", desired_internals)
    ############# From here we can compute the, initial dx = BTG-sk #######

    DX=cartesians_1-cartesians_0


    with open(file_name+".out", "a") as file:
        file.write(f"\n Initially predicted dx = BTG-sk: \n")
        file.write("\n")
        for row in DX:
            file.write(' '.join(f"{x:.6f}"' ' for x in row))
        file.write("\n")

    #while np.max(np.abs(DX))>=dx_max: #Che

        ############# Now For this cartesians we can generate a new q_(k+1)^(j+1) For this we will actualize the atoms
    i=0
    for atom in atoms:
        atom.x_coor,atom.y_coor,atom.z_coor=cartesians_1[i,0],cartesians_1[i+1,0],cartesians_1[i+2,0]
        i=i+3
        
    internal_1=np.array(compute_internal_coordinates(atoms,bonds,connectivity_data))

    with open(file_name+".out", "a") as file:
            file.write(f"\n current set of internals q_(k+1)^(j+1): \n")
            for element in internal_1:
                file.write(f"{element:.5f}  ")
            file.write("\n")

    k=0
    #while k<=1:
    
    while np.max(np.abs(DX))>=dx_max: 

        
        ############ Now we can generate a new s_(k+1)^(j+1)
        
        for i in range(np.size(s_q_k)):
                
                s_q_k[i,0]=desired_internals[i]-internal_1[i]
                if i>=len(bonds)+len(angles):
                    #time.sleep(1)
                    print(i)
                    print("desired_internals[i]  :",desired_internals[i])
                    print("internals[i]  :",internal_1[i])
                    print("sqk[i]  :",s_q_k[i,0])
                    if abs(s_q_k[i,0])>np.pi and s_q_k[i,0]>0:
                        s_q_k[i,0]=s_q_k[i,0]-2*np.pi
                        print("sqk[i]  :",s_q_k[i,0])
                    elif abs(s_q_k[i,0])>np.pi and s_q_k[i,0]<0:
                        s_q_k[i,0]=s_q_k[i,0]+2*np.pi
                        print("sqk[i]  :",s_q_k[i,0])


                
        
        with open(file_name+".out", "a") as file:
            file.write(f"\n difference between these internals and the desired internals, s_q,k^j+1: \n")
            for row in s_q_k:
                file.write(' '.join(f"{x:.6f}"' ' for x in row))
            file.write("\n")

        #print("s_q_k \n", s_q_k )
        ########### We calculate the new cartesians using the new s_(k+1)^(j+1)
        #print("helle \n",np.shape(np.transpose(wilson_B_mat)), np.shape(G_inverse) ,np.shape(s_q_k))
        cartesians_2=cartesians_1+np.transpose(wilson_B_mat) @ G_inverse @ s_q_k
        #print(np.shape(cartesians_2),cartesians_2)
        
        ########## We evaluate the change ######################3
        DX=cartesians_2-cartesians_1
        
        
        ############ we actualize the atoms ##########################
        i=0
        for atom in atoms:
            atom.x_coor,atom.y_coor,atom.z_coor=cartesians_2[i,0],cartesians_2[i+1,0],cartesians_2[i+2,0]
            i=i+3
        
        ############ with new coordinates we compute the new internals

        internal_1=compute_internal_coordinates(atoms,bonds,connectivity_data)

        with open(file_name+".out", "a") as file:
            file.write(f"\n current set of internals q_(k+1)^(j+1): \n")
            for element in internal_1:
                file.write(f"{element:.5f}  ")
            file.write("\n")
        #print(np.shape(internal_1),internal_1)
        ############ FOR THE NEW CYCLE WE SEET THE NEW CARTESIANS AS CARTESIANS_1
        cartesians_1=cartesians_2

        #print(np.max(np.abs(DX)))
        #time.sleep(5)
        

    cartesian_matrix=np.zeros((len(atoms),3))
    j=0
    for i in range(0,len(cartesians_1),3):
        cartesian_matrix[j,:]=cartesians_1[i,0],cartesians_1[i+1,0],cartesians_1[i+2,0]
        j=j+1
      
    
    return cartesian_matrix




###############################################################################################################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
###############################################################################################################################################################


def compute_step_internal(internal_hessian_inv, internal_gradient,file_name):
    
    p_q_k=-internal_hessian_inv @ internal_gradient

    with open(file_name+".out", "a") as file:
        file.write(f"\n Predicted update step in internal coordinates s_k (prior to possible scaling): \n")
        file.write("\n")
        for row in p_q_k:
            file.write(' '.join(f"{x:.6f}"' ' for x in row))
            #file.write("\n")
    
    MAX_RMS=0.020
    RMS=0                               
    for i in range(np.size(p_q_k)):
        RMS=RMS+p_q_k[i,0]*p_q_k[i,0]
    RMS=np.sqrt(RMS/np.size(p_q_k))         
    


    if RMS>MAX_RMS:
        #print("Predicted step is too long: RMS length: \n",RMS)
        p_q_k=p_q_k*(MAX_RMS/RMS)
    #print("Scaled update step in internal coordinates s_k: \n", p_q_k)
    
    return p_q_k



################# END DEFINITION OF FUNCTIONS ##########################

################# MAIN CODE ##########################################


def main():
    parser = argparse.ArgumentParser()
    molecule_file = parser.add_argument("file_path", type=str, help="introduce file path")  #get molecule file

    args = parser.parse_args()
    file_path = args.file_path

    file_name = os.path.splitext(os.path.basename(file_path))[0]

    atom_data,conectivity_data, num_atoms, num_bonds=extract_file_data(file_path)
    atoms=get_atoms(atom_data)  #extract atoms data



    bonds=get_bonds(num_bonds,conectivity_data)   #GET THE BONDS

    for bond in bonds:
        compute_bond_energies(atoms,bond)         #COMPUTE THE ENERGIES

    angles=get_angles(atoms,bonds)                #BET THE ANGLES

    for angle in angles:
        compute_angle_energies(atoms,angle)       #COMPUTE ANGLE ENERGIES
    
    torsions=get_torsions(atoms,bonds,angles)     
    
    for torsion in torsions:
        compute_torsion_energies(torsion)
    
    VDW_interactions=get_VDW_interactions(atoms,bonds,angles)

    for VDW_interaction in VDW_interactions:
        compute_VDW_energy(VDW_interaction,atoms)
    
    num_internal_coor=len(bonds)+len(angles)+len(torsions)
    internal_coor=np.zeros((num_internal_coor))             #generate array to store internal coordinates
    t=0
    for bond in bonds:
        internal_coor[t]=bond.bond_distance                 #store bonds
        t=t+1

    for angle in angles:
        internal_coor[t]=angle.angle_value                  #store angles
        t=t+1
    
    for torsion in torsions:
        internal_coor[t]=torsion.torsion_angle              #store torsions
        t=t+1
    


    with open(file_name+".out", "w") as file:
        file.write(f"ATOM DATA ({num_atoms}) \n")
        for atom in atoms:
            file.write(f"Atom {atom.atm_num}: {atom.element} {atom.x_coor}, {atom.y_coor}, {atom.z_coor} \n") 
        file.write("\n")
        file.write(f"BONDS DATA ({num_bonds}) \n")
        for bond in bonds:
            file.write(f" {atoms[bond.atom_1-1].element} {bond.atom_1} --- {atoms[bond.atom_2-1].element} {bond.atom_2} bond_type:  {bond.bond_type} bond_distance:  {bond.bond_distance} Å bond energy:  {bond.bond_energy} Kcal/mol. \n")
        file.write("\n")
        file.write(f"Angles DATA ({len(angles)}) \n")
        for angle in angles:
            file.write(f"{atoms[angle.atom_1-1].element} {angle.atom_1} --- {atoms[angle.atom_2-1].element} {angle.atom_2} --- {atoms[angle.atom_3-1].element} {angle.atom_3} angle value in radias = {angle.angle_value} angle value in degrees = {angle.angle_value *(180/np.pi)} angle energy: {angle.angle_energy}. \n")
        file.write("\n")
        file.write(f"TORSIONS DATA ({len(torsions)}) \n")
        for torsion in torsions:
            file.write(f"torsional angles. atom 1: {atoms[torsion.atom_1-1].element} {torsion.atom_1} --- {atoms[torsion.atom_2-1].element}  {torsion.atom_2} --- {atoms[torsion.atom_3-1].element} {torsion.atom_3} --- {atoms[torsion.atom_4-1].element} {torsion.atom_4} torsion angle in radians {torsion.torsion_angle} torsion angle in degrees {torsion.torsion_angle*(180/np.pi)}  torsion energy {torsion.torsion_energy}. \n")
        file.write("\n")
        file.write(f"VDW INTERACTIONS DATA ({len(VDW_interactions)}) \n")
        for VDW_interact in VDW_interactions:
            file.write(f"{atoms[VDW_interact.atom_1].element} {VDW_interact.atom_1+1} <---> {atoms[VDW_interact.atom_2].element} {VDW_interact.atom_2+1} distance= {VDW_interact.distance} Energy= {VDW_interact.VDW_energy} Kcal·mol⁻¹. \n" )
        file.close()



    



    

    ####################    Compute THE WILSON B MATRIX      ######################

    wils_B_mat=compute_wilson_B_mat(num_internal_coor, num_atoms,atoms,bonds,angles,torsions)
    with open(file_name+".out", "a") as file:
        file.write(f"\n  INITIAL WILSON B MATRIX \n")
        file.write("\n")
        for row in wils_B_mat:
            file.write(' '.join(f"{x:.5f}" for x in row))
            file.write("\n")


    ###################  FROM THE WILSON B MATRIX COMPUTE THE G MATRIX

    G_mat=compute_G_mat(wils_B_mat)

    with open(file_name+".out", "a") as file:
        file.write(f"\n INITIAL  G MATRIX \n")
        file.write("\n")
        for row in G_mat:
            file.write(' '.join(f"{x:.5f}" for x in row))
            file.write("\n")


    ################## COMPUTE THE (PSEUDO)INVERSE OF THE G MATRIX 

    G_mat_inv=np.linalg.pinv(G_mat)
    with open(file_name+".out", "a") as file:
        file.write(f"\n INITIAL G MATRIX INVERSE \n")
        file.write("\n")
        for row in G_mat_inv:
            file.write(' '.join(f"{x:.5f}" for x in row))
            file.write("\n")

    ################  COMPUTE THE GRADIENT IN CARTESIANS ##################

    cartesian_gradient=compute_gradient_cartesian(atoms, num_atoms, bonds,angles,torsions,VDW_interactions)
    #print("cartesian gradien\n",cartesian_gradient)
    

    ###############   COMPUTE THE GRADIENT IN INTERNALS  ##################

    internal_gradient=compute_gradient_internal(cartesian_gradient,wils_B_mat,G_mat_inv)  #computing the internal gradient from the cartesian

    with open(file_name+".out", "a") as file:
        file.write(f"\n INITIAL GRADIENT INTERNAL COORDINATES \n")
        file.write("\n")
        for element in internal_gradient:
            file.write(f"{element}, ")
        file.write("\n")

    RMS=0                               
    for i in range(np.size(internal_gradient)):
        RMS=RMS+internal_gradient[i]*internal_gradient[i]
        RMS=np.sqrt(RMS/np.size(internal_gradient)) 
    
    with open(file_name+".out", "a") as file:
        file.write(f"\n GRMS = {RMS}\n")
        file.write("\n")
    ##############  INITIALIZE HESSIAN INVERSE GUESS   ####################

    internal_hessian_inv=np.zeros((num_internal_coor,num_internal_coor))

    for i in range(len(bonds)):
        internal_hessian_inv[i,i]=0.00167

    for i in range(len(bonds),len(bonds)+len(angles)):
        internal_hessian_inv[i,i]=0.00667
    for i in range(len(bonds)+len(angles),len(bonds)+len(angles)+len(torsions)):
        internal_hessian_inv[i,i]=0.01250

    with open(file_name+".out", "a") as file:
        file.write(f"\n  INITIALIZE M HESSIAN MATRIX INVERSE GUESS \n")
        file.write("\n")
        for row in internal_hessian_inv:
            file.write(' '.join(f"{x:.5f}" for x in row))
            file.write("\n")

    j=0
    E1=Compute_tot_energy(atoms,bonds,angles,torsions,VDW_interactions,file_name)
    E0=0
    j=0
    #while 10^(-3)<=RMS:
    while abs(E1-E0)>10**(-8):

        E0=E1

        with open(file_name+".out", "a") as file:
            file.write(f"\n################################## \n")
            file.write("# Start of Geometry Optimization #\n")
            file.write("################################## \n")

        p_q_k=compute_step_internal(internal_hessian_inv,internal_gradient,file_name)
        alpha=1
        s_q_k=alpha*p_q_k
        #print("sqk \n",s_q_k)

        with open(file_name+".out", "a") as file:
            file.write(f"\n Scaled update step in internal coordinates s_k: \n")
            file.write("\n")
            for i in range(np.size(s_q_k)):
                file.write(f" {s_q_k[i,0]:.6f} ")
            file.write("\n")

        dx_max_cart=0.00001        

        with open(file_name+".out", "a") as file:
            file.write(f"\n  Iterative determination of optimal Cartesian coordinates: \n")
            file.write(f"  tolerated diff for Cartsian coordinates:  {dx_max_cart} \n")
       
       
        cartesian_mat=internal_to_cartesians(atoms,bonds,angles,torsions,conectivity_data,wils_B_mat,G_mat_inv,s_q_k,file_name)

        for i in range(len(atoms)):
            atoms[i].x_coor,atoms[i].y_coor,atoms[i].z_coor=cartesian_mat[i,:]

        for bond in bonds:
            compute_bond_energies(atoms,bond)

        angles=get_angles(atoms,bonds)                #GET THE ANGLES
        for angle in angles:
            compute_angle_energies(atoms,angle)       #COMPUTE ANGLE ENERGIES

        torsions=get_torsions(atoms,bonds,angles)     
        for torsion in torsions:
            compute_torsion_energies(torsion)

        VDW_interactions=get_VDW_interactions(atoms,bonds,angles)
        for VDW_interaction in VDW_interactions:
            compute_VDW_energy(VDW_interaction,atoms)

        num_internal_coor=len(bonds)+len(angles)+len(torsions)
        new_internal_coor=np.zeros((num_internal_coor))

        t=0
        for bond in bonds:
            new_internal_coor[t]=bond.bond_distance
            t=t+1
        for angle in angles:
            new_internal_coor[t]=angle.angle_value
            t=t+1
        for torsion in torsions:
            new_internal_coor[t]=torsion.torsion_angle
            t=t+1

        with open(file_name+".out", "a") as file:
            file.write(f"\n NEW ATOMIC COORDINATES ({num_atoms}) \n")
            file.write("\n")
            for atom in atoms:
                file.write(f"Atom {atom.atm_num}: {atom.element} {atom.x_coor}, {atom.y_coor}, {atom.z_coor} \n") 
            file.write("\n")

        with open(file_name+".out", "a") as file:
            file.write(f"\n New set of internals q (note these are the ones that correspond to the best fit Cartesians): \n")
            file.write("\n")
            for element in new_internal_coor:
                file.write(f"{element:.5f}, ") 
            file.write("\n")

        s_k=new_internal_coor-internal_coor
        #print("heaven",s_k,"hell")

        E1=Compute_tot_energy(atoms,bonds,angles,torsions,VDW_interactions,file_name)


        wils_B_mat=compute_wilson_B_mat(num_internal_coor,num_atoms,atoms,bonds,angles,torsions)
        with open(file_name+".out", "a") as file:
            file.write(f"\n {j} ITERATION WILSON B MATRIX \n")
            file.write("\n")
            for row in wils_B_mat:
                file.write(' '.join(f"{x:.5f}" for x in row))
                file.write("\n")

        G_mat=compute_G_mat(wils_B_mat)
        G_mat_inv=np.linalg.pinv(G_mat)
        with open(file_name+".out", "a") as file:
            file.write(f"\n {j} ITERATION  G MATRIX \n")
            file.write("\n")
            for row in G_mat:
                file.write(' '.join(f"{x:.5f}" for x in row))
                file.write("\n")

        with open(file_name+".out", "a") as file:
            file.write(f"\n {j} ITERATION G MATRIX INVERSE \n")
            file.write("\n")
            for row in G_mat_inv:
                file.write(' '.join(f"{x:.5f}" for x in row))
                file.write("\n")
        new_cartesian_gradient=compute_gradient_cartesian(atoms,num_atoms,bonds,angles,torsions,VDW_interactions)
        new_internal_gradient=compute_gradient_internal(new_cartesian_gradient,wils_B_mat,G_mat_inv)

        with open(file_name+".out", "a") as file:
            file.write(f"\n {j} ITERATION GRADIENT INTERNAL COORDINATES \n")
            file.write("\n")
            for element in internal_gradient:
                file.write(f"{element}, ")
            file.write("\n")


        RMS=0                               
        for i in range(np.size(new_internal_gradient)):
            RMS=RMS+new_internal_gradient[i]*new_internal_gradient[i]
            RMS=np.sqrt(RMS/np.size(new_internal_gradient)) 
        
        with open(file_name+".out", "a") as file:
            file.write(f"\n GRMS = {RMS}\n")
            file.write("\n")
            

        y_q_k=new_internal_gradient-internal_gradient


        with open(file_name+".out", "a") as file:
            file.write(f"\n Change in the gradient in terms of the internal coordinates: (vector yk) \n")
            file.write("\n")
            for element in y_q_k:
                file.write(f"{element}, ")
            file.write("\n")

        with open(file_name+".out", "a") as file:
            file.write(f"\n Change in the structure in terms of the internal coordinates: (vector sk) \n")
            file.write("\n")
            for element in s_k:
                file.write(f"{element}, ")
            file.write("\n")
        v_q_k=internal_hessian_inv @ y_q_k
        with open(file_name+".out", "a") as file:
            file.write(f"\n Product of the inverse Hessian and the change in gradient, vector vk: \n")
            file.write("\n")
            for element in v_q_k:
                file.write(f"{element}, ")
            file.write("\n")

        for i in range(len(s_k)):
            s_q_k[i,0]=s_k[i]     #just to give proper shape

        #print(s_q_k)

        #print("shape outer",np.shape(np.outer(s_q_k,s_q_k)))
        #print("inverse",np.linalg.pinv(s_q_k @ np.transpose(y_q_k )))
        #print("operation",((s_q_k @ np.transpose(y_q_k )+ y_q_k @ np.transpose(v_q_k) )@ np.outer(s_q_k,s_q_k))@(np.linalg.pinv(s_q_k @ np.transpose(y_q_k ))@np.linalg.pinv(s_q_k @ np.transpose(y_q_k ))))

        
        dot_s_y = np.dot(s_q_k.T, y_q_k)  # s_q_k · y_q_k
        dot_s_v = np.dot(s_q_k.T, v_q_k)  # s_q_k · v_q_k
        dot_y_v = np.dot(y_q_k.T, v_q_k)  # y_q_k · v_q_k
       
        outer_s_s = np.outer(s_q_k, s_q_k)  # s_q_k ⊗ s_q_k
        outer_v_s = np.outer(v_q_k, s_q_k)  # v_q_k ⊗ s_q_k
        outer_s_v = np.outer(s_q_k, v_q_k)  # s_q_k ⊗ v_q_k

        numerator = (dot_s_y + dot_y_v) 
        denominator = dot_s_y**2
        first_part = (numerator / denominator) * outer_s_s

        second_part = (outer_v_s + outer_s_v) / dot_s_y

        internal_hessian_inv = internal_hessian_inv + first_part - second_part

        with open(file_name+".out", "a") as file:
            file.write(f"\n {j} ITERATION M HESSIAN MATRIX INVERSE \n")
            file.write("\n")
            for row in internal_hessian_inv:
                file.write(' '.join(f"{x:.6f}" for x in row))
                file.write("\n")
        j=j+1 ########## JUST HERE TRYING TO WATCH GOING ON


        internal_coor=new_internal_coor

        Compute_tot_energy(atoms,bonds,angles,torsions,VDW_interactions,file_name)
        
        internal_gradient=new_internal_gradient
        print(f"iteration {j}")

    return

main()

