#import argparse
from Topology.linker_class import make_linked_vesicles
import itertools
import os
import numpy as np

def write_in_script(eps_attr, eps_rep, eps_rep_weak, inter_range, folderPattern, filePattern, attr_with_rep, sigma):

    #TODO introduce different size of vesicles

    filename = f"Input/Scripts/Input_{filePattern}.in"

    f = open(filename, "w")

    lj_factor = 2**(1/6)

    global_cutoff = 3

    eps = {'attr' : eps_attr, 'rep' : eps_rep, 'rep_weak' : eps_rep_weak}


    f.write("units lj\n")
    f.write("dimension 2\n")
    f.write("atom_style full\n")
    f.write("boundary p p p\n")
    f.write("read_data %s \n" %configName)
    f.write("neighbor               0.3 bin\n")
    f.write("neigh_modify           every 1 delay 1\n")

    f.write(
    '''
    ##======================================##
    ## Interactions
    ##======================================##
    \n''')

    f.write(f"pair_style cosine/squared {global_cutoff} \n")

    for i in system.types.keys():
        for j in system.types.keys():
            if (system.types[i] <= system.types[j]) and system.types[i] and system.types[j]:
                if (i=='linker' and j=='synapsin') or (j=='linker' and i =='synapsin'):
                    if attr_with_rep:
                        f.write(f"pair_coeff {system.types[i]} {system.types[j]} {eps['attr']} { lj_factor * (sigma[i] + sigma[j]) } { inter_range + lj_factor * (sigma[i] + sigma[j]) } wca \n")  ## eps sigma cutoff
                    else:
                        f.write(f"pair_coeff {system.types[i]} {system.types[j]} {eps['attr']} {lj_factor * (sigma[i] + sigma[j])} { inter_range + lj_factor * (sigma[i] + sigma[j]) } \n")  ## eps sigma cutoff
                elif (i == 'vesicle' and j == 'synapsin') or (j == 'vesicle' and i == 'synapsin'):
                    f.write(f"pair_coeff {system.types[i]} {system.types[j]} {eps['rep_weak']} {lj_factor * (sigma[i] + sigma[j])} {lj_factor * (sigma[i] + sigma[j])} wca \n")  ## eps sigma cutoff
                else:
                    f.write(f"pair_coeff {system.types[i]} {system.types[j]} {eps['rep']} {lj_factor * (sigma[i] + sigma[j]) } { lj_factor * (sigma[i] + sigma[j]) } wca \n")  ## eps sigma cutoff
    # if cutoff=sigma, only repulsive interactions more or less (between all atoms)

    f.write(
        '''
        ##======================================##
        ## Bonds
        ##======================================##
        \n''')

    f.write("bond_style harmonic \n")
    #TODO make energy big
    f.write(f"bond_coeff  1 1 {sigma['vesicle'] + sigma['linker']}\n")

    f.write(
        '''
        ##======================================##
        ## Angles
        ##======================================##
        \n''')

    f.write("angle_style harmonic \n")

    #TODO check different energy
    f.write(f"angle_coeff  1 1 90\n")

    f.write(
    '''
    ##======================================##
    ## Setting up groups
    ##======================================##
    \n''')

    if rigid_molecules:
        f.write("group \t linked_vesicle \t type \t" + str(system.types['vesicle']) + ":" + str(system.types['linker']) + "\n")
        f.write("group \t synapsin \t type \t" + str(system.types['synapsin']) + " \n")

    f.write(
    '''
    ##======================================##
    ## Computes
    ##======================================##
    \n''')

    f.write(
    '''
    ##======================================##
    ## Fixes
    ##======================================##
    \n''')

    f.write("velocity all create 1.0 1\n") # temp seed

    f.write(f"log {folderPattern}/Log_{filePattern}.dat\n")

    f.write(f"fix fLANG all langevin 1.0 1.0 1.0 {seed}\n") # temp_start temp_end damp

    if rigid_molecules:
        f.write("neigh_modify exclude molecule/intra all\n")
        f.write("fix rigidNVE linked_vesicle rigid/nve molecule\n")
        f.write("fix NVE synapsin nve 1.0 1.0 1.0\n")
    else:
        f.write("fix NVE all nve 1.0 1.0 1.0\n")

    f.write(
    '''
    ##======================================##
    ## Output
    ##======================================##
    \n''')

    f.write("thermo %i\n"%dump_step)
    f.write("thermo_style custom step temp press etotal epair\n")
    f.write("thermo_modify flush yes\n")
    f.write("timestep 0.01\n")

    f.write(
    '''
    ##======================================##
    ## Equilibration before output
    ##======================================##
    \n''')

    f.write(
    '''
    ##======================================##
    ## Run with output
    ##======================================##
    \n''')

    f.write("fix enforce_2d all enforce2d\n")
    f.write(f"dump 2 all custom {dump_step} {folderPattern}/Movie_{filePattern}.xyz id type mol x y z \n")
    f.write("run %i\n"%run_steps)
    f.close()

if __name__ == "__main__":

    rigid_molecules = False

    sim_seed = 5
    dump_step = 1000
    eq_steps = 1e3
    run_steps = 1e6

    eps_attr = 6
    eps_rep = 1
    eps_rep_weak = 1
    inter_range = 0.25
    attr_with_rep = False
    seed = 100
    np.random.seed(seed)

    num_synapsin = 4000
    num_vesicles = 100
    num_linkers = num_vesicles*4
    side_box = 200
    sigma = {'vesicle': 4, 'linker': 0.45, 'synapsin': 0.45}

    system = make_linked_vesicles(num_synapsin, num_linkers, num_vesicles, side_box, sigma)
    system.make_single_particle('synapsin')
    #system.make_single_particle('linker')
    system.make_linked_vesicles()


    filePattern = f"2d_eps_attr{eps_attr}_eps_rep{eps_rep}_eps_rep_weak{eps_rep_weak}_range{inter_range}_attr_with_rep{attr_with_rep}_n_synapsin{num_synapsin}_n_vesicles{num_vesicles}"
    folderPattern = f"Results_eps_attr{eps_attr}_eps_rep{eps_rep}_eps_rep_weak{eps_rep_weak}_range{inter_range}_attr_with_rep{attr_with_rep}_n_synapsin{num_synapsin}_n_vesicles{num_vesicles}"

    if not os.path.exists(folderPattern):
        os.makedirs(folderPattern)

    configName = "Input/Configuration/Config_%s.dat" % filePattern

    header = ["LAMMPS Description \n \n",
              "\t " + str(system.numAll) + " atoms \n \t " + str(system.numBonds) +" bonds \n \t " +
              str(system.numAngles) + " angles \n \t 0 dihedrals \n \t 0 impropers \n",
              "\n \t "+str(system.numTypes)+" atom types \n \t 1 bond types \n \t 1 angle types \n \t 0 dihedral types \n \t 0 improper types \n",
              "\n \t " + str(-system.Lx*0.5) + " " + str(system.Lx*0.5) + " xlo xhi\n \t", str(-system.Ly*0.5) + " " + str(system.Ly*0.5) + " ylo yhi \n \t",
              str(-system.Lz*0.5) + " " + str(system.Lz*0.5) + " zlo zhi\n"]


    header.append("\nMasses \n \n")

    for i in range(len(system.type_mass_list)):
        header.append("\t %i %.4f \n" % (system.type_mass_list[i][0], system.type_mass_list[i][1]))

    f = open(configName, "w")

    for item in header:
        f.write(item)

    for item in system.coords:
        f.write(item)

    for item in system.bonds:
            f.write(item)

    for item in system.angles:
            f.write(item)

    f.close()
    write_in_script(eps_attr, eps_rep, eps_rep_weak, inter_range, folderPattern, filePattern, attr_with_rep, sigma)
















