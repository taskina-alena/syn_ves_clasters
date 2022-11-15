import argparse
from Topology.linker_class import make_linked_vesicles
import itertools
import os
import numpy as np

def write_in_script(eps_attr, eps_rep, inter_range, folderPattern, filePattern, attr_with_rep, sigma):

    filename = f"Input/Scripts/Input_{filePattern}.in"

    f = open(filename, "w")

    lj_factor = 2**(1/6)

    global_cutoff = 3

    eps = {'attr' : eps_attr, 'rep' : eps_rep}
    types = {'vesicle':1, 'linker':2, 'synapsin':3}


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

    for i in types.keys():
        for j in types.keys():
            if types[i] <= types[j]:
                if i=='linker' and j=='synapsin':
                    if attr_with_rep:
                        f.write(f"pair_coeff {types[i]} {types[j]} {eps['attr']} { lj_factor * (sigma[i] + sigma[j]) } { inter_range + lj_factor * (sigma[i] + sigma[j]) } wca \n")  ## eps sigma cutoff
                    else:
                        f.write(f"pair_coeff {types[i]} {types[j]} {eps['attr']} {lj_factor * (sigma[i] + sigma[j])}\n")  ## eps sigma

                else:
                    f.write(f"pair_coeff {types[i]} {types[j]} {eps['rep']} {lj_factor * (sigma[i] + sigma[j]) } { lj_factor * (sigma[i] + sigma[j]) } wca \n")  ## eps sigma cutoff
    # if cutoff=sigma, only repulsive interactions more or less (between all atoms)

    f.write(
    '''
    ##======================================##
    ## Setting up groups
    ##======================================##
    \n''')

    f.write("group \t linked_vesicle \t type \t 1:2 \n")
    f.write("group \t synapsin \t type \t 3 \n")

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

    rigid_molecules = True

    sim_seed = 5
    dump_step = 1000
    eq_steps = 1e3
    run_steps = 1e6

    eps_attr = 1.1
    eps_rep = 0.8
    inter_range = 2
    attr_with_rep = False
    seed = 100

    num_synapsin = 400 #linker concentration
    num_vesicles = 40 # protein concentration
    num_linkers = 4 * num_vesicles # linkers concentration
    side_box = 25
    sigma = {'vesicle': 2, 'linker': 0.5, 'synapsin': 0.5}

    system = make_linked_vesicles(num_synapsin, num_linkers, num_vesicles, side_box, sigma)
    system.make_linked_vesicles()
    system.make_synapsin()

    filePattern = f"2d_eps_attr{eps_attr}_eps_rep{eps_rep}_range{inter_range}_attr_with_rep{attr_with_rep}_n_synapsin{num_synapsin}_n_vesicles{num_vesicles}"
    folderPattern = f"Results_eps_attr{eps_attr}_eps_rep{eps_rep}_range{inter_range}_attr_with_rep{attr_with_rep}_n_synapsin{num_synapsin}_n_vesicles{num_vesicles}"

    if not os.path.exists(folderPattern):
        os.makedirs(folderPattern)

    configName = "Input/Configuration/Config_%s.dat" % filePattern

    header = ["LAMMPS Description \n \n",
              "\t " + str(system.numAll) + " atoms \n \t " + str(system.numBonds) +
              " bonds \n \t " + str(system.numAngles) + " angles \n \t 0 dihedrals \n \t 0 impropers \n",
              "\n \t "+str(system.numTypes)+" atom types \n \t 2 bond types \n \t 0 angle types \n \t 0 dihedral types \n \t 0 improper types \n",
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

    if not rigid_molecules:
        for item in system.bonds:
            f.write(item)

    f.close()
    write_in_script(eps_attr, eps_rep, inter_range, folderPattern, filePattern, attr_with_rep, sigma)
















