def write_in_script(system, folderPattern, filePattern, configName, dump_step, run_steps, create_prob, break_prob, intra_range, timestep):

    seed = 100
    filename = f"Input_span/Scripts/Input_{filePattern}.in"

    f = open(filename, "w")

    f.write("units lj\n")
    f.write("dimension 2\n")
    f.write("atom_style molecular\n")
    f.write("boundary p p p\n")
    f.write(f"read_data {configName} \n")
    f.write("neighbor               0.3 bin\n")
    f.write("neigh_modify           every 1 delay 1\n")

    f.write(
    '''
    ##======================================##
    ## Interactions
    ##======================================##
    \n''')

    for item in system.interactions():
        f.write(item)

    f.write(
        '''
        ##======================================##
        ## Bonds
        ##======================================##
        \n''')

    if not system.rigid:
        f.write("\t bond_style harmonic \n")
        for item in system.bond_coeff:
            f.write(item)

    if system.newBonds:
        f.write("\t bond_style harmonic \n")
        for item in system.bond_coeff:
            f.write(item)

    f.write(
        '''
        ##======================================##
        ## Angles
        ##======================================##
        \n''')

    if system.angles_interactions:
        f.write("angle_style harmonic \n")
        f.write(f"angle_coeff  1 {system.eps['angle']} 90\n")

    f.write(
    '''
    ##======================================##
    ## Setting up groups
    ##======================================##
    \n''')

    # should be changed for different vesicle sizes
    if system.rigid:
        f.write("group \t linked_vesicle \t type \t" + str(system.types['vesicle1']) + ":" + str(system.types['linker']) + "\n")
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

    if system.newBonds:
        f.write(f"fix cr_bond all bond/create 10 {system.types['synapsin']} {system.types['linker']} {intra_range + system.sigma['linker'] + system.sigma['synapsin']} {system.bondTypes['linker-synapsin']} "
                f"iparam 1 {system.types['synapsin_b']} jparam 1 {system.types['linker']} prob {create_prob} {seed}\n")

        #f.write(f"fix br_bond all bond/break 10 {system.bondTypes['linker-synapsin']} {(system.sigma['linker'] + system.sigma['synapsin'])*1} prob {break_prob} {seed}\n")

    if system.rigid:
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

    f.write(f"thermo {dump_step}\n")
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
    f.write(f"timestep {timestep}\n")
    f.write("run %i\n"%run_steps)
    f.close()

