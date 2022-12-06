from Topology.linker_class import make_linked_vesicles
import os
import numpy as np
import input_func

# TODO try attractive interactions between synapsin
# TODO introduce clastering of synapsin before vesicles
# TODO count the sizes of clusters
# TODO make different numbers of linkers per vesicle
# TODO create bonds

if __name__ == "__main__":

    sim_seed = 5
    dump_step = 1000
    eq_steps = 1e3
    run_steps = 1e6
    create_prob = 1
    break_prob = 0
    intra_range = 0

    rigid_molecules = True
    eps = {'attr': 6, 'rep': 1, 'bond': 1, 'angle': 100}
    inter_range = 0.25
    newBonds = True
    attr_with_rep = True
    side_box = 200
    num = {'vesicle':100, 'synapsin':4000, 'linker':400}
    #side_box = 600
    #num = {'vesicle':900, 'synapsin':36000, 'linker':3600}
    sigma = {'vesicle': [3], 'linker': 0.5, 'synapsin': 0.5}

    system = make_linked_vesicles(num, eps, sigma, side_box, rigid_molecules, inter_range, attr_with_rep, newBonds)
    system.make_single_particle('synapsin')
    system.make_linked_vesicles()

    ves_sizes = ''
    for i in range(len(sigma['vesicle'])):
        ves_sizes+=str(sigma['vesicle'][i])
        ves_sizes+='_'
    filePattern = f"change_attr_ves_sizes{ves_sizes}_eps_attr{eps['attr']}_eps_rep{eps['rep']}_range{inter_range}_attr_with_rep{attr_with_rep}_angle{eps['angle']}_n_synapsin{num['synapsin']}_n_vesicles{num['vesicle']}_createProb{create_prob}_breakProb{break_prob}_bond{eps['bond']}"
    folderPattern = f"Results_change_attr_ves_sizes{ves_sizes}_eps_attr{eps['attr']}_eps_rep{eps['rep']}_range{inter_range}_attr_with_rep{attr_with_rep}_angle{eps['angle']}_n_synapsin{num['synapsin']}_n_vesicles{num['vesicle']}_createProb{create_prob}_breakProb{break_prob}_bond{eps['bond']}"

    if not os.path.exists(folderPattern):
        os.makedirs(folderPattern)

    configName = "Input/Configuration/Config_%s.dat" % filePattern

    header = ["LAMMPS Description \n \n",
              "\t " + str(system.numAll) + " atoms \n \t " + str(system.numBonds) + " bonds \n \t " +
              str(system.numAngles) + " angles \n \t 0 dihedrals \n \t 0 impropers \n",
              "\n \t " + str(system.numTypes) + " atom types \n \t " + str(
                  system.numBondTypes) + " bond types \n \t" + str(
                  system.numAngleTypes) + " angle types \n \t 0 dihedral types \n \t 0 improper types \n",
              "\n \t 10 extra bond per atom  \n \t 10000 extra special per atom",
              "\n \t " + str(-system.Lx * 0.5) + " " + str(system.Lx * 0.5) + " xlo xhi\n \t",
              str(-system.Ly * 0.5) + " " + str(system.Ly * 0.5) + " ylo yhi \n \t",
              str(-system.Lz * 0.5) + " " + str(system.Lz * 0.5) + " zlo zhi\n", "\nMasses \n \n"]

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
        for item in system.angles:
            f.write(item)

    f.close()

    input_func.write_in_script(system, folderPattern, filePattern, configName, dump_step, run_steps, create_prob, break_prob, intra_range)
















