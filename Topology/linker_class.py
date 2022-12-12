import numpy as np
import random
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


class make_linked_vesicles(object):

    def __init__(self, num, eps, sigma, side, rigid, inter_range, attr_with_rep, newBonds, linker, angles_interactions, ves_fraction, l_per_ves):

        self.rigid = rigid

        self.num = num
        self.eps = eps
        self.sigma = sigma
        self.inter_range = inter_range
        self.attr_with_rep = attr_with_rep
        self.newBonds = newBonds
        self.linker = linker
        self.angles_interactions = angles_interactions
        self.ves_fraction = ves_fraction
        self.l_per_ves = l_per_ves

        # shows whether this type is in the system or not (yet)
        self.types = {'linker': False, 'synapsin': False}
        for i in range(len(self.sigma['vesicle'])):
            self.types['vesicle' + str(i + 1)] = False

        self.side = side
        self.lattice_constant = self.sigma['synapsin'] * 2.5  # resolution of lattice
        self.Lx = self.side
        self.Ly = self.side
        self.Lz = 0.1

        self.lj_factor = 2 ** (1 / 6)

        self.coords = ["\nAtoms \n \n"]
        self.bonds = ["\nBonds \n \n"]
        self.angles = ["\nAngles \n \n"]

        self.bond_coeff = []

        # m is number of individual and k is number of "molecules"
        self.numAll = 0
        self.k = 0

        self.numBonds = 0
        self.bondId = 0
        self.bondTypes = {}
        self.numBondTypes = 0

        self.numAngles = 0
        self.angleId = 0
        if self.rigid:
            self.numAngleTypes = 0
        else:
            self.numAngleTypes = 1

        self.numTypes = 0
        self.type_mass_list = []

        # return coordinates of lattice
        self.num_x = int(self.Lx / self.lattice_constant)
        self.num_y = int(self.Ly / self.lattice_constant)

        def init_lattice_2d():
            lattice_out = np.zeros(shape=(self.num_x * self.num_y, 3))
            counter = 0
            for i in range(self.num_y):
                for j in range(self.num_x):
                    lattice_out[counter, 0] = self.lattice_constant * i - self.lattice_constant * self.num_x * 0.5
                    lattice_out[counter, 1] = self.lattice_constant * j - self.lattice_constant * self.num_y * 0.5
                    lattice_out[counter, 2] = 0.0
                    counter += 1
            return lattice_out

        self.lattice_sites = init_lattice_2d()

        # assert self.num_synapsin + self.num_linkers + self.num_vesicles <= self.lattice_sites.shape[0], "ERROR: not enough lattice sites, increase box size."

        # numerate lattice coordinates
        lattice_inds = np.arange(self.lattice_sites.shape[0])

        # choose randomly point at the lattice (its number)
        self.start_inds = {}
        lattice_mask = np.ones(len(lattice_inds), bool)
        lattice_mask_ves = np.zeros(len(lattice_inds), bool)

        vesicle_size = (max(self.sigma['vesicle']) + 2 * self.sigma['linker']) / self.lattice_constant
        vesicle_step = round(vesicle_size + 0.5) * 2
        num_x_vesicles = self.num_x // vesicle_step
        num_y_vesicles = self.num_y // vesicle_step
        assert (num_y_vesicles > 0)
        assert (num_x_vesicles > 0)
        # shift_x = ((self.num_x - 1) % vesicle_step) // 2 + 1
        # shift_y = ((self.num_y -1) % vesicle_step) // 2 + 1
        for i in range(num_y_vesicles):
            for j in range(num_x_vesicles):
                lattice_mask_ves[(int((vesicle_step / 2 + i * vesicle_step) * self.num_x + (
                            vesicle_step / 2 + j * vesicle_step)))] = True

        assert (sum(lattice_mask_ves) >= self.num['vesicle'])

        seed = 100
        np.random.seed(seed)

        for i in range(len(self.sigma['vesicle'])):
            n = round( self.num['vesicle'] * self.ves_fraction[ 'vesicle_' + str(sigma['vesicle'][i]) ])
            self.start_inds['vesicle' + str(i + 1)] = np.random.choice(lattice_inds[lattice_mask_ves], size=n,
                                                                       replace=False)
            lattice_mask_ves[self.start_inds['vesicle' + str(i + 1)]] = False
            for ind in self.start_inds['vesicle' + str(i + 1)]:
                for i in range(round(-vesicle_size - 0.5), round(vesicle_size + 0.5)):
                    for j in range(round(-vesicle_size - 0.5), round(vesicle_size + 0.5)):
                        ind_occupied = int(ind + i * self.num_x + j)
                        lattice_mask[ind_occupied] = False
        # which lattice inds remained free

        self.start_inds['synapsin'] = np.random.choice(lattice_inds[lattice_mask], size=self.num['synapsin'],
                                                       replace=False)

        # for i in range(len(self.vesicle_sigma)):
        # plt.scatter(self.lattice_sites[self.start_inds['vesicle' + str(i+1)]][:, 0], 1+self.lattice_sites[self.start_inds['vesicle' + str(i+1)]][:, 1 ], s=self.vesicle_sigma[i])
        # plt.scatter( self.lattice_sites[self.start_inds['synapsin']][:, 0 ], self.lattice_sites[self.start_inds['synapsin']][:, 1 ], c ='g', s=0.5)
        # plt.scatter( self.lattice_sites[self.start_inds['linker']][:, 0 ], -1+self.lattice_sites[self.start_inds['linker']][:, 1 ], c ='r', s=5)
        # plt.show()

        # test_list = list(self.start_inds['vesicle']) + list(self.start_inds['synapsin']) + list(self.start_inds['linker'])
        # assert len(test_list) == self.num_vesicles + self.num_linkers + self.num_synapsin

    def interactions(self, global_cutoff=3):
        output = [f"\t pair_style cosine/squared {global_cutoff} \n"]

        for i in self.types.keys():
            for j in self.types.keys():
                if (self.types[i] <= self.types[j]) and self.types[i] and self.types[j]:

                    if (self.newBonds and (i == 'synapsin_b' and 'j == synapsin_b')) or (j == 'linker' and i == 'synapsin')\
                            or (not self.linker and (i[:8]=='synapsin' and j[:8]=='synapsin')):
                        epsilon = self.eps['attr']
                        add_cutoff = self.inter_range
                        if self.attr_with_rep:
                            potential = 'wca'
                        else:
                            potential = ''

                    else:
                        epsilon = self.eps['rep']
                        add_cutoff = 0
                        potential = 'wca'

                    if i[:7] == 'vesicle':
                        sigma1 = self.sigma['vesicle'][int(i[7]) - 1]
                    elif j[:8] == 'synapsin':
                        sigma1 = self.sigma['synapsin']
                    else:
                        sigma1 = self.sigma[i]

                    if j[:7] == 'vesicle':
                        sigma2 = self.sigma['vesicle'][int(j[7]) - 1]
                    elif j[:8] == 'synapsin':
                        sigma2 = self.sigma['synapsin']
                    else:
                        sigma2 = self.sigma[j]

                    output.append(
                        f"\t pair_coeff {self.types[i]} {self.types[j]} {epsilon} {self.lj_factor * (sigma1 + sigma2)} {add_cutoff + self.lj_factor * (sigma1 + sigma2)} {potential} \n")  ## eps sigma cutoff

        return output

    def make_single_particle(self, p_type):
        if not self.types[p_type]:
            self.numTypes += 1
            self.types[p_type] = self.numTypes
            self.type_mass_list.append([self.numTypes, 1])

        for i in range(len(self.start_inds[p_type])):  ## go through chains
            indBuf = self.start_inds[p_type][i]
            self.k += 1  # counter of molecules
            self.numAll += 1  # counter of particles
            x = self.lattice_sites[indBuf, 0]
            y = self.lattice_sites[indBuf, 1]
            z = self.lattice_sites[indBuf, 2]

            self.coords.append(
                "\t " + str(self.numAll) + " " + str(self.k) + " " + str(self.types[p_type]) + " " + str(
                    x) + " " + str(y) + " " + str(z) + " 0 0 0 \n")

    def dist_on_circle(self, sigma, num):
        cf = 2 * np.pi * sigma
        max_num = np.floor(cf / (2*self.sigma[linked]))
        if num >= max_num:
            raise RuntimeError('Too many linkers on circumference, reduce number of linkers.')
        linker_angle = 2.0 * np.pi / num
        positions = np.zeros(shape=(num, 2))
        for n in range(num):
            angle_arg = n * linker_angle
            positions[n, 0] = sigma * np.cos(angle_arg)
            positions[n, 1] = sigmas[i] * np.sin(angle_arg)
        return positions

    def make_linked_vesicles(self):

        def dist_on_circle(sigma, num):
            cf = 2 * np.pi * sigma
            max_num = np.floor(cf / self.sigma[linked])
            if num >= max_num:
                raise RuntimeError('Too many linkers on circumference, reduce number of linkers.')
            linker_angle = 2.0 * np.pi / num
            positions = np.zeros(shape=(num, 2))
            for n in range(num):
                angle_arg = n * linker_angle
                positions[n, 0] = sigma * np.cos(angle_arg)
                positions[n, 1] = sigma * np.sin(angle_arg)
            return positions

        for t in ['vesicle' + str(i + 1) for i in range(len(self.sigma['vesicle']))]:
            if not self.types[t]:
                self.numTypes += 1
                self.types[t] = self.numTypes
                self.type_mass_list.append([self.numTypes, 1])

        if self.linker:
            linked = 'linker'
        else:
            linked = 'synapsin'

        if linked=='linker':
            if not self.types['linker']:
                self.numTypes += 1
                self.types['linker'] = self.numTypes
                self.type_mass_list.append([self.numTypes, 1])
        else:
            for i in range(self.num['vesicle']):
                self.numTypes += 1
                self.types['synapsin' + str(i)] = self.numTypes
                self.type_mass_list.append([self.numTypes, 1])

        for sigma_index in range(len(self.sigma['vesicle'])):

            if not self.rigid:
                self.numBondTypes += 1
                self.bondTypes['vesicle' + str(sigma_index + 1) + '-' + linked] = self.numBondTypes
                self.bond_coeff.append(
                    "\t bond_coeff " + str(self.numBondTypes) + " " + str(self.eps['bond']) + " " + str(
                        self.sigma['vesicle'][sigma_index] + self.sigma[linked]) + " \n")
                # bond_type energy eq_distance

            position_on_cirlce = dist_on_circle(self.sigma['vesicle'][sigma_index], self.l_per_ves['vesicle_' + str(self.sigma['vesicle'][sigma_index])])

            for i in range(len(self.start_inds['vesicle' + str(sigma_index + 1)])):  ## go through chains
                indBuf = self.start_inds['vesicle' + str(sigma_index + 1)][i]  # index of lattice occupied by vesicle

                ## counter of vesicles
                self.k += 1

                ## coordinates of vesicle position
                x = self.lattice_sites[indBuf, 0]
                y = self.lattice_sites[indBuf, 1]
                z = self.lattice_sites[indBuf, 2]

                bonded_num = 0

                for n in range(self.l_per_ves['vesicle_' + str(self.sigma['vesicle'][sigma_index])]+1):
                    ## counter of all particles
                    self.numAll += 1
                    '''
                    write info about particles
                    particle_number  molecular_number inner_mol_number charge x y z 0 0 0
                    '''
                    ## place vesicle
                    if n == 0:
                        self.coords.append("\t " + str(self.numAll) + " " + str(self.k) + " " + str(
                            self.types['vesicle' + str(sigma_index + 1)]) + " " + str(x) + " " + str(y) + " " + str(
                            z) + " 0 0 0 \n")
                        bonded_num = self.numAll


                    if n > 0:

                        if linked=='linker':
                            self.coords.append(
                                "\t " + str(self.numAll) + " " + str(self.k) + " " + str(self.types['linker']) + " " +
                                str(x - position_on_cirlce[n - 1, 0]) + " " + str(y - position_on_cirlce[n - 1, 1]) + " " + str(z) + " 0 0 0 \n")
                        else:
                            self.coords.append(
                                "\t " + str(self.numAll) + " " + str(self.k) + " " + str(self.types['synapsin' + str(self.k-1)]) + " " +
                                str(x - position_on_cirlce[n - 1, 0]) + " " + str(
                                    y - position_on_cirlce[n - 1, 1]) + " " + str(z) + " 0 0 0 \n")
                        # particle_number molecular_number type_number x y z 0 0 0

                        if not self.rigid:
                            self.bondId += 1
                            self.bonds.append("\t" + str(self.bondId) + " " + str(self.numBondTypes) + " " + str(
                                self.numAll) + " " + str(bonded_num) + " \n")
                            self.numBonds += 1
                            # bond_number bond_type particle_number particle_bonded_number

                        if self.angles_interactions:
                            self.angleId += 1
                            if n > 1:
                                self.angles.append("\t " + str(self.angleId) + " 1 " + str(self.numAll) + " " + str(
                                    bonded_num) + " " + str(self.numAll - 1) + "\n")
                            else:
                                self.angles.append("\t " + str(self.angleId) + " 1 " + str(self.numAll) + " " + str(
                                    bonded_num) + " " + str(self.numAll + 3) + "\n")
                            self.numAngles += 1
                            # angle_number angle_type linker1 vesicle linker2


        if self.newBonds:
            self.numBondTypes += 1
            self.bondTypes['linker-synapsin'] = self.numBondTypes
            self.bond_coeff.append(
                f"\t bond_coeff {self.numBondTypes} {self.eps['bond']} {self.sigma['linker'] + self.sigma['synapsin']} \n")
            self.numTypes += 1
            self.types['synapsin_b'] = self.numTypes
            self.type_mass_list.append([self.numTypes, 1])
