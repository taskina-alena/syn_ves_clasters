import numpy as np
import random
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


class make_linked_vesicles(object):


    def __init__(self, num_synapsin, num_linkers, num_vesicles, side, sigma):

        self.num_synapsin = num_synapsin
        self.num_linkers = num_linkers
        self.num_vesicles = num_vesicles

        self.vesicle_sigma = sigma["vesicle"]  # size of vesicle
        self.linker_sigma = sigma["linker"]  # size of synapsin
        self.synapsin_sigma = sigma["synapsin"]

        self.side = side
        self.lattice_constant = 2.5*self.vesicle_sigma
        self.Lx = self.side*self.lattice_constant
        self.Ly = self.side*self.lattice_constant
        self.Lz = 0.1

        self.lj_factor = 2**(1/6)

        self.coords = ["\nAtoms \n \n"]
        self.bonds = ["\nBonds \n \n"]
        self.angles = ["\nAngles \n \n"]

        # m is number of individual and k is number of "molecules"
        self.numAll = 0
        self.k = 0

        self.numBonds = 0
        self.bondId = 0

        self.numAngles = 0
        self.angleId = 0

        self.numTypes = 0
        self.type_mass_list =[]

        # return coordinates of lattice
        def init_lattice_2d(lat_con):
            num_x = int(self.Lx/lat_con)
            num_y = int(self.Ly/lat_con)

            lattice_out = np.zeros(shape=(num_x*num_y, 3))
            counter = 0
            for i in range(num_x):
                for j in range(num_y):
                    lattice_out[counter, 0] = lat_con * i - lat_con * num_x * 0.5
                    lattice_out[counter, 1] = lat_con * j - lat_con * num_y * 0.5
                    lattice_out[counter, 2] = 0.0
                    counter += 1
            return lattice_out

        self.lattice_sites = init_lattice_2d(self.lattice_constant)

        assert self.num_synapsin + self.num_linkers + self.num_vesicles <= self.lattice_sites.shape[0], "ERROR: not enough lattice sites, increase box size."

        # numerate lattice coordinates
        lattice_inds = np.arange(self.lattice_sites.shape[0])

        # choose randomly point at the lattice (its number)
        self.start_inds_vesicles = np.random.choice(lattice_inds, size = self.num_vesicles, replace = False)

        # which lattice inds remained free
        self.lattice_inds_remain = np.array( [item for item in lattice_inds if item not in set(list(self.start_inds_vesicles))] )

        self.start_inds_synapsin = np.random.choice( self.lattice_inds_remain, size = self.num_synapsin, replace = False )

        self.lattice_inds_remain = np.array([item for item in lattice_inds if item not in set(list(self.start_inds_vesicles) +  list(self.start_inds_synapsin))])

        self.start_inds_linker = np.random.choice(self.lattice_inds_remain, size=self.num_linkers, replace=False)

        # plt.scatter( self.lattice_sites[self.start_inds_l1][:, 0 ], 1+self.lattice_sites[self.start_inds_l1][:, 1 ])
        # plt.scatter( self.lattice_sites[self.start_inds_l2][:, 0 ], self.lattice_sites[self.start_inds_l2][:, 1 ], c ='g')
        # plt.scatter( self.lattice_sites[self.start_inds_p][:, 0 ], -1+self.lattice_sites[self.start_inds_p][:, 1 ], c ='r')
        # plt.show()
        #
        # import sys
        # sys.exit()
        test_list = list(self.start_inds_vesicles) + list(self.start_inds_synapsin) + list(self.start_inds_linker)
        assert len(test_list) == self.num_vesicles + self.num_linkers + self.num_synapsin

    def make_linked_vesicles(self):

        particles_list = [1,2]
        particle_types = len(particles_list)
        self.numTypes += len(np.unique(particles_list))

        # for i in range(particle_types):
        #     to_add = [particles_list[i], 1]
        #     if to_add not in self.type_mass_list:
        #         self.type_mass_list.append(to_add)
        self.type_mass_list.append([1,1])
        self.type_mass_list.append([2,1])

        for i in range(len(self.start_inds_vesicles)):  ## go through chains
            indBuf = self.start_inds_vesicles[i] # index of lattice occupied by vesicle

            ## counter of vesicles
            self.k += 1

            ## coordinates of vesicle position
            x = self.lattice_sites[indBuf, 0]
            y = self.lattice_sites[indBuf, 1]
            z = self.lattice_sites[indBuf, 2]

            # print(x,y,z)
            bonded_num = 0
            for n in range(5):
                ## counter of all particles
                self.numAll += 1
                '''
                write info about particles
                particle_number  molecular_number inner_mol_number charge x y z 0 0 0 ???
                '''
                ## place vesicle
                if n == 0:
                    self.coords.append(
                        "\t " + str(self.numAll) + " " + str(self.k) + " " + str(particles_list[0]) + " 0 " +
                        str(x) + " " + str(y) + " " + str(z) + " 0 0 0 \n")
                    bonded_num = self.numAll

                radial_dist = 0.5 * self.lj_factor * self.vesicle_sigma
                position_on_cirlce = np.array([[radial_dist,0], [0,-radial_dist], [-radial_dist,0], [0,radial_dist]])

                if n > 0:
                    self.coords.append(
                        "\t " + str(self.numAll) + " " + str(self.k) + " " +  str(particles_list[1]) + " 0 " +
                        str(x - position_on_cirlce[n-1,0]) + " " + str(y - position_on_cirlce[n-1,1]) + " " + str(z) + " 0 0 0 \n")

                    #self.bondId += 1
                    #self.bonds.append("\t" + str(self.bondId) + " 1 " + str(self.numAll) + " " + str(bonded_num) + " \n")
                    #self.numBonds += 1

                    # self.angleId += 1
                    # self.angles.append("\t " + str(self.angleId) + " 1 " + str(self.numAll - 1) + " " + str(
                    #     self.numAll - 2) + " " + str(self.numAll) + "\n")
                    # self.numAngles += 1

    def make_synapsin(self):

        particles_list = [3]  ## first one is central body
        particle_types = len(particles_list)
        self.numTypes += particle_types
        self.type_mass_list.append([3, 1])

        for i in range(len(self.start_inds_synapsin)):  ## go through chains
            indBuf = self.start_inds_synapsin[i]
            self.k += 1 # counter of molecules
            self.numAll += 1 # counter of particles
            x = self.lattice_sites[indBuf, 0]
            y = self.lattice_sites[indBuf, 1]
            z = self.lattice_sites[indBuf, 2]

            self.coords.append("\t " + str(self.numAll) + " " + str(self.k) + " " + '3' + " 0 " + str(x) + " " + str(y) + " " + str(z) + " 0 0 0 \n")

    def make_vesicles(self):

        particles_list = [1]  ## first one is central body
        particle_types = len(particles_list)
        self.numTypes += particle_types
        self.type_mass_list.append([1, 1])

        for i in range(len(self.start_inds_vesicles)):  ## go through chains
            indBuf = self.start_inds_vesicles[i]
            self.k += 1  # counter of molecules
            self.numAll += 1  # counter of particles
            x = self.lattice_sites[indBuf, 0]
            y = self.lattice_sites[indBuf, 1]
            z = self.lattice_sites[indBuf, 2]

            self.coords.append("\t " + str(self.numAll) + " " + str(self.k) + " " + '1' + " 0 " + str(x) + " " + str(y) + " " + str(z) + " 0 0 0 \n")

    def make_linkers(self):

        particles_list = [2]  ## first one is central body
        particle_types = len(particles_list)
        self.numTypes += particle_types
        self.type_mass_list.append([2, 1])

        for i in range(len(self.start_inds_linker)):  ## go through chains
            indBuf = self.start_inds_linker[i]
            self.k += 1  # counter of molecules
            self.numAll += 1  # counter of particles
            x = self.lattice_sites[indBuf, 0]
            y = self.lattice_sites[indBuf, 1]
            z = self.lattice_sites[indBuf, 2]

            self.coords.append("\t " + str(self.numAll) + " " + str(self.k) + " " + '2' + " 0 " + str(x) + " " + str(y) + " " + str(z) + " 0 0 0 \n")