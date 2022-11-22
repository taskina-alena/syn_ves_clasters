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

        self.types =  {'vesicle' : False, 'linker' : False, 'synapsin' : False}

        self.side = side
        self.lattice_constant = self.synapsin_sigma * 2.5
        self.Lx = self.side
        self.Ly = self.side
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

        #assert self.num_synapsin + self.num_linkers + self.num_vesicles <= self.lattice_sites.shape[0], "ERROR: not enough lattice sites, increase box size."

        # numerate lattice coordinates
        lattice_inds = np.arange(self.lattice_sites.shape[0])

        # choose randomly point at the lattice (its number)
        self.start_inds = {}
        lattice_mask = np.ones(len(lattice_inds), bool)

        lattice_inds_vesicle = []
        vesicle_size = (self.vesicle_sigma/self.lattice_constant)
        vesicle_step = round(vesicle_size + 0.5) * 2
        num_x_vesicles = self.num_x // vesicle_step
        num_y_vesicles = self.num_y // vesicle_step
        assert(num_y_vesicles>0)
        assert (num_x_vesicles > 0)
        #shift_x = ((self.num_x - 1) % vesicle_step) // 2 + 1
        #shift_y = ((self.num_y -1) % vesicle_step) // 2 + 1
        for i in range(num_y_vesicles):
            for j in range(num_x_vesicles):
                lattice_inds_vesicle.append(int((vesicle_step/2 + i * vesicle_step) * self.num_x + (vesicle_step/2 + j * vesicle_step)))

        assert(len(lattice_inds_vesicle) >= self.num_vesicles)

        self.start_inds['vesicle'] = np.random.choice(lattice_inds_vesicle, size = self.num_vesicles, replace = False)


        for ind in self.start_inds['vesicle']:
            for i in range(round(-vesicle_size-0.5), round(vesicle_size+0.5)):
                for j in range(round(-vesicle_size-0.5), round(vesicle_size+0.5)):
                    ind_occupied = int(ind + i * self.num_x + j)
                    lattice_mask[ind_occupied] = False
        # which lattice inds remained free

        self.start_inds['synapsin'] = np.random.choice(lattice_inds[lattice_mask], size = self.num_synapsin, replace = False )

        #plt.scatter(self.lattice_sites[self.start_inds['vesicle']][:, 0], 1+self.lattice_sites[self.start_inds['vesicle']][:, 1 ], s=40)
        #plt.scatter( self.lattice_sites[self.start_inds['synapsin']][:, 0 ], self.lattice_sites[self.start_inds['synapsin']][:, 1 ], c ='g', s=5)
        #plt.scatter( self.lattice_sites[self.start_inds['linker']][:, 0 ], -1+self.lattice_sites[self.start_inds['linker']][:, 1 ], c ='r', s=5)
        #plt.show()


        #test_list = list(self.start_inds['vesicle']) + list(self.start_inds['synapsin']) + list(self.start_inds['linker'])
        #assert len(test_list) == self.num_vesicles + self.num_linkers + self.num_synapsin
    def make_single_particle(self, p_type):
        if self.types[p_type]:
            pass
        else:
            self.numTypes += 1
            self.types[p_type] = self.numTypes
            self.type_mass_list.append([self.numTypes, 1])

        for i in range(len(self.start_inds[p_type])):  ## go through chains
            indBuf = self.start_inds[p_type][i]
            self.k += 1 # counter of molecules
            self.numAll += 1 # counter of particles
            x = self.lattice_sites[indBuf, 0]
            y = self.lattice_sites[indBuf, 1]
            z = self.lattice_sites[indBuf, 2]

            self.coords.append("\t " + str(self.numAll) + " " + str(self.k) + " " + str(self.types[p_type]) + " 0 " + str(x) + " " + str(y) + " " + str(z) + " 0 0 0 \n")

    def make_linked_vesicles(self):

        for t in ['vesicle', 'linker']:
            if self.types[t]:
                pass
            else:
                self.numTypes += 1
                self.types[t] = self.numTypes
                self.type_mass_list.append([self.numTypes, 1])

        for i in range(len(self.start_inds['vesicle'])):  ## go through chains
            indBuf = self.start_inds['vesicle'][i] # index of lattice occupied by vesicle

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
                particle_number  molecular_number inner_mol_number charge x y z 0 0 0
                '''
                ## place vesicle
                if n == 0:
                    self.coords.append(
                        "\t " + str(self.numAll) + " " + str(self.k) + " " + str(self.types['vesicle']) + " 0 " +
                        str(x) + " " + str(y) + " " + str(z) + " 0 0 0 \n")
                    bonded_num = self.numAll

                radial_dist = self.vesicle_sigma + self.linker_sigma

                position_on_cirlce = np.array([[radial_dist,0], [0,-radial_dist], [-radial_dist,0], [0,radial_dist]])

                if n > 0:
                    self.coords.append(
                        "\t " + str(self.numAll) + " " + str(self.k) + " " +  str(self.types['linker']) + " 0 " +
                        str(x - position_on_cirlce[n-1,0]) + " " + str(y - position_on_cirlce[n-1,1]) + " " + str(z) + " 0 0 0 \n")
                    #particle_number molecular_number inner_mol_number charge x y z 0 0 0

                    self.bondId += 1
                    self.bonds.append("\t" + str(self.bondId) + " 1 " + str(self.numAll) + " " + str(bonded_num) + " \n")
                    self.numBonds += 1
                    #bond_number 1 particle_number particle_bonded_number


                    self.angleId += 1
                    if n>1:
                        self.angles.append("\t " + str(self.angleId) + " 1 " + str(self.numAll) + " " + str(bonded_num) + " " + str(self.numAll - 1) + "\n")
                    else:
                        self.angles.append("\t " + str(self.angleId) + " 1 " + str(self.numAll) + " " + str(bonded_num) + " " + str(self.numAll + 3) + "\n")
                    self.numAngles += 1
                    #angle_number 1 linker_1 vesicle linker_2


