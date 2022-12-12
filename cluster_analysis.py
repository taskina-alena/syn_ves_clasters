import warnings
from tools import lmp_tools as lt
import freud as fr
import numpy as np
from matplotlib import pyplot as plt
import time

start_time = time.time()
import os


class vesicle_analysis(object):
    def __init__(self, movie_name):
        self.movie_name = movie_name
        self.vesicle_type = [1, 2, 3]

    def sv_cluster_size(self, oligo_file_name):
        '''
        This function computes cluster properties of synaptic vesicle clusters.

        '''
        movie_name = self.movie_name

        sv_type = self.vesicle_type

        r_cut = 1.5 * 8  ## cut-off for cluster computes
        max_k_mer = 400

        movie_format = ['id', 'type', 'mol', 'x', 'y', 'z']  # structure of the movie data to be loaded
        box_data = np.genfromtxt(movie_name, skip_header=5, max_rows=3)  # grab box data from movie file x_min-x_max y_min-y_max z_min-z_max

        # workingConfig = np.array(pd.read_csv(configName, sep="\s+", names=self.configFormat, header=None, dtype=object))

        read_lines = lt.header_lines(movie_name)  ## select or slice from read_lines to load single frames ToDo to be implemented
        work_movie = lt.split_movie(movie_name, movie_format, read_lines)

        box_lens = box_data[:, 1] - box_data[:, 0]
        box = fr.box.Box(box_lens[0], box_lens[1], is2D=True)  # define the box that the freud package can use
        # box.periodic_z = False

        frames_to_process = [-10, -9, -8, -7, -6, -5, -4, -3, -2, -1]

        oligo_dist_time = np.zeros(shape=(len(frames_to_process), max_k_mer + 2))

        for frame in frames_to_process:  # frame to be looked into

            config = work_movie[frame]
            coords = config[:, 3:6]

            r = [i in sv_type for i in config[:, 1]]
            lines_sv = np.where(r)[0]

            coords_sv = coords[lines_sv]

            aq = fr.AABBQuery(box, coords_sv)  ## Neighbour query for cluster computation
            # query_result = aq.query(pore_coords, dict(r_max=r_cut))
            # n_list = query_result.toNeighborList()

            cl = fr.cluster.Cluster()  ## create cluster compute object
            cl.compute(aq, neighbors={'r_max': r_cut})
            num_clusters = cl.num_clusters

            cl_prop = fr.cluster.ClusterProperties()
            cl_prop.compute(aq, cl.cluster_idx)

            oligo_dist = np.zeros(shape=max_k_mer + 1)
            mean_oligo_size = np.mean(cl_prop.sizes)
            oligo_dist[0] = mean_oligo_size

            size_dist = np.unique(cl_prop.sizes, return_counts=True)

            np.testing.assert_array_equal(size_dist[0], np.sort(size_dist[0]))
            assert np.max(size_dist[0]) <= max_k_mer
            count_clust = 0
            for i in size_dist[0]:
                oligo_dist[i] = size_dist[1][count_clust]
                count_clust += 1

            max_slice = max_k_mer + 2
            oligo_dist_time[frame, 0] = frame
            oligo_dist_time[frame, 1:max_slice] = oligo_dist

            # fig, ax = plt.subplots(1, 1, figsize=(9, 6))
            # for cluster_id in range(cl.num_clusters):
            #     cluster_system = fr.AABBQuery(aq.box, aq.points[cl.cluster_keys[cluster_id]])
            #     cluster_system.plot(ax=ax, s=100, label=f"Cluster {cluster_id}")
            #     print(
            #         f"There are {len(cl.cluster_keys[cluster_id])} points in cluster {cluster_id}."
            #     )
            # print(mean_oligo_size)
            # plt.show()

            np.savetxt(oligo_file_name, oligo_dist_time, delimiter=' ')


if __name__ == '__main__':
    for eps_attr in range(1, 7):
        for d in range(6):
            density = np.round(0.08 + d * 0.1, 2)
            file_pattern = f"without_linker_ves_sizes3_4_5_density{density}_frac0.25_0.5_0.25_eps_attr{eps_attr}_eps_rep1_range0.25_attr_with_repTrue_n_vesicles100"
            movie_pattern = f"without_linker_ves_sizes3_4_5_frac0.25_0.5_0.25_density{density}_eps_attr{eps_attr}_eps_rep1_range0.25_attr_with_repTrue_n_vesicles100"
            scratch_loc = 'Results_span/Results_%s/' % file_pattern

            movie_name = scratch_loc + 'Movie_%s.xyz' % movie_pattern
            oligo_name = scratch_loc + 'Oligo_%s.dat' % file_pattern

            print(movie_name)
            print(os.path.exists(movie_name))

            if os.path.exists(movie_name):
                va = vesicle_analysis(movie_name)
                va.sv_cluster_size(oligo_name)

            print("--- %s seconds ---" % (time.time() - start_time))

