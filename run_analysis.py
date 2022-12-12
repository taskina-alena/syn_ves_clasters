from cluster_analysis import vesicle_analysis
import argparse
import os

if __name__ == '__main__':

	scratch_loc = "../"
	parser = argparse.ArgumentParser(description="Script for Synaptic vesicle simulations.")
	parser.add_argument('--eps', '-eps', dest='eps', action='store', type=float, default=1.0, help='range of attraction')
	parser.add_argument('--num_v', '-num_v', dest='num_v', action='store', type=int, default=1, help='range of attraction')
	parser.add_argument('--num_l', '-num_l', dest='num_l', action='store', type=int, default=1, help='range of attraction')
	parser.add_argument('--real', '-real', dest='real', action='store', type=int, default=0, help='range of attraction')

	args = parser.parse_args()
	eps = args.eps
	num_vesicles = args.num_v
	num_linkers = args.num_l
	real = args.real


	file_pattern = f"2d_eps{eps:.4f}_nv{num_vesicles:n}_nl{num_linkers:n}_{real:n}"
	movie_name = scratch_loc + 'Results/Movie_%s.xyz'%file_pattern
	oligo_name = scratch_loc + 'Results/Oligo_%s.dat'%file_pattern


	va = vesicle_analysis(movie_name)
	if os.path.exists(movie_name):
		va = vesicle_analysis(movie_name)
		va.sv_cluster_size(oligo_name)

