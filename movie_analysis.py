from tools import lmp_tools as lt
import numpy as np
import pandas as pd

data={}
for eps_attr in range(1, 11):
    data[eps_attr] = {}
    for d in range(10):
        density = np.round(0.08 + d*0.1, 2)
        movie_name = f"Results_span/Results_without_linker_ves_sizes3_4_5_density{density}_frac0.25_0.5_0.25_eps_attr{eps_attr}_eps_rep1_range0.25_attr_with_repTrue_n_vesicles100/Movie_without_linker_ves_sizes3_4_5_frac0.25_0.5_0.25_density{density}_eps_attr{eps_attr}_eps_rep1_range0.25_attr_with_repTrue_n_vesicles100.xyz"
        movie_format = ['id', 'type', 'mol', 'x', 'y', 'z']
        read_lines = lt.header_lines(movie_name)
        try:
            work_movie = lt.split_movie(movie_name, movie_format, read_lines)
        except:
            data[eps_attr][density] = 'error'
        else:
            data[eps_attr][density] = len(work_movie)
            config = work_movie[0]
            coords = config[:, 3:6]
            r = [ i in [1] for i in config[:, 1] ]
            lines_sv = np.where(r)

print(data)
print('final')