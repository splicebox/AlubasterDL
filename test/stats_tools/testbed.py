import pybedtools
from exalu.data_preprocess.curate_data import curate_data

# curate_data(infer_set='Gencode', strand=True, mode='none', single_side_pad_len=0)
# curate_data('all_human_analysis', True, 'none', work_dir='/home/zhe29/Projects/eXAlu/data/all_human_analysis/known_datasets')

# a = {1: ['1','2','3','4','5','6']}
# import random
# b = {}
# for k, v in a.items():
#     random.seed(42)
#     random.shuffle(v)
#     print(v)
#     a[1] = v[0:2] * 3
#     b[0] = v[2:] * 3
# print(a)
# print(b)

# def t1():
#     print(a)

# a = 1000000
# t1()



# R: purine, = A or G
# Y: pyrimidine, = C or T
# M is A or C
# base_dict = {
#             'A': [1., 0., 0., 0.], 'a': [1., 0., 0., 0.],
#             'C': [0., 1., 0., 0.], 'c': [0., 1., 0., 0.],
#             'G': [0., 0., 1., 0.], 'g': [0., 0., 1., 0.],
#             'T': [0., 0., 0., 1.], 't': [0., 0., 0., 1.],
#             # 'N': [0.25, 0.25, 0.25, 0.25], 'n': [0.25, 0.25, 0.25, 0.25], 
#             'N': [0., 0., 0., 0.], 'n': [0., 0., 0., 0.], 
#             'R': [0.5, 0., 0.5, 0.], 'r': [0.5, 0., 0.5, 0.], 
#             'Y': [0., 0.5, 0., 0.5], 'y': [0., 0.5, 0., 0.5], 
#             'M': [0.5, 0.5, 0, 0], 'm': [0.5, 0.5, 0, 0]}

# def seq_encoding_onehot(seq):
#     import numpy as np
#     matrix = [base_dict[base] for base in seq]
#     return np.array(matrix)

# if __name__ == "__main__":
#     matrix = seq_encoding_onehot('AGTGAAGGTAAGGTCTACCAGGTTAGAT')
#     for i in matrix.T:
#         for j in i:
#             print(int(j), sep='', end='')
#         print('')

import numpy as np
# a = [
# [0.        ,  0.09758019,  0.05638206,  0.03606001],
# [0.00835928,  0.        ,  0.00122505, -0.02645838],
# [0.05264246,  0.12638748,  0.11498618,  0.        ],
# [0.00077811,  0.04577684,  0.        ,  0.02092943],
# [0.00302896,  0.134877  ,  0.01410973,  0.        ],
# [0.00224656,  0.10479909,  0.        ,  0.01092133],
# [0.        ,  0.08041632,  0.0413374 ,  0.01704729],
# [0.        ,  0.06007487,  0.01055062,  0.00825566],
# [0.        ,  0.01251188,  0.04790562,  0.00021157],
# [0.        ,  0.04836807,  0.03520691,  0.02105251]
# ]
# a = [
# 0.        ,  0.09758019,  0.05638206,  0.03606001,
# 0.00835928,  0.        ,  0.00122505, -0.02645838,
# 0.05264246,  0.12638748,  0.11498618,  0.        ,
# 0.00077811,  0.04577684,  0.        ,  0.02092943,
# 0.00302896,  0.134877  ,  0.01410973,  0.        ,
# 0.00224656,  0.10479909,  0.        ,  0.01092133,
# 0.        ,  0.08041632,  0.0413374 ,  0.01704729,
# 0.        ,  0.06007487,  0.01055062,  0.00825566,
# 0.        ,  0.01251188,  0.04790562,  0.00021157,
# 0.        ,  0.04836807,  0.03520691,  0.02105251
# ]

b = [[-0.02430326,  0.00332928,  0.        , -0.01447126],
 [ 0.02990192,  0.03854769,  0.05445281,  0.        ],
 [ 0.02308571,  0.02879864,  0.        , -0.00048441],
 [-0.00087896, -0.00280699,  0.        , -0.00825965],
 [ 0.03798553,  0.        ,  0.01937383,  0.02415267],
 [ 0.07931906, -0.0062674 ,  0.01854062,  0.        ],
 [-0.00485647,  0.        , -0.00604931, -0.01572374],
 [ 0.        ,  0.01758298,  0.00485659,  0.03130439],
 [ 0.00165585,  0.00225231,  0.0095605 ,  0.        ],
 [ 0.40617102, -0.00623754,  0.        , -0.02355921]]
b = np.array(b)
print(b.mean())
print(b.var())

b = [[-0.02430326,  0.00332928,  0.        , -0.01447126],
 [ 0.02990192,  0.03854769,  0.05445281,  0.        ],
 [ 0.02308571,  0.02879864,  0.        , -0.00048441],
 [-0.00087896, -0.00280699,  0.        , -0.00825965],
 [ 0.03798553,  0.        ,  0.01937383,  0.02415267],
 [ 0.07931906, -0.0062674 ,  0.01854062,  0.        ],
 [-0.00485647,  0.        , -0.00604931, -0.01572374],
 [ 0.        ,  0.01758298,  0.00485659,  0.03130439],
 [ 0.00165585,  0.00225231,  0.0095605 ,  0.        ],
 [ 0.00617102, -0.00623754,  0.        , -0.02355921]]
b = np.array(b)
print(b.mean())
print(b.var()) 
# m = sum(a)/len(a)
# v = sum([(i - m)**2 for i in a])/len(a)
# print(m)
# print(v)
# print()
# a = np.array(a)
# print(a)
# print(a.mean())
# print(a.var())