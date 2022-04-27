import sys, os
import argparse
import pdb
sys.path.append('./Models/')
sys.path.append('../src/ccscp/')
sys.path.append('../src/utils/')
# import astrobee_plot
from src.ccscp.cc_ocp import CCOCP
from Models.vehicle import Model
from arguments import add_arguments
# from astro_iss_plot import plot
# from astrobee_mc import monte_carlo
from Models.vehicle_plot import plot
import numpy as np

argparser = argparse.ArgumentParser(description='CARLA CILQR')
add_arguments(argparser)
args = argparser.parse_args()

# Initialize Astrobe model
m = Model(args)
N = m.N

# par               = dict()
# par['X_last']     = np.empty(shape=[m.n_x,       N  ])
# par['U_last']     = np.empty(shape=[m.n_u,       N-1])

# par['X_last'], par['U_last'] = m.initialize_trajectory(N)
# f_all, A_all, B_all    = m.compute_dynamics( par['X_last'], par['U_last'])
# Vars_all, Vars_dxu_all = m.propagate_variances(par['X_last'], par['U_last'],
#                                                        A_all, B_all)
# all_X, all_U, all_V = [], [], []
# all_X.append(par['X_last'])
# all_U.append(par['U_last'])
# all_V.append(Vars_all)


# Create chance-constrained problem
problem = CCOCP(m)

# Solve problem using CCSCP
problem.solve_ccscp(m)

# (Optional) Verify with Monte-Carlo
# X_sol, U_sol = problem.get_XU_solution_CCSCP(m)
# Xs_true, Us_true, nb_in_obs = monte_carlo(X_sol, U_sol, m, N_MC=1)

# Plot results
# print("X_last: {}".format(par['X_last']))
# print("X_last shape: {}".format(par['X_last'].shape))
# plot(np.stack(all_X), np.stack(all_U), np.stack(all_V), m)
plot(problem.all_X, problem.all_U, problem.all_V, m)
# plot(problem.all_X, problem.all_U, problem.all_V, m, Xs_true, Us_true)

# View solutions at each SCP iteration
# astrobee_plot.plot(problem.all_X, problem.all_U, problem.all_V, m)