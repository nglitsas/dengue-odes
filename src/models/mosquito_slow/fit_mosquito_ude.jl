module CapacityNN

"""
This module optimizes the weights on a neural network that takes temperature 
and precipitation in order to minimize the SSE between an ODE system
(defined in model_dynamics.jl) that uses the output of the NN to model
the carrying capacity of the mosquito population at time t (C(t)). 
"""
