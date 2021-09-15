from GaussianInput import *
from SynapticConstants import *
from ConnectivityMatrix import *
from misc import *

import numpy as np
from math import pi

class RingPING:

    def __init__(self, nr_excit, nr_inhibit, simulation_time, dt):
        self.nr_excit = nr_excit
        self.nr_inhibit = nr_inhibit
        self.nr_neurons = self.nr_excit + self.nr_inhibit

        # timescale of recover variable u
        self.a = [0.02 for _ in range(self.nr_excit)] + [0.1 for _ in range(self.nr_inhibit)]
        # sensitivity of  u to sub-threshold oscillations
        self.b = [0.2 for _ in range(self.nr_excit)] + [0.2 for _ in range(self.nr_inhibit)]
        # membrane voltage after spike (reset)
        self.c = [-65 for _ in range(self.nr_excit)] + [-65 for _ in range(self.nr_inhibit)]
        # spike reset of recover variable u
        self.d = [8 for _ in range(self.nr_excit)] + [2 for _ in range(self.nr_inhibit)]

        # initial values of v = voltage
        self.v = [-65 for _ in range(self.nr_excit + self.nr_inhibit)]
        # initial values of u = membrane recovery variable
        self.u = np.multiply(self.b, self.v)

        # spike times
        self.firings = []
        self.simulation_time = simulation_time
        self.dt = dt

        # gaussian input
        self.gaussian_input = GaussianInput(
            excit_input=1.5,
            inhibit_input=1.5
        )

        # synaptic constants
        self.synaptic_constraints = SynapticConstants(
            nr_excit=self.nr_excit,
            nr_inhibit=self.nr_inhibit
        )

        # connectivity matrix
        self.connect_matrix = ConnectivityMatrix(
            nr_neurons=self.nr_neurons,
            nr_excit=self.nr_excit,
            nr_inhibit=self.nr_inhibit
        )

    def create_main_input_stimulus(self):
        # sinusoidal spatial modualtion of input strength
        amplitude = 1
        # mean input level to RS cells
        mean_input_lvl_RS = 7

        step = 2*pi / (self.nr_excit - 1)
        stim_input = mean_input_lvl_RS + amplitude * np.sin(
            crange(-pi, pi, step)
        )
        # additional mean input to FS cells
        stim_input += 3.5 * np.ones(self.nr_inhibit)

        return stim_input






