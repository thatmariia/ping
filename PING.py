from GaussianInput import *
from SynapticConstants import *
from ConnectivityMatrixRing import *
from ConnectivityMatrixGrid import *
from PINGTypes import *
from misc import *

from tqdm import tqdm
import numpy as np
from math import pi

class PING:

    def __init__(self, simulation_time, dt, ping_type, nr_excit, nr_inhibit, nr_oscillators=0):
        self.nr_excit = nr_excit
        self.nr_inhibit = nr_inhibit
        self.nr_neurons = self.nr_excit + self.nr_inhibit

        # timescale of recover variable u
        self.a = np.array([0.02 for _ in range(self.nr_excit)] + [0.1 for _ in range(self.nr_inhibit)])
        # sensitivity of  u to sub-threshold oscillations
        self.b = np.array([0.2 for _ in range(self.nr_excit)] + [0.2 for _ in range(self.nr_inhibit)])
        # membrane voltage after spike (reset)
        self.c = np.array([-65 for _ in range(self.nr_excit)] + [-65 for _ in range(self.nr_inhibit)])
        # spike reset of recover variable u
        self.d = np.array([8 for _ in range(self.nr_excit)] + [2 for _ in range(self.nr_inhibit)])

        # initial values of v = voltage
        self.v = np.array([-65 for _ in range(self.nr_excit + self.nr_inhibit)])
        # initial values of u = membrane recovery variable
        self.u = np.multiply(self.b, self.v)

        # spike timings
        self.firings = []

        # spike times
        self.simulation_time = simulation_time
        self.dt = dt

        # gaussian input
        self.gaussian_input = GaussianInput(
            excit_input=1.5,
            inhibit_input=1.5
        )

        self.stim_input = self._create_main_input_stimulus()

        # synaptic constants
        self.synaptic_constants = SynapticConstants(
            nr_excit=self.nr_excit,
            nr_inhibit=self.nr_inhibit
        )

        # connectivity matrix
        if ping_type == PINGTypes.RING:
            self.connect_matrix = ConnectivityMatrixRing(
                nr_neurons=self.nr_neurons,
                nr_excit=self.nr_excit,
                nr_inhibit=self.nr_inhibit
            )
        if ping_type == PINGTypes.GRID:
            self.connect_matrix = ConnectivityMatrixGrid(
                nr_excit=self.nr_excit,
                nr_inhibit=self.nr_inhibit,
                nr_oscillators=nr_oscillators
            )

    def run(self):
        print("Simulation started")

        for t in tqdm(range(self.simulation_time)):

            # thalamic input
            I = self._get_new_thalamic_input()

            # indices of spikes
            fired = np.argwhere(self.v > 30).flatten()
            self.firings = [self.firings, [t for _ in range(len(fired))] + fired]
            for f in fired:
                self.v[f] = self.c[f]
                self.u[f] += self.d[f]

            # synaptic potentials
            self.synaptic_constants.gampa = self._get_new_gampa()
            self.synaptic_constants.gaba = self._get_new_gaba()
            gsyn = np.append(self.synaptic_constants.gampa, self.synaptic_constants.gaba)

            # defining input to eah neuron as the summation of all synaptic input
            # form all connected neurons
            I = np.add(I, np.matmul(self.connect_matrix.S, gsyn))

            self.v = np.add(self.v, self._addon_to_v(I=I))
            self.v = np.add(self.v, self._addon_to_v(I=I))
            self.u = np.add(self.u, self._addon_to_u())

        print("Simulation ended")

    def _addon_to_u(self):
        return np.multiply(
            self.a,
            np.multiply(self.b, self.v) - self.u
        )

    def _addon_to_v(self, I):
        return 0.5 * (0.04 * self.v**2 + 5 * self.v + 140 - self.u + I)

    def _get_new_gaba(self):
        alpha = self.v[self.nr_excit:] / 10.0 + 2
        z = np.tanh(alpha)
        comp1 = (z + 1) / 2.0
        comp2 = (1 - self.synaptic_constants.gaba) / self.synaptic_constants.rise_gaba
        comp3 = self.synaptic_constants.gaba / self.synaptic_constants.decay_gaba
        new_comp = self.dt * 0.3 * (np.multiply(comp1, comp2) - comp3)
        return np.add(
            self.synaptic_constants.gaba,
            new_comp
        )

    def _get_new_gampa(self):
        alpha = self.v[:self.nr_excit] / 10.0 + 2
        z = np.tanh(alpha)
        comp1 = (z + 1) / 2.0
        comp2 = (1 - self.synaptic_constants.gampa) / self.synaptic_constants.rise_ampa
        comp3 = self.synaptic_constants.gampa / self.synaptic_constants.decay_ampa
        new_comp = self.dt * 0.3 * (np.multiply(comp1, comp2) - comp3)
        return np.add(
            self.synaptic_constants.gampa,
            new_comp
        )

    def _get_new_thalamic_input(self):
        return np.add(
            self.stim_input,
            np.append(
                self.gaussian_input.excit * np.random.randn(self.nr_excit),
                self.gaussian_input.inhibit * np.random.randn(self.nr_inhibit)
            )
        )


    def _create_main_input_stimulus(self):
        # sinusoidal spatial modualtion of input strength
        amplitude = 1
        # mean input level to RS cells
        mean_input_lvl_RS = 7

        step = 2*pi / (self.nr_excit - 1)
        stim_input = mean_input_lvl_RS + amplitude * np.sin(
            crange(-pi, pi, step)
        )
        # additional mean input to FS cells
        stim_input = np.append(stim_input, 3.5 * np.ones(self.nr_inhibit))

        return stim_input






