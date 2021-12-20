from Oscillator import *
from NeuronTypes import *

import numpy as np
from misc import *


class ConnectivityMatrixGrid():

    def __init__(self, nr_excit, nr_inhibit, nr_oscillators):

        assert nr_excit % nr_oscillators == 0, "Cannot allocated equal number of excitatory neurons to each oscillator. Make sure the number of oscillators divides the number of excitatory neurons."
        assert nr_inhibit % nr_oscillators == 0, "Cannot allocated equal number of inhibitory neurons to each oscillator. Make sure the number of oscillators divides the number of inhibitory neurons."
        assert int(math.sqrt(nr_oscillators)) == math.sqrt(nr_oscillators), "The oscillators should be arranged in a square grid. Make sure the number of oscillators is a perfect square."

        self.nr_oscillators = nr_oscillators
        self.oscillators = []

        # number of neurons of each type in each oscillator
        self.nr_excit_per_oscillator = nr_excit // nr_oscillators
        self.nr_inhibit_per_oscillator = nr_inhibit // nr_oscillators

        # maps the neuron ID to an oscillator it belongs to
        self.neuron_oscillator_map = {
            NeuronTypes.E: {},
            NeuronTypes.I: {}
        }

        self.grid_size = int(math.sqrt(nr_oscillators)) # now assuming the grid is square

        self.assign_oscillators()

        # excitatory  to excitatory
        self.EE = 0.004
        self.sEE = 0.4

        # excitatory  to inhibitory
        self.EI = 0.07
        self.sEI = 0.3

        # inhibitory to excitatory
        self.IE = -0.04
        self.sIE = 0.3

        # inhibitory to inhibitory
        self.II = -0.015
        self.sII = 0.3

        self.KEE, self.KII, self.KEI, self.KIE = self.get_KXXs(nr_excit=nr_excit, nr_inhibit=nr_inhibit)

        # TODO:: compute self.S
        # self.S = self.get_S(
        #     nr_neurons=nr_neurons,
        #     nr_excit=nr_excit,
        #     nr_inhibit=nr_inhibit
        # )

    def get_KXXs(self, nr_excit, nr_inhibit):
        dist_EE = self._get_neurons_dist(
            X1=NeuronTypes.E,
            X2=NeuronTypes.E,
            nr1=nr_excit,
            nr2=nr_excit
        )
        dist_II = self._get_neurons_dist(
            X1=NeuronTypes.I,
            X2=NeuronTypes.I,
            nr1=nr_inhibit,
            nr2=nr_inhibit
        )
        dist_EI = self._get_neurons_dist(
            X1=NeuronTypes.E,
            X2=NeuronTypes.I,
            nr1=nr_excit,
            nr2=nr_inhibit
        )
        dist_IE = self._get_neurons_dist(
            X1=NeuronTypes.I,
            X2=NeuronTypes.E,
            nr1=nr_inhibit,
            nr2=nr_excit
        )

        return self._compute_KXX(dist=dist_EE, XX=self.EE, sXX=self.sEE), \
               self._compute_KXX(dist=dist_II, XX=self.II, sXX=self.sII), \
               self._compute_KXX(dist=dist_EI, XX=self.EI, sXX=self.sEI), \
               self._compute_KXX(dist=dist_IE, XX=self.IE, sXX=self.sIE)

    def _compute_KXX(self, dist, XX, sXX):
        KXX = XX * np.exp(np.true_divide(-dist, sXX))
        return KXX

    def _get_neurons_dist(self, X1, X2, nr1, nr2):
        """
        Computes the matrix of Euclidian distances between each pair of neurons of given types.
        :param X1: neurons type 1
        :param X2: neurons type 2
        :param nr1: number of neurons of type 1
        :param nr2: number of neurons of type 2
        :return: dist: R^(nr1, nr2)
        """
        dist = np.zeros((nr1, nr2))

        print(f"\nMatrix of distances between {X1.value} and {X2.value} neurons:")

        for id1 in range(nr1):
            print()
            for id2 in range(nr2):
                # finding to which oscillators the neurons belong
                oscillator1 = self.oscillators[self.neuron_oscillator_map[X1][id1]]
                oscillator2 = self.oscillators[self.neuron_oscillator_map[X2][id2]]

                # computing the distance between the found oscillators
                # (which = the distance between neurons in those oscillators)
                # assuming unit distance for now
                dist[id1][id2] = euclidian_dist_R2(
                    x=oscillator1.location[0] - oscillator2.location[0],
                    y=oscillator1.location[1] - oscillator2.location[1]
                )

                print(round(dist[id1, id2], 1), end=" ")

        print()

        return dist

    def assign_oscillators(self):
        """
        Creates oscillators and assigns grid locations and the same number of each type of neurons to them.
        See Oscillator.py
        """
        for i in range(self.nr_oscillators):
            x = i // self.grid_size
            y = i % self.grid_size

            excit_ids = []

            for id in range(i * self.nr_excit_per_oscillator, (i + 1) * self.nr_excit_per_oscillator + 1):
                excit_ids.append(id)
                self.neuron_oscillator_map[NeuronTypes.E][id] = i

            inhibit_ids = []

            for id in range(i * self.nr_inhibit_per_oscillator, (i + 1) * self.nr_inhibit_per_oscillator + 1):
                inhibit_ids.append(id)
                self.neuron_oscillator_map[NeuronTypes.I][id] = i

            oscillator = Oscillator(
                location=(x, y),
                excit_ids=excit_ids,
                inhibit_ids=inhibit_ids
            )
            self.oscillators.append(oscillator)






