import numpy as np
from math import pi, exp
import cmath
from misc import *

class ConnectivityMatrixRing():

    def __init__(self, nr_neurons, nr_excit, nr_inhibit):
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

        self.S = self.get_S(
            nr_neurons=nr_neurons,
            nr_excit=nr_excit,
            nr_inhibit=nr_inhibit
        )

    def get_S(self, nr_neurons, nr_excit, nr_inhibit):
        S = np.zeros((nr_neurons, nr_neurons))
        S[:nr_excit, :nr_excit] = self.KEE
        S[nr_excit:nr_neurons, nr_excit:nr_neurons] = self.KII
        S[:nr_excit, nr_excit:nr_neurons] = self.KIE.T
        S[nr_excit:nr_neurons, :nr_excit] = self.KEI.T
        S = np.nan_to_num(S)
        return S

    def get_KXXs(self, nr_excit, nr_inhibit):
        dist_EE = self._get_dist(
            nr1=nr_excit,
            nr2=nr_excit,
            op1=self._get_op(nr=nr_excit),
            op2=self._get_op(nr=nr_excit)
        )
        dist_II = self._get_dist(
            nr1=nr_inhibit,
            nr2=nr_inhibit,
            op1=self._get_op(nr=nr_inhibit),
            op2=self._get_op(nr=nr_inhibit)
        )
        dist_EI = self._get_dist(
            nr1=nr_excit,
            nr2=nr_inhibit,
            op1=self._get_op(nr=nr_excit),
            op2=self._get_op(nr=nr_inhibit)
        )
        dist_IE = self._get_dist(
            nr1=nr_inhibit,
            nr2=nr_excit,
            op1=self._get_op(nr=nr_inhibit),
            op2=self._get_op(nr=nr_excit)
        )

        return self._compute_KXX(dist=dist_EE, XX=self.EE, sXX=self.sEE), \
               self._compute_KXX(dist=dist_II, XX=self.II, sXX=self.sII), \
               self._compute_KXX(dist=dist_EI, XX=self.EI, sXX=self.sEI), \
               self._compute_KXX(dist=dist_IE, XX=self.IE, sXX=self.sIE)

    def _compute_KXX(self, dist, XX, sXX):
        KXX = XX * np.exp(np.true_divide(-dist, sXX))
        return KXX

    def _get_op(self, nr):
        step = 2 * pi / (nr - 1)
        op = crange(-pi, pi, step)
        return op

    def _get_dist(self, nr1, nr2, op1, op2):
        dist = np.zeros((nr1, nr2))

        for i in range(nr1):
            zs_nom = [cmath.exp(complex(real=0.0, imag=op1[i]))] * len(op2)
            zs_denom = [cmath.exp(complex(real=0.0, imag=j)) for j in op2]
            zs = np.true_divide(zs_nom, zs_denom)
            angles = abs(np.angle(z=zs))
            dist[i] = angles

        dist[dist < 0.001] = None
        return dist




