from NeuronTypes import *

class Oscillator:

    def __init__(self, location, excit_ids, inhibit_ids):
        self.location = location
        self.ids = {
            NeuronTypes.E: excit_ids,
            NeuronTypes.I: inhibit_ids
        }