import numpy as np

class GaussianInput:

    def __init__(self, excit_input, inhibit_input):
        
        # gaussian input to excit. neurons
        self.excit = excit_input

        # gaussian input to inhibit. neurons
        self.inhibit = inhibit_input