from RingPING import *

if __name__ == "__main__":

    ringPING = RingPING(
        nr_excit=4,
        nr_inhibit=2,
        simulation_time=8,
        dt=1
    )
    ringPING.run()

