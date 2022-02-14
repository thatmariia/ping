from PING import *

if __name__ == "__main__":

    ringPING = PING(
        simulation_time=8,
        dt=1,
        nr_excit=4,
        nr_inhibit=2,
    )
    ringPING.run()

