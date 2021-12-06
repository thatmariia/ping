from PING import *
from PINGTypes import *


if __name__ == "__main__":

    ringPING = PING(
        nr_excit=4,
        nr_inhibit=2,
        simulation_time=8,
        dt=1,
        ping_type=PINGTypes.RING
    )
    ringPING.run()

