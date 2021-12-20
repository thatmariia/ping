from PING import *
from PINGTypes import *


if __name__ == "__main__":

    # ringPING = PING(
    #     simulation_time=8,
    #     dt=1,
    #     ping_type=PINGTypes.RING,
    #     nr_excit=4,
    #     nr_inhibit=2,
    # )
    # ringPING.run()

    gridPING = PING(
        simulation_time=8,
        dt=1,
        ping_type=PINGTypes.GRID,
        nr_excit=8,
        nr_inhibit=4,
        nr_oscillators=4
    )
    # gridPING.run()

