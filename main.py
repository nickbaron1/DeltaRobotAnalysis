from delta_robot import DeltaRobot
import numpy as np


def main():
    dr = DeltaRobot()
    x, y, z = 0, 0, -900
    angles = dr.compute_ik(x, y, z)
    print(angles)
    x, y, z = dr.forward_ik(angles)
    print(x, y, z)
    J = dr.compute_jacobian(angles)
    print(J)


if __name__ == "__main__":
    main()