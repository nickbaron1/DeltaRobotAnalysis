from inverse_kinematics import InverseKinematics

def main():
    ik = InverseKinematics()
    angles = ik.compute_ik([0, 0, -50])
    print(angles)


if __name__ == "__main__":
    main()