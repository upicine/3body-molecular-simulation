import random
import argparse

def generate_test(filename, lines_n):
    with open(filename, 'w+') as f:
        for l in range(lines_n):
            x = random.uniform(0.0, 2.0)
            y = random.uniform(0.0, 2.0)
            z = random.uniform(0.0, 2.0)
            vx = random.uniform(0.0, 2.0)
            vy = random.uniform(0.0, 2.0)
            vz = random.uniform(0.0, 2.0)
            print(x, y, z, vx, vy, vz, file=f)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('file', help='File with molcules')
    parser.add_argument('n', type=int, help='Number of lines')
    args = parser.parse_args()

    generate_test(args.file, args.n)


if __name__ == '__main__':
    main()