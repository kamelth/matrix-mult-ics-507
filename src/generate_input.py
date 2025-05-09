#!/usr/bin/env python3
"""
generate_input.py

Generates a matrix-mul input file with
  • first line: n (power of 2)
  • next n lines: n random ints in [-9,9] (matrix A)
  • next n lines: n random ints in [-9,9] (matrix B)

Usage:
    python3 generate_input.py <n> <output_file>
"""
import sys
import random


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <n> <output_file>")
        sys.exit(1)

    n = int(sys.argv[1])
    # check power of two
    if n < 1 or (n & (n - 1)) != 0:
        print("Error: n must be a power of 2 and ≥1")
        sys.exit(1)

    out_fname = sys.argv[2]
    with open(out_fname, 'w') as f:
        # write dimension
        f.write(f"{n}\n")
        # write A
        for _ in range(n):
            row = [str(random.randint(-9, 9)) for __ in range(n)]
            f.write(" ".join(row) + "\n")
        # write B
        for _ in range(n):
            row = [str(random.randint(-9, 9)) for __ in range(n)]
            f.write(" ".join(row) + "\n")

    print(f"Generated `{out_fname}` with n={n}")

if __name__ == "__main__":
    main()

