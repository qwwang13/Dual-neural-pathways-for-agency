"""Generate the paired participant draws used for the ITPC bootstrap."""

import argparse
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True)
    parser.add_argument("--participants", type=int, required=True)
    parser.add_argument("--replicates", type=int, required=True)
    parser.add_argument("--seed", type=int, required=True)
    args = parser.parse_args()
    rng = np.random.Generator(np.random.PCG64(args.seed))
    draws = rng.integers(0, args.participants,
                         size=(args.replicates, args.participants), dtype=np.int64)
    # R uses one-based indices. The draw order matches the original analysis.
    (draws + 1).astype("<i4").tofile(args.output)


if __name__ == "__main__":
    main()
