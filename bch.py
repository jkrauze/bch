#!/usr/bin/env python3
"""BCH v0.1

Usage:
  bch.py [options] enc CODE_FILE [FILE]
  bch.py [options] dec CODE_FILE [FILE]
  bch.py [options] gen N B D CODE_FILE
  bch.py (-h | --help)
  bch.py --version

Options:
  -i, --poly-input   Interpret input as polynomial
                       represented by integer array.
  -o, --poly-output  Interpret output as polynomial
                       represented by integer array.
  -h, --help         Show this screen.
  --version          Show version.
  -d, --debug        Debug mode.
  -v, --verbose      Verbose mode.

"""
from docopt import docopt
from sympy.abc import x, alpha
from sympy import ZZ, Poly
from bch.bchcodegenerator import BchCodeGenerator
from bch.bchcoder import BchCoder
import numpy as np
import sys
import logging

log = logging.getLogger("bch")


def generate(n, b, d, code_file):
    bch_gen = BchCodeGenerator(n, b, d)
    r, g = bch_gen.gen()
    np.savez_compressed(code_file, n=n, b=b, d=d, r=r.all_coeffs()[::-1], g=g.all_coeffs()[::-1])
    log.info("BCH code saved to {} file".format(code_file))


def encode(code_file, input_arr):
    code = np.load(code_file)
    bch = BchCoder(int(code['n']), int(code['b']), int(code['d']), Poly(code['r'][::-1], alpha),
                   Poly(code['g'][::-1], x))

    if len(input_arr) > bch.k:
        raise Exception("Input is too large for current BCH code (max: {})".format(bch.k))
    return bch.encode(Poly(np.pad(input_arr, (0, bch.k - len(input_arr)), 'constant')[::-1], x)).all_coeffs()[::-1]


def decode(code_file, input_arr):
    code = np.load(code_file)
    bch = BchCoder(int(code['n']), int(code['b']), int(code['d']), Poly(code['r'][::-1], alpha),
                   Poly(code['g'][::-1], x))
    if (len(input_arr) + 1 == bch.n):
        log.info("Trimming input to {}".format(bch.n))
        input_arr = input_arr[:-1]
    elif len(input_arr) > bch.n:
        raise Exception("Input is too large for current BCH code (max: {})".format(bch.n))
    return bch.decode(Poly(np.pad(input_arr, (0, bch.n - len(input_arr)), 'constant')[::-1], x))[::-1]


if __name__ == '__main__':
    args = docopt(__doc__, version='BCH v0.1')
    root = logging.getLogger()
    root.setLevel(logging.DEBUG)
    ch = logging.StreamHandler(sys.stdout)
    if args['--debug']:
        ch.setLevel(logging.DEBUG)
    elif args['--verbose']:
        ch.setLevel(logging.INFO)
    else:
        ch.setLevel(logging.WARN)
    root.addHandler(ch)

    log.debug(args)
    poly_input = bool(args['--poly-input'])
    poly_output = bool(args['--poly-output'])
    input_arr, output = None, None
    if not args['gen']:
        if args['FILE'] is None or args['FILE'] == '-':
            input = sys.stdin.read() if poly_input else sys.stdin.buffer.read()
        else:
            with open(args['FILE'], 'rb') as file:
                input = file.read()
        log.info("---INPUT---")
        log.info(input)
        log.info("-----------")
        if poly_input:
            input_arr = eval(input)
        else:
            input_arr = np.unpackbits(np.frombuffer(input, dtype=np.uint8))
        log.info("POLYNOMIAL LENGTH: {}".format(len(input_arr)))
        log.debug("BINARY: {}".format(input_arr))

    if args['gen']:
        generate(int(args['N']), int(args['B']), int(args['D']), args['CODE_FILE'])
    elif args['enc']:
        output = encode(args['CODE_FILE'], input_arr)
    elif args['dec']:
        output = decode(args['CODE_FILE'], input_arr)

    if not args['gen']:
        if poly_output:
            print(list(output))
        else:
            sys.stdout.buffer.write(np.packbits(np.array(output).astype(np.int)).tobytes())
