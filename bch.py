#!/usr/bin/env python3
"""BCH v0.1

Usage:
  bch.py [options] enc CODE_FILE [FILE]
  bch.py [options] dec CODE_FILE [FILE]
  bch.py [options] gen N B D CODE_FILE
  bch.py (-h | --help)
  bch.py --version

Options:
  -b, --block        Interpret input/output as
                       block stream.
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
from sympy import Poly
from bch.bchcodegenerator import BchCodeGenerator
from bch.bchcoder import BchCoder
from padding.padding import *
import numpy as np
import sys
import logging

log = logging.getLogger("bch")


def generate(n, b, d, code_file):
    bch_gen = BchCodeGenerator(n, b, d)
    r, g = bch_gen.gen()
    np.savez_compressed(code_file, n=n, b=b, d=d, r=r.all_coeffs()[::-1], g=g.all_coeffs()[::-1])
    log.info("BCH code saved to {} file".format(code_file))


def encode(code_file, input_arr, block=False):
    code = np.load(code_file, allow_pickle=True)
    bch = BchCoder(int(code['n']), int(code['b']), int(code['d']), Poly(code['r'][::-1], alpha),
                   Poly(code['g'][::-1], x))
    if not block:
        if len(input_arr) > bch.k:
            raise Exception("Input is too large for current BCH code (max: {})".format(bch.k))
        return bch.encode(Poly(input_arr[::-1], x))[::-1]

    input_arr = padding_encode(input_arr, bch.k)
    input_arr = input_arr.reshape((-1, bch.k))
    output = np.array([])
    block_count = input_arr.shape[0]
    for i, b in enumerate(input_arr, start=1):
        log.info("Processing block {} out of {}".format(i, block_count))
        next_output = np.array(bch.encode(Poly(b[::-1], x))[::-1])
        if len(next_output) < bch.n:
            next_output = np.pad(next_output, (0, bch.n - len(next_output)), 'constant')
        output = np.concatenate((output, next_output))
    return output


def decode(code_file, input_arr, block=False):
    code = np.load(code_file, allow_pickle=True)
    bch = BchCoder(int(code['n']), int(code['b']), int(code['d']), Poly(code['r'][::-1], alpha),
                   Poly(code['g'][::-1], x))
    if not block:
        if len(input_arr) > bch.n:
            raise Exception("Input is too large for current BCH code (max: {})".format(bch.n))
        return bch.decode(Poly(input_arr[::-1], x))[::-1]

    input_arr = input_arr.reshape((-1, bch.n))
    output = np.array([])
    block_count = input_arr.shape[0]
    for i, b in enumerate(input_arr, start=1):
        log.info("Processing block {} out of {}".format(i, block_count))
        next_output = np.array(bch.decode(Poly(b[::-1], x))[::-1])
        if len(next_output) < bch.k:
            next_output = np.pad(next_output, (0, bch.k - len(next_output)), 'constant')
        output = np.concatenate((output, next_output))

    return padding_decode(output, bch.k)


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
    block = bool(args['--block'])
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
            input_arr = np.array(eval(input))
        else:
            input_arr = np.unpackbits(np.frombuffer(input, dtype=np.uint8))
        input_arr = np.trim_zeros(input_arr, 'b')
        log.info("POLYNOMIAL DEGREE: {}".format(max(0, len(input_arr) - 1)))
        log.debug("BINARY: {}".format(input_arr))

    if args['gen']:
        generate(int(args['N']), int(args['B']), int(args['D']), args['CODE_FILE'])
    elif args['enc']:
        output = encode(args['CODE_FILE'], input_arr, block=block)
    elif args['dec']:
        output = decode(args['CODE_FILE'], input_arr, block=block)

    if not args['gen']:
        if poly_output:
            print(list(output.astype(np.int)))
        else:
            sys.stdout.buffer.write(np.packbits(np.array(output).astype(np.int)).tobytes())
