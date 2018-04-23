# bch
**bch** is a simple implementation of binary BCH error-correcting code, written in Python 3.6.
Polynomial operations are implemented using [SymPy](http://www.sympy.org) library.
It was made as a homework project for "Error-Correcting Codes and Cryptography" workshops on [Faculty of Mathematics and Information Science of Warsaw University of Technology](http://www.mini.pw.edu.pl).

**This package was made as an excercise. It shouldn't be use anywhere to secure data.**

## How to run?

### Install Python 3.x
You should have Python 3.x installed on your system. To do this on [Fedora OS](https://getfedora.org/ "Get Fedora") you need to execute
```
sudo dnf install python3
```

### Install dependencies
You should have [SymPy](http://www.sympy.org) and [docopt](http://www.docopt.org/) package installed.
```
pip3 install --user sympy
pip3 install --user docopt
```

## How to use it?
To print help screen execute `bch.py` script with `-h` argument.
```
$ ./bch.py -h
BCH v0.1

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
```

As you can see there are few commands:
* `enc` - encode file using BCH code
* `dec` - decode file using BCH code
* `gen` - generate BCH code with given parameters
