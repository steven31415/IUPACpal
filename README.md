# IUPACpal

An inverted repeat is a DNA sequence followed downstream by its reverse complement, potentially with a gap in the centre. 

IUPACpal is an exact tool for efficient identification of inverted repeats in IUPAC-encoded DNA sequences as substrings of a large text, allowing also for potential mismatches and gaps.

## Project Structure

| FILE / FOLDER | PURPOSE |
| :--- | :--- |
| `Makefiles` | Folder of Makefiles enabling alternative compilation parameters
| `sdsl-lite` | Succinct Data Structure Library (unzipped folder)
| `test_data` | Folder containing test data that may be used as input for IUPACpal
| `test_results` | Folder containing results from test runs
| `tools` | Folder containg various Python tools for testing purposes
| `Makefile.gcc` | Makefile for compiling the project
| `README.md` | This file
| `main.cc` | IUPACpal main code
| `main.h` | IUPACpal header code
| `pre-install.sh` | Pre-installation script
| `sdsl-lite.tar.gz` | Succinct Data Structure Library (zipped folder)
| `timing_tests.sh` | Execute timing tests (requires timing_tests.cfg)

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

To run IUPACpal requires the pre-installation of:
- cmake
- C++ compiler
- libdivsufsort (https://github.com/y-256/libdivsufsort)
- sdsl (https://github.com/simongog/sdsl-lite)

To build testing configuration files requires pre-installation of:
- Python (https://www.python.org/)

The necessary files for libdivsufsort and sdsl are already include in this repository.

### Installation

Pre-installation of the above tools and libraries may be easily completed by running the following commands on a Ubuntu system within the project directory:

```
$ sudo apt-get install cmake
$ sudo apt-get install gcc
$ sudo apt-get install python3
$ sudo ./pre-install.sh
```

The IUPACpal software may then be compiled:

```
$ make -f Makefile.gcc
```

This will compile the program for 64-bit integers using the recommended type of Range Minimum Queries (**rmq1**). This requires double the amount of memory as 32-bit integers.

There are alternative makefiles available in the `Makefiles` folder, which provide an alternative implementation of Range Minimum Queries (**rmq2**) and a choice between 32-bit and 64-bit integers .

## Running IUPACpal

After compilation the binary file `IUPACpal` will be created in the working
directory, you may run it within this directory via the command:

```
$ ./IUPACpal <arguments>
```

To see the full list of required inputs, run the command with no arguments:

```
$ ./IUPACpal
```

### Input Parameters

| FLAG | PARAMETER | TYPE | DEFAULT | DESCRIPTION |
| :--- | :--- | :--- | :--- | :--- |
| -f | input_file | string | input.fasta | Input filename (FASTA). |
| -s | seq_name | string | seq0 | Input sequence name. |
| -m | min_len | integer | 10 | Minimum length. |
| -M | max_len | integer | 100 | Maximum length. |
| -g | max_gap | integer | 100 | Maximum permissible gap. |
| -x | mismatches | integer | 0 | Maximum permissible mismatches. |
| -o | output_file| string | IUPACpal.out | Output filename. |

### Examples

#### Example 1

Run on the input file `test_data/d00596.fasta` with default parameters:
```
$ ./IUPACpal -f test_data/d00596.fasta
```

#### Example 2

Run on the input file `test_data/test2.fasta` with parameters:
- sequence name in fasta file: seq2
- Maximum permissible gap: 3
- Minimum length: 5
- all other parameters default
```
$ ./IUPACpal -f test_data/test2.fasta -s seq2 -g 3 -m 5
```

#### Example 3

Run on the input file `test_data/rand1000000.fasta` with parameters:
- Minimum length: 20
- Maximum length: 25
- Maximum permissible gap: 8
- Maximum permissible mismatches: 5
- Output filename: output.txt
```
$ ./IUPACpal -f test_data/rand1000000.fasta -s seq0 -m 20 -M 25 -g 8 -x 5 -o output.txt
```

## Testing

Automated testing can be performed through the use of `timing_tests.sh` and a configuration file `timing_tests.cfg`.

To run the default series of test on IUPACpal, first create the `timing_tests.cfg` file by running:

```
$ python3 tools/timing_tests_builder.py
```

This creates the default configruation file `timing_tests.cfg`. Testing of IUPACpal can then be performed by running:

```
$ ./timing_tests.sh
```

Ensure `timing_tests.sh` and `timing_tests.cfg` are in the same directory. To modify the test parameters, you can either directly modify `timing_tests.cfg` or alternatively modify the Python script `tools/timing_tests_builder.py` to enable the script to build a different set of tests over the desired range of input parameters.

Test results are written to an output file `test_results_<CURRENT_DATE_AND_TIME>.csv`.

## Author

**Steven Watts** - [steven31415](https://github.com/steven31415)

## License

This project is licensed under the MIT License.
