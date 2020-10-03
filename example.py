from IUPACpal import find_inverted_repeats, config
import re

# PRINT CONFIG HELP
# config()

# EXAMPLE RUN
inverted_repeats = find_inverted_repeats(
	input_file = 'test_data/rand1000.fasta',
	seq_name = 'seq0',
	min_len = 10,
	max_len = 20,
	max_gap = 100,
	mismatches = 2,
	output_file = 'IUPACpal.out',
	)

# PRINT FOUND INVERTED REPEATS / ERROR MESSAGE
if isinstance(inverted_repeats , str):
	print inverted_repeats
else:
	print("FORMAT: (left_strand_start, left_strand_end), (right_strand_start, right_strand_end)")
	print("FOUND INVERTED REPEATS:")
	for ir in inverted_repeats:
		print(ir)