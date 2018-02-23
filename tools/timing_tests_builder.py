input_file_array = ["test_data/rand1000.fasta",
					"test_data/rand10000.fasta",
					"test_data/rand100000.fasta",
					"test_data/rand1000000.fasta",
					"test_data/randIUPAC1000.fasta",
					"test_data/randIUPAC10000.fasta",
					"test_data/randIUPAC100000.fasta",
					"test_data/randIUPAC1000000.fasta",
					"test_data/randrandIUPAC_N1000.fasta",
					"test_data/randrandIUPAC_N10000.fasta",
					"test_data/randrandIUPAC_N100000.fasta",
					"test_data/randrandIUPAC_N1000000.fasta"]
seq_name_array = ["seq0"]
min_len_array = [10]
max_len_array = [100]
max_gap_array = range(0, 1001, 100)
mismatches_array = range(0, 10)

output_file = open("timing_tests.cfg", "w")
output = ""

for input_file in input_file_array:
	for seq_name in seq_name_array:
		for min_len in min_len_array:
			for max_len in max_len_array:
				for max_gap in max_gap_array:
					for mismatches in mismatches_array:
						output += "file {}\ns {}\nm {}\nM {}\ng {}\nx {}\n\n".format(input_file, seq_name, min_len, max_len, max_gap, mismatches)

output += "END\n"
output_file.write(output)
output_file.close()