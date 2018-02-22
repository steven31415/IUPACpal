import sys

# Make matching set
match_set = {}
match_set['a'] = ['a']
match_set['c'] = ['c']
match_set['g'] = ['g']
match_set['t'] = ['t']
match_set['u'] = ['t']
match_set['r'] = ['a', 'g']
match_set['y'] = ['c', 't']
match_set['s'] = ['g', 'c']
match_set['w'] = ['a', 't']
match_set['k'] = ['g', 't']
match_set['m'] = ['a', 'c']
match_set['b'] = ['c', 'g', 't']
match_set['d'] = ['a', 'g', 't']
match_set['h'] = ['a', 'c', 't']
match_set['v'] = ['a', 'c', 'g']
match_set['n'] = ['a', 'c', 'g', 't']
match_set['*'] = ['a', 'c', 'g', 't']
match_set['-'] = ['a', 'c', 'g', 't']

# Store alphabet complements
complement = {}
complement['a'] = 't'
complement['c'] = 'g'
complement['g'] = 'c'
complement['t'] = 'a'
complement['u'] = 'a'
complement['r'] = 'y'
complement['y'] = 'r'
complement['s'] = 's'
complement['w'] = 'w'
complement['k'] = 'm'
complement['m'] = 'k'
complement['b'] = 'v'
complement['d'] = 'h'
complement['h'] = 'd'
complement['v'] = 'b'
complement['n'] = 'n'
complement['-'] = 'n'
complement['*'] = 'n'

# Check match of two characters
def checkMatch(x, y):
	x_set = match_set[x]
	y_set = match_set[y]

	for x_set_member in x_set:
		if complement[x_set_member] in y_set:
			return True

	return False

# Extract list of palindromes from a file
def getPalindromes(file):
	palindromes = []

	line = ""

	while ("Palindromes:" not in line):
		line = file.readline()

	while (line != ""):
		line = file.readline() # Get left hand side values
		left = [int(s) for s in line.split() if s.isdigit()]

		line = file.readline() # Skip a line

		line = file.readline() # Get right hand side values
		right = [int(s) for s in line.split() if s.isdigit()]

		line = file.readline() # Skip a line

		if (len(left) == 2 and len(right) == 2):
			palindromes.append(left + right)

	return palindromes

# Check correct number of arguments given
if (len(sys.argv)) != 8:
	print(
"""Incorrect number of arguments. Must provide 7 arguments in order:
1: original_input_file
2: original_input_sequence
3: output_to_be_checked
4: min_length
5: max_length
6: max_gap
7: mismatches"""
	)
	sys.exit(-1)

# Store arguments
filename_input = sys.argv[1]
seq_name = sys.argv[2]
filename_output = sys.argv[3]
min_length = int(sys.argv[4])
max_length = int(sys.argv[5])
max_gap = int(sys.argv[6])
max_mismatches = int(sys.argv[7])

# Extract and store list of palindromes
file_output = open(filename_output, 'r')
p = getPalindromes(file_output)
file_output.close()

# Extract contents of file and find sequence to be extracted
with open (filename_input, "r") as file_input:
	content = file_input.readlines()

data = ""

i = 0;
for line in content:
	i+=1

	if line[0] == '>':
		name = line[1:].split()[0]

		if name == seq_name:
			break;

if (i == len(content)):
	print("Error: Sequence '" + seq_name +  "' not found in file.")
	exit(-1)

# Extract sequence data
for line in content[i:]:
	line = line.replace('\n', '')

	if line == "" or line[0] == " " or line[0] == ";" or line[0] == ">":
		break

	data += line

# Make sequence data lower case
data = data.lower()

# Check correctness of palindromes
correct = 0
incorrect = 0

for palindrome in p:
	outer_left = palindrome[0] - 1
	inner_left = palindrome[1] - 1
	outer_right = palindrome[2] - 1
	inner_right = palindrome[3] - 1

	valid_find = True
	mismatches = 0
	gap = inner_right - inner_left - 1
	length = outer_right - inner_right + 1
	error = ""

	for i in range(0, outer_right - inner_right + 1):
		if not checkMatch(data[outer_left + i], data[outer_right - i]):
			mismatches += 1

			if mismatches > max_mismatches:
				error = "TOO MANY MISMATCHES"
				break

	if gap > max_gap:
		error = "GAP TOO LARGE"

	if length > max_length:
		error = "TOO LONG"

	if length < min_length:
		error = "TOO SHORT"

	if (error == ""):
		correct += 1
	else:
		print("BAD (" + error + "): " + "[" + str(outer_left + 1) + ", " + str(inner_left + 1) + "]-[" + str(inner_right + 1) + ", " + str(outer_right + 1) + "]")

		bad_palindrome = ""
		for i in range(outer_left, inner_left + 1):
			bad_palindrome += data[i]

		bad_palindrome += '\n'

		for i in range(0, inner_left - outer_left + 1):
			if ( checkMatch(data[outer_left + i], data[outer_right - i]) ):
				bad_palindrome += "|"
			else:
				bad_palindrome += " "

		bad_palindrome += "\n"

		for i in range(outer_right, inner_right - 1, -1):
			bad_palindrome += data[i]

		print(bad_palindrome)
		print("\n")

		incorrect += 1

# Output results of correctness checks
print("file_input: " + filename_input)
print("seq_name: " + seq_name)
print("file_output: " + filename_output)
print("min_len: " + str(min_length))
print("max_len: " + str(max_length))
print("gap: " + str(max_gap))
print("mismatches: " + str(max_mismatches))
print("total_palindromes: " + str(len(p)))
print("correct_palindromes: " + str(correct))
print("incorrect_palindromes: " + str(incorrect))
