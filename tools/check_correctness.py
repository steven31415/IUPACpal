import sys

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



def checkMatch(x, y):
	x_set = match_set[x]
	y_set = match_set[y]

	for x_set_member in x_set:
		if complement[x_set_member] in y_set:
			return True

	return False



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


if (len(sys.argv)) != 4:
	print("Incorrect number of arguments: Must provide 2 filesnames and 1 integer (input_text, output_to_be_checked, max_mismatches)")
	sys.exit(-1)

filename_input = sys.argv[1]
filename_output = sys.argv[2]
max_mismatches = int(sys.argv[3])
print(max_mismatches)

file_output = open(filename_output, 'r')
p = getPalindromes(file_output)
file_output.close()

with open (filename_input, "r") as file_input:
    data = file_input.read().replace('\n', '')

correct = 0
incorrect = 0

for palindrome in p:
	outer_left = palindrome[0] - 1
	inner_left = palindrome[1] - 1
	outer_right = palindrome[2] - 1
	inner_right = palindrome[3] - 1

	valid_find = True
	mismatches = 0

	for i in range(0, outer_right - inner_right + 1):
		if not checkMatch(data[outer_left + i], data[outer_right - i]):
			mismatches += 1
			
			if mismatches > max_mismatches:
				valid_find = False
				break

	if (valid_find):
		correct += 1
	else:
		incorrect += 1


print("File Input: " + filename_input)
print("File Output: " + filename_output)
print("")
print("Total palindromes: \t" + str(len(p)))
print("Correct palindromes: \t" + str(correct))
print("Incorrect palindromes: \t" + str(incorrect))
