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

if (len(sys.argv)) != 7:
	print("Incorrect number of arguments: Must provide 2 filenames and 4 integers (input_text, output_to_be_checked, min_length, max_length, max_gap, mismatches)")
	sys.exit(-1)

filename_input = sys.argv[1]
filename_output = sys.argv[2]
min_length = int(sys.argv[3])
max_length = int(sys.argv[4])
max_gap = int(sys.argv[5])
max_mismatches = int(sys.argv[6])

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
		if (false):
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

print("Min Length: " + str(min_length))
print("Max Length: " + str(max_length))
print("Gap: " + str(max_gap))
print("Mismatches: " + str(max_mismatches))
print("File Input: " + filename_input)
print("File Output: " + filename_output)
print("")
print("Total Palindromes: \t" + str(len(p)))
print("Correct Palindromes: \t" + str(correct))
print("Incorrect Palindromes: \t" + str(incorrect))
