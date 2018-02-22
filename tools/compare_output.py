import sys
import os.path
from operator import itemgetter

def comparePals(a, b):
	i = 0

	while (a[i] == b[i]):
		i += 1

		if i == 4:
			return "="

	if a[i] < b[i]:
		return "<"
	else:
		return ">"

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

if (len(sys.argv)) != 3:
	print("Incorrect number of arguments: Must provide 2 filenames.")
	sys.exit(-1)

filename1 = sys.argv[1]
filename2 = sys.argv[2]

if not os.path.exists(filename1):
	print("File '" + filename1 + "' was not found.")
	sys.exit(-1)

if not os.path.exists(filename2):
	print("File '" + filename2 + "' was not found.")
	sys.exit(-1)

file1 = open(filename1, 'r')
p1 = getPalindromes(file1)
file1.close()

file2 = open(filename2, 'r')
p2 = getPalindromes(file2)
file2.close()

p1_sorted = sorted(p1, key=itemgetter(0, 1, 2, 3))
p2_sorted = sorted(p2, key=itemgetter(0, 1, 2, 3))

common = 0
one_only = 0
two_only = 0

i = 0;
j = 0;

while (i < len(p1_sorted) or j < len(p2_sorted)):

	if (i >= len(p1_sorted)):
		two_only += 1
		j += 1
	elif (j >= len(p2_sorted)):
		one_only += 1
		i += 1
	else:
		p1_pal = p1_sorted[i]
		p2_pal = p2_sorted[j]

		if comparePals(p1_pal, p2_pal) == "=":
			common += 1
			i += 1
			j += 1
		elif comparePals(p1_pal, p2_pal) == "<":
			one_only += 1
			i += 1
		else:
			two_only += 1
			j += 1

print("file_1: " + filename1)
print("file_2: " + filename2)
print("file_1_total: " + str(len(p1)))
print("file_2_total: " + str(len(p2)))
print("in_both: " + str(common))
print("only_in_file_1: " + str(one_only))
print("only_in_file_2: " + str(two_only))
