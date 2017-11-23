import sys

if (len(sys.argv)) != 2:
	print("Incorrect number of arguments: Must provide 1 filename")
	sys.exit(-1)

filename = sys.argv[1]

with open(filename, 'r') as file:
    seq = file.readlines()

seq = "".join(seq)

output = ""

output += "LOCUS\n"
output += "ORIGIN"

for i in range(0, len(seq)):
	if ((i % 40) == 0):
		output += "\n" + str(i + 1)

	if ((i % 10) == 0):
		output += " "

	output += seq[i]

output += "\n//\n"

output_file = open(filename + ".gb", 'w')
output_file.write(output)
output_file.close()

file.close()