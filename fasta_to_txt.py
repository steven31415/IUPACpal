import sys

output = ""
in_file = open("pure.gb", "r")
out_file = open("pure.parse.gb", "w")

for line in in_file:
	line = line[10:]
	line = line.split(" ")
	output += "".join(line).rstrip()

out_file.write(output)


