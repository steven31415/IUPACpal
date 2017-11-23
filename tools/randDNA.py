import random

alphabet = "acgt"
iupac = "uryswkmbdhv"
n = "n"

length = input('Enter sequence length: ')
include_iupac = input('Include IUPAC characters (Y/N): ')
include_n = input('Include N character (Y/N): ')

if include_iupac.lower() != 'n':
	alphabet += iupac

if include_n.lower() != 'n':
	alphabet += n

print("Used alphabet:", alphabet)

output = ""
for i in range(0, int(length)):
	output += random.choice(alphabet)

F = open('randDNA.out', 'w')
F.write(output)
F.close()

