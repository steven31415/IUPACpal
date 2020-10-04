import re
import subprocess

CONFIG = """
PARAMETER       TYPE      DEFAULT         DESCRIPTION
input_file      <str>     input.fasta     Input filename (FASTA).
seq_name        <str>     seq0            Input sequence name.
min_len         <int>     10              Minimum length.
max_len         <int>     100             Maximum length.
max_gap         <int>     100             Maximum permissible gap.
mismatches      <int>     0               Maximum permissible mismatches.
output_file     <str>     IUPACpal.out    Output filename.
"""

def _run(cmd):
	proc = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
	stdout, stderr = proc.communicate()
	return proc.returncode, stdout, stderr

def _extract_locs(data):
	x = data.split(' ')
	a = int(x[0])
	b = int(x[-1])
	return a, b

def config():
	print(CONFIG)

def find_inverted_repeats(	input_file='input.fasta',
							seq_name='seq0',
							min_len='10',
							max_len='100',
							max_gap='100',
							mismatches='0',
							output_file='IUPACpal.out'):

	code, out, err = _run([
							'./IUPACpal',
							'-f', input_file,
							'-s', seq_name,
							'-m', str(min_len),
							'-M', str(max_len),
							'-g', str(max_gap),
							'-x', str(mismatches),
							'-o', output_file,
							])

	inverted_repeats = []

	valid_run = (not 'Error' in str(out))

	if valid_run:
		with open(output_file) as f_in:
			lines = list(line for line in (l.strip() for l in f_in) if line)

		lines = lines[lines.index('Palindromes:') + 1:]

		chunks = [lines[i : i + 3] for i in range(0, len(lines), 3)]

		for chunk in chunks:
			left_start, left_end = _extract_locs(chunk[0])
			right_end, right_start = _extract_locs(chunk[2])
			inverted_repeats.append(( (left_start, left_end), (right_start, right_end) ))

		return inverted_repeats
	else:
		return str(out)
