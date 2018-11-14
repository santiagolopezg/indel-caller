# mk_cigars
import glob, os

for file in glob.glob('aligned_S*'):
	if '.txt' not in file:
		os.rename(file, file+'.txt')
		
	with open(file, 'r') as f:
		data = [i for i in f.read().split('\r') if len(i)>1]
				
	fname = file.split('_')[1]

with open('refseq_guide_cigars_{0}'.format(fname), 'w') as g:
	for read in data:
		read = read.split('\t')
		gene, cigar = read[2], read[5]
		
		print gene, cigar
		
		g.write(gene+'-'+cigar+'\n')
		