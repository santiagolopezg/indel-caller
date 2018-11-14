#filter_tallies.py

'''
input is tallies files and alignment file (from samtools)
'''

import glob, os
import math
from operator import itemgetter
from itertools import groupby

sites_of_interest = {}
with open('refseq_guide.txt', 'r') as f:
	data = f.read().split('\r')
	
	for d in data:
		d = d.split('\t')
		ref, guide, seq, pam = d[0], d[1], d[2], d[3]
		if 'CC' in pam:
			pos = seq.find(guide) + 3
			
		else:
			pos = seq.find(guide) + 20 - 3 

		sites_of_interest[ref] = pos		

name_refseq = {}

with open('refseq_guide.txt', 'r') as f:
	sqs = f.read().split('\r')
	for sq in sqs:
		sq = sq.split('\t')
		name_refseq[sq[0]]=sq[2]

for file in glob.glob('aligned_S*.txt'):

	reads_dict = {}
	ref_tally = {}
	reinge = []	
	nm = file.split('_')[1].split('.')[0]
	
	with open(file, 'r') as w:
		reads = w.read().split('\n')
		reads = [i.strip('\n') for i in reads[:len(reads)-1]]
		for read in reads:
			read = read.split('\t')
			id = read[2]+'-'+read[5]
			reads_dict[id] = read[9]

	for acname in name_refseq.keys():
		with open('tally_refseq_guide_cigars_{0}_{1}.txt'.format(nm, acname), 'r') as f:
			tallies = f.read().split('\n')
			tallies = tallies[:len(tallies)-1]
			tally_nucl = {}

			for tally in tallies:
				tal = tally.split('\t')

				if ',' in tal[1]:
					v = [] #this is the split list
					y = [] #this is a list to store corresponding nucleotides
					pos, f = tal[1].split(','), tal[2]
					pos = [int(i) for i in pos[:len(pos)-1]]
					for k, g in groupby(enumerate(pos), lambda (i,x):i-x):
						v.append(map(itemgetter(1), g))
				
					c = ''
					for i in xrange(len(v)):
						for j in v[i]:
							c+= name_refseq[acname][j-1]

						c+=','
					tally_nucl[tally] = c	

				elif '-' in tal[1]:						
				
					id = tal[0]	
					tag = id.split('-')[1]
				
					pos_l, pos_u = tal[1].split('-')[0], tal[1].split('-')[1].split(' ')[0]
					d = tal[1].split('(')[1].split(')')[0]				
					readseq = reads_dict[id]
					refseq = name_refseq[acname]

					st = ''
					p = 0
					p_ref = 0
				
					for i in tag: # for each element in string:
		
						if i not in ['1','2','3','4','5','6','7','8','9','0']:
						
							l = int(st)
							if i=='S':
								p += l
							
								st = ''
							
							elif i=='M':
								p+= l
								p_ref +=l
								st = ''
						
							elif i=='D':
								st=''
								p_ref += l+1
						
							elif i=='I':
						
								if p_ref == int(pos_l) and p_ref+1 == int(pos_u):
									ins = readseq[p:p+l]
							
								p+=l
								st=''		
						
						else:
							st += i 
				
					tally_nucl[tally] = ins		
					
		with open('tally_w_del_seq_fset_guide_cigars_{0}_{1}.txt'.format(nm, acname), 'w') as g:
				c = 0
				for l in tally_nucl.keys():
					g.write('{0}\t\t{1}\n'.format(l, tally_nucl[l]))
					c+=int(l.split('\t')[2])
				
				ref_tally[acname]=c
										
	ref_total = {}
	with open('refseq_guide_cigars_{0}.txt'.format(nm), 'r') as h:
		data = [i for i in h.read().split('\n') if len(i)>0]
		for l in data:
			l = l.split('-')[0]
			if l in ref_total.keys():
				ref_total[l] += 1
			else:
				ref_total[l]=1
		
	print ref_tally	
	print ref_total
	
	with open('tally_all_seqs_{0}.txt'.format(nm), 'w') as v:
		v.write('ref_seq\tindel events\ttotal events\tpercent\n')
		for g in ref_total.keys():
			events, total = ref_tally[g], ref_total[g]
			p = float(events)/float(total)*100.0
			v.write('{0}\t{1}\t{2}\t{3}\n'.format(g, events, total, p))
		




				
				