'''
Input is file with ref name and cigar
'''

import re
import glob, os
import math
from Bio.Seq import Seq

# add .txt to files

for file in glob.glob('aligned_S*'):
	if '.txt' not in file:
		os.rename(file, file+'.txt')

for file in glob.glob('refseq_guide_cigars_S*'):
	if '.txt' not in file:
		os.rename(file, file+'.txt')


for file in glob.glob('refseq_guide'):
	if '.txt' not in file:
		os.rename(file, file+'.txt')

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
		
		print pos
		
		sites_of_interest[ref] = pos
		
print sites_of_interest
	
names = []

for file in glob.glob("refseq_guide_cigars_*.txt"):
	name = file.split(".")[0]
	names.append(name)


for name in names:
	ref = []
	reinge = []
	event = []
	events_dict = {}
	## make an event dict per ac number
	for acnumber in sites_of_interest.keys():
		events_dict[acnumber] = {}
	
	with open('{0}.txt'.format(name), 'r') as f:
		data = [i for i in f.read().split('\n') if len(i)>0]

	for line in data:
		
		key = line.split('-')[0]
		soi = sites_of_interest[key] #for endogenous locus
		p = 0 #position on refseq
		del_pos = []
		l_range = 0
		u_range = 0
		
		ins =''
		e = ' '
		d = 100
		d_u = 100
		d_l = 100
		best_close_l, best_close_u = ['NA', 'NA']
		best_dels = ['NA']

		tag = line.split('-')[1]
		st = ''
		tag = tag.strip('\t')
		tag = tag.strip(' ')
		tag = tag.strip('\r')
		
		print tag
		
		search_here = [i for i in xrange(soi-8, soi+12)]
		
		for i in tag: # for each element in string:
		
			if i not in ['1','2','3','4','5','6','7','8','9','0']:
				l = int(st)
		#		print l
	
		#eg:
		#	49S 89M		20D			14M			5I24M
		#	  	89 		89-109	   	103-113		site_of_interest	102-103		102-115
	
				if i =='S':
					pass
				if i == 'M':
					p+=l
				if i == 'D':
					for k in xrange(l):
						del_pos.append((p+k+1))
						
					p += l+1	

					for delpos in del_pos:
						if delpos in search_here:
							if math.fabs(delpos-soi) <= d:
								d = math.fabs(delpos-soi)
								best_dels = del_pos
								e = 'deletion'
							
					st = ''

				if i =='I':
					l_range = p
					u_range = p+1
					if u_range in search_here or l_range in search_here:
						if math.fabs(u_range-soi) <= d or math.fabs(l_range-soi) <= d:
							best_close_u = u_range
							d_u = math.fabs(u_range-soi)
							d_l = math.fabs(l_range-soi)
							d = min(d_u, d_l)
							best_close_l = l_range
							best_close_u = u_range
							ins+=str(st)
							e = 'insertion'
						st = ''
								
				st = ''
			else:
				st+=i

		if e == 'deletion':
			reinge.append(best_dels)
		else:
			reinge.append([best_close_l, best_close_u, ins])
		event.append(e)
		ref.append(line)

	with open('tagged+ranges_ch_range_{0}.txt'.format(name), 'w') as g:
		g.write('\tseq and cigar \t\t\tdeletion/insertion range\t\t event \n')
		for i in xrange(len(ref)):
			layn = ref[i]
			layn =layn.strip('\r \t')
			e = event[i]
			r = ''
			if e == 'insertion':
				r += str(reinge[i][0])+'-'+str(reinge[i][1])+' ({0})'.format(str(reinge[i][2]))
			elif e == 'deletion':
				for j in reinge[i]:
					r += str(j)+','
			g.write('{0}\t\t\t{1}\t\t\t{2}\n'.format(layn, r,e))
		
	for acnumber in events_dict.keys():
		with open('tally_{0}_{1}.txt'.format(name, acnumber), 'w') as z:
			for i in xrange(len(ref)):
				layn = ref[i]
				if acnumber == layn.split('-'[0])[0]:
					layn =layn.strip('\r \t')
					e = event[i]
					r = ''
				
					if e == 'insertion':
						r += '{0}\t'.format(layn)+str(reinge[i][0])+'-'+str(reinge[i][1])+' ({0})'.format(str(reinge[i][2]))
					elif e == 'deletion':
						for j in reinge[i]:
							r +=str(j)+','
						r = '{0}\t'.format(layn)+r
					else:
						r = '{0}_no_event\t '.format(layn.split('-'[0])[0])
						
					if r in events_dict[acnumber].keys():
						events_dict[acnumber][r]+=1
					elif r not in events_dict.keys():
						events_dict[acnumber][r] = 1	
									
			for i in events_dict[acnumber].keys():
				z.write('{0}\t{1}\n'.format(i, events_dict[acnumber][i]))
	
"""

			102  ...   ...	108   109 10 11 12 13 14 15 16 17 18 19 20
RefPos:     1  2  3  4  5  6  7     8  9 10 11 12 13 14 15 16 17 18 19
Reference:  C  C  A  T  A  C  T     G  A  A  C  T  G  A  C  T  A  A  C
Read:                   A  C  T  A  G  A  A     T  G  G  C  T	  A	 C

POS: 5
CIGAR: 3M	1I	3M	1D	5M 1D 2M
	  108M	1I	3M	1D	5M 1D 2M

		I: 108-109; D: 112, 118

	
"""