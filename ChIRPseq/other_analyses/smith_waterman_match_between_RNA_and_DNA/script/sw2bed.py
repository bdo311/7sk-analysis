'''
parse output of smith_waterman to bed file


'''
# Process string like 'mm9_ct_UserTrack_3545_MACS_peak_803 range=chr1:37485779-37487816 5'pad=0 3'pad=0 strand=+ repeatMasking=none' to get 'mm9', 'chr1', '37485779', '37487816', '+'
def processSubject(subject):
	info = subject.split(' ')
	genome = info[0].split('_')[0]
	
	posInfo = info[1][6:]
	chrom = posInfo.split(':')[0]
	
	startEnd = posInfo.split(':')[1]
	start = startEnd.split('-')[0]
	end = startEnd.split('-')[1]
	
	strand = info[4][-1]
	
	return genome, chrom, start, end, strand
	
def parseOutput(filename,outfile,rc):
	f = open(filename)
	out = open(outfile,'w')
	#out.write('chrom\tgenomic_start\tgenomic_end\tname\tscore\tstrand\tquery_start\tquery_end\tquery_len\tquery_seq\tsubject_start\tsubject_end\tsubject_len\tsubject_seq\n')
	records = []
	n = -1
	for line in f:
		if ',' in line:
			n = 0
			records = [line.strip()]
		elif n > -1:
			n = n + 1
			records.append(line.strip())
		if n == 4:
			# an entire record
			query,query_len,subject,subject_len = records[0].split(',')
			# get target genomic position
			#print subject

			#genome,chrom,start,end,strand = subject.split('_')

			# above line only works for Hs so edit by BTD for MM 3/7/14
			genome, chrom, start, end, strand = processSubject(subject)
			
			score,scorevalue,i,ivalue,j,jvalue = records[1].split('\t')
			seq1 = records[2]
			seq2 = records[4]
			if not rc:
				strand = '-'
				start1 = str(int(ivalue)-len(seq1.replace('-','')))
				end1 = ivalue 
			else:
				strand= '+'
				start1 = str(int(query_len) - int(ivalue))
				end1 = str(int(start1) + len(seq1.replace('-','')))
			out.write('\t'.join([chrom,str(int(start)+int(jvalue)-len(seq2.replace('-',''))),str(int(start)+int(jvalue)),query,scorevalue,strand,subject,start1,end1,query_len,seq1,str(int(jvalue)-len(seq2.replace('-',''))),jvalue,subject_len,seq2])+'\n')
	f.close()
	out.close()
import sys
parseOutput(sys.argv[1],sys.argv[1]+'.bed',sys.argv[2]=='rc')
