#!/usr/bin/python3

import os, re
import pdb
import subprocess as sp
from collections import Counter

""" Prefix for type checking
s:string
i:int
f:float
l:list
d:dictionary
"""

sSamtools        = "/path/to/samtools-1.6/samtools"
sBam_output      = "/path/to/Result/BWA/Sorted_sample.bam"
sMutCount_output = "/path/to/Result/MutCount_output/MutCount_output.txt"

sBam_output_dir = "/path/to/Result/BWA"
sMutFunc_output = "/path/to/Result/MutFunc_output"

dCodontable = {
	'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
	'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
	'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
	'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
	'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
	'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
	'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
	'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
	'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
	'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
	'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
	'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
	'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
	'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
	'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
	'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
}


def Count_mut():
	
	cnt = Counter()
	iTot = 0

	for i, sRow in enumerate(sp.Popen(sSamtools +" view %s" % sBam_output, stdout=sp.PIPE, shell=True).stdout.readlines()):
		lCol   = sRow.decode("utf-8").replace("\n", "").split('\t')
		iStart = int(lCol[3]) - 1
		sSeq   = lCol[9]

		sTarget_seq = sSeq[31-iStart:41-iStart]

		if len(sTarget_seq.strip()) == 10:
			cnt[sTarget_seq] += 1
			iTot += 1

	with open(sMutCount_output, 'w') as Output:
		Output.write("Total reads: " + str(iTot)+'\n')
		for sMut_type, iCount in cnt.most_common():
			Output.write(sMut_type+"\t"+str(iCount)+"\n")


def Count_mut_func():
	"""
				Without indels					With indels
	Samples	Sample explantation	WT	Original GOF	Secondary GOF	Unclassfied point mutations	Without indel total (C+D+E+F)	inframe indels	out of frame indels	With indel total (H+I)	total read (G+J)
	"""

	with open(sMutFunc_output + '/MutFunc_output.txt', 'w') as Output:
		Output.write('\t'.join(['Sample', 'WT', 'Original_GOF', 'Secondary_GOF', 'Unclassified_point_mutations',
								'Without_indel_total', 'inframe_indels', 'out_of_frame_indels',
								'With_indel_total', 'total_read', 'AA_cnt', '6bp_Sequence_cnt', 'Total_unclass_syn', 'unclass_syn_AA', 'unclass_syn_seq',
								'Total_unclass_nonsyn', 'unclass_nonsyn_AA', 'unclass_nonsyn_seq',]) + '\n')

		for sFilename in os.listdir(sBam_output_dir):
			if sFilename.split('.')[-1] != 'bam':continue
			if sFilename[:6] != 'Sorted': continue
			# print(sFilename)

			dCnt = {'WT':0, 'Original_GOF':0, 'Secondary_GOF':0, 'Unclassfied_point_mutations':0, 'Without_indel_total':0,
					   'inframe_indels':0, 'out_of_frame_indels':0, 'With_indel_total':0, 'total_read':0}
			AA_cnt    = Counter()
			Codon_cnt = Counter()

			AA_uncla_nonsyn_cnt = Counter()
			AA_uncla_syn_cnt    = Counter()

			Codon_uncla_nonsyn_cnt = Counter()
			Codon_uncla_syn_cnt = Counter()

			for i, sRow in enumerate(sp.Popen(sSamtools + " view %s" % sBam_output_dir+'/'+sFilename, stdout=sp.PIPE, shell=True).stdout.readlines()):
				lCol = sRow.decode("utf-8").replace("\n", "").split('\t')

				iStart = int(lCol[3]) - 1
				sCigar = lCol[5]
				sSeq   = lCol[9]

				iSeq_len = sum(map(int, re.findall(r'\d+', sCigar)))

				if sSeq[:3] == 'ATG' and iSeq_len >= 70:

					if len(sSeq) % 3 in [0, 2]:
						dCnt['out_of_frame_indels'] += 1
						dCnt['With_indel_total'] += 1
						#print(sSeq, sCodon, sFunc_verdict)

					elif len(sSeq) % 3 and (sCigar.find('I') != -1 or sCigar.find('D') != -1):
						dCnt['inframe_indels'] += 1
						dCnt['With_indel_total'] += 1
						# print(sSeq, sCodon, sFunc_verdict)

					else:
						dCnt['Without_indel_total'] += 1
						sSix_codon_seq    = sSeq[33 - iStart:39 - iStart].upper()
						
						if len(sSix_codon_seq) == 6:
							sDouble_codon = dCodontable[sSix_codon_seq[:3]] + dCodontable[sSix_codon_seq[3:]]

							if sSix_codon_seq == 'GGTGGC':
								dCnt['WT'] += 1

							else:

								if sSix_codon_seq in ['GATGGC', 'GTTGGC', 'TGTGGC', 'GGTGAC', 'CGTGGC', 'GCTGGC', 'AGTGGC']:
									sFunc_verdict = 'Original_GOF'

								elif sSix_codon_seq in ['GACGGC', 'GTCGGC', 'GTAGGC', 'GTGGGC', 'TGCGGC', 'CGCGGC',
														'CGAGGC', 'CGGGGC', 'AGAGGC', 'AGGGGC', 'GCCGGC', 'GCAGGC',
														'GCGGGC', 'TCTGGC', 'TCCGGC', 'TCAGGC', 'TCGGGC', 'AGCGGC',
														'GGTGAT']:
									sFunc_verdict = 'Secondary_GOF'

									AA_cnt[sDouble_codon] += 1
									Codon_cnt[sSix_codon_seq + '(%s)' % sDouble_codon] += 1

								elif '*' in sDouble_codon:
									sFunc_verdict = 'out_of_frame_indels'

								else:
									sFunc_verdict = 'Unclassfied_point_mutations'
									if sDouble_codon == 'GG':
										AA_uncla_syn_cnt[sDouble_codon] += 1
										Codon_uncla_syn_cnt[sSix_codon_seq + '(%s)' % sDouble_codon] += 1
									else:
										AA_uncla_nonsyn_cnt[sDouble_codon] += 1
										Codon_uncla_nonsyn_cnt[sSix_codon_seq + '(%s)' % sDouble_codon] += 1

								dCnt[sFunc_verdict] += 1

			iWithout_indel_total = dCnt['WT'] + dCnt['Original_GOF'] + dCnt['Secondary_GOF'] + dCnt['Unclassfied_point_mutations']
			iWith_indel_total    = dCnt['inframe_indels'] + dCnt['out_of_frame_indels']
			iTotal_read          = dCnt['Without_indel_total'] + dCnt['With_indel_total']

			sAA_cnt    = ','.join([i[0]+':'+str(i[1]) for i in AA_cnt.most_common()])
			sCodon_cnt = ','.join([i[0]+':'+str(i[1]) for i in Codon_cnt.most_common()])

			sAA_uncla_syn_cnt       = ','.join([i[0] + ':' + str(i[1]) for i in AA_uncla_syn_cnt.most_common()])
			sCodon_uncla_syn_cnt    = ','.join([i[0] + ':' + str(i[1]) for i in Codon_uncla_syn_cnt.most_common()])
			iTotal_AA_uncla_syn_cnt = sum(i[1] for i in AA_uncla_syn_cnt.most_common())

			sAA_uncla_nonsyn_cnt       = ','.join([i[0] + ':' + str(i[1]) for i in AA_uncla_nonsyn_cnt.most_common()])
			sCodon_uncla_nonsyn_cnt    = ','.join([i[0] + ':' + str(i[1]) for i in Codon_uncla_nonsyn_cnt.most_common()])
			iTotal_AA_uncla_nonsyn_cnt = sum(i[1] for i in AA_uncla_nonsyn_cnt.most_common())

			Output.write('\t'.join(map(str, [sFilename.replace('Sorted_', '').replace('.bam', ''),
											 dCnt['WT'], dCnt['Original_GOF'], dCnt['Secondary_GOF'], dCnt['Unclassfied_point_mutations'],
											 iWithout_indel_total, dCnt['inframe_indels'], dCnt['out_of_frame_indels'], iWith_indel_total,
											 iTotal_read, sAA_cnt, sCodon_cnt, iTotal_AA_uncla_syn_cnt, sAA_uncla_syn_cnt, sCodon_uncla_syn_cnt,
											 iTotal_AA_uncla_nonsyn_cnt, sAA_uncla_nonsyn_cnt, sCodon_uncla_nonsyn_cnt])) + "\n")


def Main():
	Count_mut()
	Count_mut_func()

if __name__ == '__main__':
	Main()