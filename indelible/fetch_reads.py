import pysam

BASE_QUALITIES = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~"

def hard_clip(seq,qual,threshold=10):
	res_seq = ""
	res_qual = ""

	
	offset = 0
	q = BASE_QUALITIES.index(qual[offset])
	while q < threshold and offset < len(seq):
		offset += 1
		if offset < len(seq):
			q = BASE_QUALITIES.index(qual[offset])

	return [seq[offset:],qual[offset:]]

def average_quality(qual):
	if len(qual) != 0:
		s = 0
		for q in qual:
			s += BASE_QUALITIES.index(q)
		return float(s)/len(qual)
	else:
		return 0


def print_if_ok(sr,config):
	outputfmt = "%(chr)s\t%(split_position)s\t%(prime)s\t%(length)s\t%(seq)s\t%(qual)s\t%(mapq)s\t%(avg_sr_qual).2f\t%(strand)s"
	if sr["length"] >= config["MINIMUM_LENGTH_SPLIT_READ"] and sr["mapq"] >= config["MINIMUM_MAPQ"] and sr["avg_sr_qual"] >= config["MININUM_AVERAGE_BASE_QUALITY_SR"]:
		return (outputfmt % sr)+"\n"
	else:
		return ""



def fetch_reads(input_path,output_path,config):
	infile = pysam.Samfile(input_path,'rb')
	outfile = open(output_path,'w')
	outfile.write("chr\tsplit_position\tprime\tsplit_length\tseq\tqual\tmapq\tavg_sr_qual\treverse_strand\n")

	for s in infile:
		cigar = s.cigar
		if len(cigar) == 2:
			#5' Split Reads
			if cigar[0][0] == 4 and cigar[1][0] == 0: # 4 = SOFT CLIP/ 0 = MATCH
				sr = {}
				sr["chr"] = infile.getrname(s.tid)
				sr["split_position"] = s.pos
				sr["prime"] = 5
				sr["seq"],sr["qual"] = hard_clip(s.seq[0:cigar[0][1]],s.qual[0:cigar[0][1]],config["HC_THRESHOLD"])
				sr["length"] = len(sr["seq"])
				sr["mapq"] = s.mapq
				sr["avg_sr_qual"] = average_quality(sr["qual"])
				sr["strand"] = s.is_reverse
				outfile.write(print_if_ok(sr,config))
			#3' Split Reads
			if cigar[0][0] == 0 and cigar[1][0] == 4:
				sr = {}
				sr["chr"]  = infile.getrname(s.tid)
				sr["split_position"] = s.pos + cigar[0][1]
				sr["prime"] = 3
				seq = s.seq[-cigar[1][1]:]
				qual = s.qual[-cigar[1][1]:]
				seq,qual = hard_clip(seq[::-1],qual[::-1],config["HC_THRESHOLD"])
				sr["seq"] = seq[::-1]
				sr["qual"] = qual[::-1]
				sr["length"] = len(sr["seq"])
				sr["mapq"] = s.mapq
				sr["avg_sr_qual"] = average_quality(sr["qual"])
				sr["strand"] = s.is_reverse
				outfile.write(print_if_ok(sr,config))
		elif len(cigar) == 3:

			# These are reads with Soft-clips on both sides, likely due to dropped quality at the end
			if cigar[0][0] == 4 and cigar[1][0] == 0 and cigar[2][0] == 4:
				#1st split-segment
				sr = {}
				sr["chr"] = infile.getrname(s.tid)
				sr["split_position"] = s.pos
				sr["prime"] = 5
				sr["seq"],sr["qual"] = hard_clip(s.seq[0:cigar[0][1]],s.qual[0:cigar[0][1]],config["HC_THRESHOLD"])
				sr["length"] = len(sr["seq"])
				sr["mapq"] = s.mapq
				sr["avg_sr_qual"] = average_quality(sr["qual"])
				sr["strand"] = s.is_reverse
				outfile.write(print_if_ok(sr,config))
				#2nd split segment
				sr = {}
				sr["chr"]  = infile.getrname(s.tid)
				sr["split_position"] = s.pos + cigar[1][1] #alignment_start + length of matching segment = start of 3' split segment
				sr["prime"] = 3
				seq = s.seq[-cigar[2][1]:]
				qual = s.qual[-cigar[2][1]:]
				seq,qual = hard_clip(seq[::-1],qual[::-1],config["HC_THRESHOLD"])
				sr["seq"] = seq[::-1]
				sr["qual"] = qual[::-1]
				sr["length"] = len(sr["seq"])
				sr["mapq"] = s.mapq
				sr["avg_sr_qual"] = average_quality(sr["qual"])
				sr["strand"] = s.is_reverse
				outfile.write(print_if_ok(sr,config))

		
