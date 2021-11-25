import pysam
import uuid
import numpy

class CCS:
    def __init__(self):
        self.cigarstring = None
        self.cigartuples = None # list of tuples, [(operation, length), ...]
        self.flag = None

        self.query_alignment_start = None # 0-based, including.
        self.query_alignment_end = None # 0-based, excluding.
        self.query_alignment_length = None  # equal to `query_alignment_end - query_alignment_start`
        self.query_alignment_sequence = None
        self.query_alignment_qualities = None # an array of base (read) qualities

        self.query_name = None
        self.query_length = None # equal to `infer_read_length()`
        self.query_qualities = None
        self.query_sequence = None  # including soft clipped bases

        self.reference_name = None
        self.reference_start = None # 0-based, including. reference position corresponding to the first aligned base (excluding soft clipped bases)
        self.reference_end = None # 0-based, excluding. reference position corresponding to the first aligned base (excluding soft clipped bases)
        self.reference_length = None
        self.tags = None # corresponding to `read.get_tags()`

        self.flag = None
        self.is_duplicate = None
        self.is_paired = None
        self.is_proper_pair = None
        self.is_qcfail = None
        self.is_read1 = None
        self.is_read2 = None
        self.is_reverse = None
        self.is_secondary = None
        self.is_supplementary = None
        self.is_unmapped = None
    
    def infer_read_length(self): # for hard clipped alignment
        pass

    def get_tags(self, with_value_type=False):
        pass

    def get_reference_sequence(self):
        pass

    def get_reference_positions(self, full_length=False): # returned list will be of the same length as the read when setting `full_length=True`.
        pass

    def get_overlap(self, ref_start, ref_end): # excluding clipping; excluding 'D', 'N', 'P'
        pass

    def get_forward_sequence(self): # return the original read sequence. That is 'reverse complementary' --> 'original orientation'
        pass

    def get_forward_qualities(self):
        pass

    def get_cigar_stats(self):
        pass

    def get_blocks(self): # a list of (reference)start and end positions of aligned gapless blocks.
        pass

    def get_aligned_pairs(self, matches_only=False, with_seq=False): # return list of tuples; [(q_pos, r_pos), ...]
        pass


class InsertFragment:
    """ 用于表征一条alignment当中，insert或clip的片段(CIGAR operation为'I'或'S') """
    def __init__(self, read, q_is, q_ie, type, r_is, r_ie):
        self.query_name = read.query_name
        self.query_seq = read.query_sequence
        self.query_qual = read.qual
        self.query_ins_start = q_is # 0-based, including. Left-most position of the insert/clip fragment on the query sequence.
        self.query_ins_end = q_ie # 0-based, excluding. Left-most position of the insert/clip fragment on the query sequence.
        self.ins_length = self.query_ins_end - self.query_ins_start
        self.type = type

        self.ref_name = read.reference_name
        self.ref_ins_start = r_is
        self.ref_ins_end = r_ie
        self.mapping_quality = read.mapping_quality
        self.is_reverse = read.is_reverse

        self.query_seq_trimmed = None
        self.query_qual_trimmed = None
        self.query_ins_start_trimmed = None
        self.query_ins_end_trimmed = None

        self.te_alignment = None
        self.te_name = None
        self.te_mapping_quality = None
        self.te_ref_start = None
        self.te_ref_end = None
        self.te_query_start = None
        self.te_query_end = None
        self.te_orient = None

    def trim(self, flanksize=200):
        """ 将insert/clip片段从原本的query序列中剪切出来，两端各延伸flanksize长度。"""
        trimmed_start = self.query_ins_start - flanksize
        trimmed_end = self.query_ins_end + flanksize
        
        if trimmed_start < 0:
            trimmed_start = 0
        if trimmed_end > len(self.query_seq):
            trimmed_end = len(self.query_seq)
        
        self.query_seq_trimmed = self.query_seq[trimmed_start:trimmed_end]
        self.query_qual_trimmed = self.query_qual[trimmed_start:trimmed_end]
        self.query_ins_start_trimmed = self.query_ins_start - trimmed_start
        self.query_ins_end_trimmed = self.query_ins_start_trimmed + self.ins_length

    def write_fastq(self):
        return '@%s\n%s\n+\n%s\n' % (self.query_name, self.query_seq_trimmed, self.query)

    def write_fasta(self):
        return '>%s_%s\n%s\n' % (self.query_name, self.type, self.query_seq_trimmed)


class InsCandidate:
    def __init__(self, id):
        self.id = id
        self.tes = None
        self.reads = []
        self.break_points = None
        self.consensus = None
        self.tsd = None
    
    def make_consensus(self):
        """ define consensus sequence """
        pass

    def find_bp(self):
        """ find break points """
        pass

    def find_tsd():
        """ find TSD sequence """
        pass












samfile = pysam.AlignmentFile("test.bam","rb")
reads=[]
for read in samfile.fetch():
    reads.append(read)

read = reads[1]
