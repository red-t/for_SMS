import pysam
import uuid

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





samfile = pysam.AlignmentFile("test.bam","rb")
reads=[]
for read in samfile.fetch():
    reads.append(read)

read = reads[0]
