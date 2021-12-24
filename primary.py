from pysam import AlignmentFile
from concurrent.futures import ThreadPoolExecutor, as_completed
from read_alignment import build_cluster, process_cluster

def main(bam='test.bam'):
    chrom2clusters = dict()
    with ThreadPoolExecutor(max_workers=5) as executor:
        bam_file = AlignmentFile(bam, 'rb')
        chroms = list(bam_file.references)

        future2chrom = {executor.submit(build_cluster, bam, chrom):chrom for chrom in chroms}
        for future in as_completed(future2chrom):
            chrom = future2chrom[future]
            try:
                chrom2clusters[chrom] = future.result()
            except Exception as exc:
                print('building cluster for %r generated an exception: %s' % (chrom, exc))
        
        future2chrom = {executor.submit(process_cluster, chrom2clusters[chrom], chrom, "ALUL1SVA.mmi"):chrom for chrom in chroms}
        for future in as_completed(future2chrom):
            chrom = future2chrom[future]
            try:
                chrom2clusters[chrom] = future.result()
            except Exception as exc:
                print('processing cluster for %r generated an exception: %s' % (chrom, exc))
    
    return chrom2clusters


chrom2c = main()





# def t_return(s):
#     return s

# def tt():
#     a2b = {}
#     pb2a = {}
#     with ThreadPoolExecutor(max_workers=5) as executor:
#         pa2a = {executor.submit(t_return, chr(ord('A')+a)):a for a in range(4)}
#         for pa in as_completed(pa2a):
#             a = pa2a[pa]
#             b = pa.result()
#             a2b[a] = b
#             pb2a[executor.submit(t_return, b)] = a
        
#         print(pa2a)
        
#         for pb in as_completed(pb2a):
#             a = pb2a[pb]
#             b = pb.result()
#             print(a, b)
        
#         print(pb2a)


# def t_return(s):
#     return s

# def tt():
#     a2b = {}
#     with ThreadPoolExecutor(max_workers=5) as executor:
#         future2a = {executor.submit(t_return, chr(ord('A')+a)):a for a in range(4)}
#         for future in as_completed(future2a):
#             a2b[future2a[future]] = future.result()
    
#         print(future2a, a2b)

#         future2a = {executor.submit(t_return, a2b[a]):a for a in range(4)}
#         for future in as_completed(future2a):
#             a = future2a[future]
#             b = future.result()
#             print(a,b)
        
#         print(future2a)


# ref_n = list(bam_file.references)
# from pysam import AlignmentFile
# samfile = AlignmentFile("test.bam","rb")
# reads=[]
# for read in samfile.fetch():
#     reads.append(read)

# read = reads[0]

# tree = Intersecter()
# d={}
# segs = [15, 20, 30, 40]
# for seg in segs:
#     c_id = None
#     c_ids = [x.value for x in tree.find(seg-1, seg+1)]
#     print(seg, c_ids, type(c_id))
#     if c_ids:
#         for c_id in c_ids:
#             d[c_id].append(seg)
#     else:
#         c_id = str(uuid4())
#         d[c_id] = [seg]
#         tree.add_interval(Interval(seg-1, seg+1, value=c_id))

# a=[(x.start, x.end, x.value) for x in tree.find(0,100)]


# from uuid import uuid4
# # from importlib import reload
# from pysam import AlignmentFile
# from read_alignment import build_cluster, process_cluster
# from concurrent.futures import ThreadPoolExecutor, as_completed

# bam = "nanovar_ds.bam"
# clusters = {}
# with ThreadPoolExecutor(max_workers=5) as executor:
#     bam_file = AlignmentFile(bam, 'rb')
#     chroms = list(bam_file.references)
#     tasks = {executor.submit(build_cluster, bam, chrom):chrom for chrom in chroms}
#     for future in as_completed(tasks):
#         chrom = tasks[future]
#         clusters[chrom] = future.result()








# class CCS:
#     def __init__(self):
#         self.cigarstring = None
#         self.cigartuples = None # list of tuples, [(operation, length), ...]
#         self.flag = None

#         self.query_alignment_start = None # 0-based, including.
#         self.query_alignment_end = None # 0-based, excluding.
#         self.query_alignment_length = None  # equal to `query_alignment_end - query_alignment_start`
#         self.query_alignment_sequence = None
#         self.query_alignment_qualities = None # an array of base (read) qualities

#         self.query_name = None
#         self.query_length = None # equal to `infer_read_length()`
#         self.query_qualities = None
#         self.query_sequence = None  # including soft clipped bases

#         self.reference_name = None
#         self.reference_start = None # 0-based, including. reference position corresponding to the first aligned base (excluding soft clipped bases)
#         self.reference_end = None # 0-based, excluding. reference position corresponding to the first aligned base (excluding soft clipped bases)
#         self.reference_length = None
#         self.tags = None # corresponding to `read.get_tags()`

#         self.flag = None
#         self.is_duplicate = None
#         self.is_paired = None
#         self.is_proper_pair = None
#         self.is_qcfail = None
#         self.is_read1 = None
#         self.is_read2 = None
#         self.is_reverse = None
#         self.is_secondary = None
#         self.is_supplementary = None
#         self.is_unmapped = None
    
#     def infer_read_length(self): # for hard clipped alignment
#         pass

#     def get_tags(self, with_value_type=False):
#         pass

#     def get_reference_sequence(self):
#         pass

#     def get_reference_positions(self, full_length=False): # returned list will be of the same length as the read when setting `full_length=True`.
#         pass

#     def get_overlap(self, ref_start, ref_end): # excluding clipping; excluding 'D', 'N', 'P'
#         pass

#     def get_forward_sequence(self): # return the original read sequence. That is 'reverse complementary' --> 'original orientation'
#         pass

#     def get_forward_qualities(self):
#         pass

#     def get_cigar_stats(self):
#         pass

#     def get_blocks(self): # a list of (reference)start and end positions of aligned gapless blocks.
#         pass

#     def get_aligned_pairs(self, matches_only=False, with_seq=False): # return list of tuples; [(q_pos, r_pos), ...]
#         pass
























# Read name = 31e7dde4-0f61-46de-bc44-ec349ffc5639_chrTEST_90777_0_4628_r
# Read length = 4,828bp
# Flags = 16
# ----------------------
# Mapping = Primary @ MAPQ 60
# Reference span = LINE1:4-3,422 (-) = 3,419bp
# Cigar = 203S5M1D37M2I11M2I5M1D7M2I63M1I29M1I4M2I11M3I1M1I31M2D19M6I6M1I28M2I46M3D6M1I7M3I13M2I17M4I5M4D37M1I34M2D7M3I4M3I3M3D2M1I8M1I15M1D19M1I4M1D6M3D24M1I10M2I10M1I2M1I20M1I28M1I11M1I6M2D43M4D19M1I24M4I1M1D12M1I8M2D8M1D13M5I3M1I18M2I3M3D15M1D8M2I5M1D5M3I3M1I45M1I17M2I30M1I1M1I28M1I7M1D1M1D21M1D1M1D9M1I11M1D23M1D36M1I2M1I29M1D2M1D18M2I11M1D26M3D4M3D2M1D12M1I18M1I19M1I6M3I2M1I15M1I19M2I5M1I24M4D5M1I3M1D5M2D7M1D51M1I23M1I45M3I2M1D3M1D6M1I18M3D42M3I16M1I13M4D45M1D13M1D56M1I4M2D16M1D9M2D8M4D3M1I1M2I2M2I11M...7M3D5M1D7M1D14M1D9M1D7M1D15M1I12M1D8M7D10M1I26M1D1M1D2M1I6M5D5M1I6M4D1M1D38M1I5M1D11M1D4M1D3M2D4M1D41M1I6M1D23M4D8M1I20M4D2M1D5M2I5M1D8M1D4M1D10M1I5M2D4M1D6M1I10M1D4M1D12M1D5M2I5M1I19M1I15M1D6M1I1M2D4M3D10M2I6M2I8M2D5M1D9M3I6M1D14M3D6M1D4M1D5M2D16M3I18M2D9M5D3M3D22M1D18M1I8M1I23M2D10M1D3M1D10M4D2M9D8M1D3M1I7M1D3M4D1M10D12M3I8M5D13M1D7M2D8M2D2M1I1M2D11M1I7M2I1M1D12M2D8M10D8M1D13M2I7M5D7M1I10M1D6M2D15M1I6M1D19M1I2M1D1M1D7M1D3M1D3M1D1M3I12M1D1M1D23M1D2M1D3M1D29M1D11M3D4M1D13M1D4M2D4M3D6M2I18M1293S
# Clipping = Left 203 soft; Right 1,293 soft
# ----------------------
# s1 = 556
# s2 = 0
# NM = 781
# AS = 2412
# de = 0.1694
# rl = 0
# cm = 63
# nn = 0
# tp = P
# ms = 2629Location = LINE1:1,967
# Base = C @ QV 126



