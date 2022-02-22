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