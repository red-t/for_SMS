import pysam
import numpy as np
from datetime import datetime
import time
import sys

filename = sys.argv[1]
samfile = pysam.AlignmentFile(filename,"rb")

###################################
# process_start = time.process_time()
# perf_start = time.perf_counter()
# for read in samfile.fetch():
#     positions_1 = []
#     cigar_array = np.array(read.cigartuples)
#     idx_I_op = np.where(np.isin(cigar_array[:,0], [1,4]) & (cigar_array[:, 1] >= 200)) # tuple of arrary
#     if len(idx_I_op[0]):
#         for i in range(len(idx_I_op[0])):
#             idx_current = idx_I_op[0][i]

#             if i > 0:
#                 idx_prior = idx_I_op[0][i-1]
#                 tmp_array = cigar_array[idx_prior+1:idx_current,:]
#                 start = positions_1[i-1][1] + np.sum(tmp_array[np.isin(tmp_array[:,0], [2,3,5], invert=True),1])
#                 end = start + cigar_array[idx_current][1]
#             else:
#                 tmp_array = cigar_array[0:idx_current,:]
#                 start = np.sum(tmp_array[np.isin(tmp_array[:,0], [2,3,5], invert=True),1])
#                 end = start + cigar_array[idx_current][1]
            
#             positions_1.append((start, end))

# process_end = time.process_time()
# perf_end = time.perf_counter()
# print('Proc time1: ', (process_end - process_start)*1000)
# print('Perf time1: ', (perf_end - perf_start)*1000, '\n')


###################################
process_start = time.process_time()
perf_start = time.perf_counter()
read_counts = 0
frag_counts = 0
for read in samfile.fetch():
    positions_2 = []
    b = 0
    add_read = 0

    for cig in read.cigartuples:
        if cig[1] >= 200 and cig[0] in (1,4):
            
            add_read = 1
            q_start = b-1
            if q_start < 0:
                q_start = 0
            
            q_end = b+cig[1]
            positions_2.append((q_start, q_end))
        
        if cig[0] in (0,1,4):
            b += cig[1]
    
    read_counts += add_read
    frag_counts += len(positions_2)


process_end = time.process_time()
perf_end = time.perf_counter()
print('Read counts: ', read_counts)
print('Fragments counts: ', frag_counts)
print('Proc time2: ', (process_end - process_start)*1000)
print('Perf time2: ', (perf_end - perf_start)*1000, '\n')


##################################
# process_start = time.process_time()
# perf_start = time.perf_counter()
# for read in samfile.fetch():
#     positions_3 = []
#     cigar_array = np.array(read.cigartuples)
#     idx_I_op = np.where(np.isin(cigar_array[:,0], [1,4]) & (cigar_array[:, 1] >= 200)) # tuple of arrary
#     if(len(idx_I_op[0])):
#         b = 0
#         for cig in read.cigartuples:
#             if cig[1] >= 200 and cig[0] in (1,4):
                
#                 q_start = b-1
#                 if q_start < 0:
#                     q_start = 0
                
#                 q_end = b+cig[1]
#                 positions_3.append((q_start, q_end))
            
#             if cig[0] in (0,1,4):
#                 b += cig[1]

# process_end = time.process_time()
# perf_end = time.perf_counter()
# print('Proc time3: ', (process_end - process_start)*1000)
# print('Perf time3: ', (perf_end - perf_start)*1000, '\n')





# process_start = time.process_time()
# perf_start = time.perf_counter()
# for read in samfile.fetch():
#     # read.cigartuples
#     cigar_array = np.array(read.cigartuples)
#     idx_I_op = np.where(np.isin(cigar_array[:,0], [1,4]) & (cigar_array[:, 1] >= 200)) # tuple of arrary
#     if(len(idx_I_op[0])):
#         continue

# process_end = time.process_time()
# perf_end =  time.perf_counter()
# # execute_end = time.time()
# print('Proc time4: ', (process_end - process_start)*1000000)
# print('Perf time4: ', (perf_end - perf_start)*1000000)
# print('Exec time4: ', (execute_end - execute_start)*1000)

# frac=$(awk 'BEGIN {total=17473032}')