
import re
import pysam
from operator import itemgetter

data_file = '/data/tusers/boxu/lrft2/result/40/SMS/temp/dm3.transposon_for_simulaTE_chr2L.tmp.bam'
te_reference = 'FBgn0000155_roo'
te_reference = 'FBgn0067381_hopper2'
from collections import Counter

bam_obj = pysam.AlignmentFile(data_file, 'rb')
consensus_file = open('/data/tusers/boxu/lrft2/result/40/SMS/temp/chr2L/polished.fa','w')


con_sequence = ''
for pileup_col in bam_obj.pileup(te_reference):
    print(pileup_col.pos)
    base_col = pileup_col.get_query_sequences(mark_matches=False, mark_ends=False, add_indels=True)

    anchor_base_list = []
    indel_base_list = []
    for base in base_col:
        anchor_base_list.append(base[0].upper())
        if '+' in base:
            # indel_base_list.append(base.split('+')[1])

            indel_base_list.append(re.findall( '[ATCG]+', base.split('+')[1].upper() )[0] )


    print(">>>")
    print(pileup_col.pos)
    # print(base_col)
    print(anchor_base_list)
    print(indel_base_list)

    # base in this position
    # anchor_base_list
    map_base_list_dict = Counter(anchor_base_list)
    map_base_list = [(x,map_base_list_dict[x]) for x in map_base_list_dict ]
    map_base_list_sorted = sorted(map_base_list, key = itemgetter(1), reverse=True)
    if map_base_list_sorted[0][0] != "*":
        con_sequence = con_sequence + map_base_list_sorted[0][0].split('-')[0]
    
    # indel after this position
    # indel_base_list
    if len(indel_base_list) == 0:
        continue
    map_base_list_dict = Counter(indel_base_list)
    map_base_list = [(x,map_base_list_dict[x]) for x in map_base_list_dict ]
    map_base_list_sorted = sorted(map_base_list, key = itemgetter(1), reverse=True)
    if map_base_list_sorted[0][0] != "*" and map_base_list_sorted[0][1] >= 3:
        con_sequence = con_sequence + map_base_list_sorted[0][0].split('-')[0]



consensus_file.write('>consensus\n' + con_sequence )
consensus_file.close()



