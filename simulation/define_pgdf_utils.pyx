import random
import re

cpdef set get_germline_pos(int ngermline, int nsomatic):
    cdef:
        set ret = set()
        int ntotal = ngermline + nsomatic - 1
        int x
    
    while len(ret) < ngermline:
        x = random.randint(0, ntotal)
        if x in ret:
            continue
        ret.add(x)
    
    return ret
#
#
#
cpdef define_nest_ins(dict idx2te, int nte, int parent_te_len, int parent_trunc_len, str te_type, float d_rate):
    # Chose TE id
    cdef:
        str te_id
        int ccs_idx, te_len
    if te_type == "s":
        # ccs_idx = random.choice([72, 15, 54, 44, 100, 106]) # for fly
        ccs_idx = random.choice([1, 2, 3, 4]) # for human
    else:
        ccs_idx = random.randint(1, nte)
    te_id = idx2te[ccs_idx][0]
    te_len = idx2te[ccs_idx][1]

    # TSD
    cdef int tsd
    tsd = random.randint(6, 8)

    # Nest position
    cdef int pos
    pos = random.randint(tsd, parent_te_len - parent_trunc_len)

    # Strand
    cdef str strand, l_strand
    if random.random() < 0.5:
        strand = "-"
        l_strand = "R"
    else:
        strand = "+"
        l_strand = "F"

    # Truncations
    cdef:
        str trunc_id = "U"
        str trunc_dsl = ""
        int trunc_len = 0
        int trunc_start
    if random.random() <= 0.01:
        trunc_id = "T"
        trunc_len = int(random.uniform(0.1, 0.6)*te_len)
        trunc_start = random.randint(1, te_len - trunc_len)
        trunc_dsl = "[{0}..{1}]".format(trunc_start, trunc_start + trunc_len - 1)

    # Nest
    cdef:
        str nest_dsl = ""
        str nest_id = ""
        str left, right
    if random.random() <= 0.01:
        nest_id, nest_dsl = define_nest_ins(idx2te, nte, te_len, trunc_len, te_type, d_rate)
        nest_id = "(" + nest_id + ")"
        nest_dsl = "{" + nest_dsl + "}"

    left = "{0}~{1}{2}{3}".format(te_id, trunc_id, l_strand, nest_id)
    if d_rate > 0:
        right = "{0}:${1}{2}{3}{4}{5}%{6}bp".format(pos, ccs_idx, trunc_dsl, strand, nest_dsl, d_rate, tsd)
    else:
        right = "{0}:${1}{2}{3}{4}{5}bp".format(pos, ccs_idx, trunc_dsl, strand, nest_dsl, tsd)
    return left, right
#
#
#
cpdef define_ins(int idx, dict idx2te, int nte, str te_type="g", float d_rate=0):
    # Chose TE id
    cdef:
        str te_id
        int ccs_idx, te_len
    if te_type == "s":
        # ccs_idx = random.choice([72, 15, 54, 44, 100, 106]) # for fly
        ccs_idx = random.choice([1, 2, 3, 4]) # for human
    else:
        ccs_idx = random.randint(1, nte)
    te_id = idx2te[ccs_idx][0]
    te_len = idx2te[ccs_idx][1]

    # TSD
    cdef int tsd
    tsd = random.randint(6,8)

    # Strand
    cdef str strand, l_strand
    if random.random() < 0.5:
        strand = "-"
        l_strand = "R"
    else:
        strand = "+"
        l_strand = "F"

    # Truncations
    cdef:
        str trunc_id = "U"
        str trunc_dsl = ""
        int trunc_len = 0
        int trunc_start
    if random.random() <= 0.1:
        trunc_id = "T"
        trunc_len = int(random.uniform(0.1, 0.6)*te_len)
        trunc_start = random.randint(1, te_len - trunc_len)
        trunc_dsl = "[{0}..{1}]".format(trunc_start, trunc_start + trunc_len - 1)

    # Nest
    cdef:
        str nest_dsl = ""
        str nest_id = ""
        str left, right
    if random.random() <= 0.1:
        nest_id, nest_dsl = define_nest_ins(idx2te, nte, te_len, trunc_len, te_type, d_rate)
        nest_id = "(" + nest_id + ")"
        nest_dsl = "{" + nest_dsl + "}"

    left = "{0}~{1}~{2}{3}{4}".format(idx, te_id, trunc_id, l_strand, nest_id)
    if d_rate > 0:
        right = "${0}{1}{2}{3}{4}%{5}bp".format(ccs_idx, trunc_dsl, strand, nest_dsl, d_rate, tsd)
    else:
        right = "${0}{1}{2}{3}{4}bp".format(ccs_idx, trunc_dsl, strand, nest_dsl, tsd)

    return left, right
#
#
#
cpdef define_header(str ref_fa, str te_fa, int germline_count, int somatic_count, set germline_pos, float d_rate):
    cdef:
        str ref_idx = ref_fa+".fai"
        str te_idx = te_fa+".fai"
        str contig_id
        dict idx2te = {}
        dict id2dsl = {}
        int genome_size

    # reference information
    for l in open(ref_idx, "r"):
        l = l.split()
        contig_id = l[0]
        genome_size = int(l[1])
        with open("{0}.tmp.pgd.header".format(contig_id), "w") as fout:
            fout.write("# Chasis {0}; Length {1} nt\n".format(contig_id, l[1]))
    
    # parse TE information
    cdef int counter = 1
    for l in open(te_idx, "r"):
        l = l.split()
        if l[0] != "":
            idx2te[counter] = (l[0], int(l[1]))
            counter += 1

    # simualte total (germline_count + somatic_count) insertions
    cdef:
        list ins_ids = []
        int ins_count = germline_count + somatic_count
        int idx, nte = len(idx2te)
        str left, right

    fout = open("{0}.tmp.pgd.header".format(contig_id), "a")
    for idx in range(1, ins_count+1):
        if (idx-1) in germline_pos:
            left, right = define_ins(idx, idx2te, nte, d_rate=d_rate)
        else:
            left, right = define_ins(idx, idx2te, nte, "s", d_rate=d_rate)
        
        ins_ids.append(left)
        id2dsl[left] = right
        fout.write("{0}={1}\n".format(left, right))
    
    fout.close()
    return id2dsl, contig_id, genome_size, ins_ids
#
#
#
cpdef define_body(dict id2dsl, int popsize, int ntotal, str contig_id, int mindist, int maxdist, set germline_pos, list ins_ids):
    cdef:
        int i, j
        int dist
        int counter = 1
        int count_te, count_empty
        float popfreq
        str id, tows, dsl, strand
        list toshuf
        int tsd_len, tsd_start

    fout = open("{0}.tmp.pgd.body".format(contig_id), "w")
    summary = open("{0}.ins.summary".format(contig_id), "w")
    # write insertions
    for i in range(ntotal):
        dist = random.randint(mindist, maxdist)
        counter += dist # in other words, "pos"
        if i in germline_pos:
            # popfreq = random.uniform(0.1, 1) # for fly
            popfreq = random.choice((0.5, 1)) # for human
            count_te = int(popfreq * popsize)
        else:
            popfreq = float(1) / popsize
            count_te = 1
        count_empty = popsize - count_te
        id = ins_ids[i]
        toshuf = [id for j in range(count_te)]
        toshuf.extend("*" * count_empty)
        random.shuffle(toshuf)
        toshuf.insert(0, str(counter))
        tows = " ".join(toshuf)
        # write out body
        fout.write(tows+"\n")

        # write out insertion information
        dsl = id2dsl[id]
        strand = re.search(r"([+-])", dsl).group(1)
        tsd_len = int(re.search("([\d\.]+)bp$", dsl).group(1))
        tsd_start = counter - tsd_len + 1
        summary.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format(contig_id, counter, counter+1, id, dsl, strand, popfreq, tsd_start, counter))
    
    fout.close()
    summary.close()