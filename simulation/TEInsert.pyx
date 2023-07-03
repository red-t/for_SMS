cdef class TESequence:
    def __init__(self, str sequence, int id, int tsd):
        self.sequence = sequence
        self.id = id
        self.tsd = tsd


cdef str insertSequence(str ref, int pos, TESequence toinsert, fout="", str seqid=""):
    cdef:
        str seq = toinsert.sequence
        int tsd = toinsert.tsd
        int tlp = pos - tsd

    if(tlp < 0):
        raise Exception("Invalid position of TE; Insertion position minus TSD must be larger than zero")

    cdef:
        str left = ref[:tlp] ## hmm guess we need tlp instead of pos
        int trp = tlp + tsd
        str tsdseq = ref[tlp:trp]
        str right = ref[trp:] ## guess we need trp instead of pos
        int p_left, p_right
        str t_left, t_right

    # write out simulated insertion sequence
    if fout:
        p_left = 0
        p_right = len(ref)-1
        
        if tlp - 2000 > 0:
            p_left = tlp - 2000
        if trp + 2000 < len(ref) - 1:
            p_right = trp + 2000
        
        t_left = ref[p_left:tlp]
        t_right = ref[trp:p_right]
        fout.write("{0}\t{1}\n".format(seqid, t_left + tsdseq + seq + tsdseq + t_right))
    
    return left + tsdseq + seq + tsdseq + right


cdef str insertSequences_1(str ref, list posinstuples, str fout_name="", list seqidtups=[]):
    """
    ref: reference genome
    posinstuples
    [(1,TESequence), (100, TESEquence)....]
    """
    cdef:
        list tmp, idtmp
        str seq
        int i, pos
        TESequence toins
    
    if fout_name:
        fout = open(fout_name, "a")
        tmp = sorted(posinstuples, key=lambda i:-i[0])
        idtmp = sorted(seqidtups, key=lambda i:-i[0])
        seq = ref
        i = 0
        for pos, toins in tmp:
            seq = insertSequence(seq, pos, toins, fout, idtmp[i][1])
            i += 1
        
        fout.close()
    else:
        tmp = sorted(posinstuples, key=lambda i:-i[0])
        seq = ref
        for pos, toins in tmp:
            seq = insertSequence(seq, pos, toins)
    
    return seq


cpdef str insertSequences(str ref, list posinstuples, str fout_name="", list seqidtups=[]):
    return insertSequences_1(ref, posinstuples, fout_name, seqidtups)