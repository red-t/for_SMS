cpdef build_popg(str te_fasta,
                 str pgdf,
                 str ref_fasta,
                 int sub_idx,
                 int sub_size,
                 str output,
                 str outseqf):
    cdef:
        tuple t
        list tmp = readAllTuples(te_fasta)
        list tetuples = [t[1] for t in tmp]
        SequenceContainer sc = SequenceContainer(tetuples)
        PopGenDefinitionReader pgdr = PopGenDefinitionReader(pgdf, sc)
        list tedeftuples = pgdr.read_transposed()
        str crap, chasis
        int counter = 1 + sub_idx * sub_size
        FastaWriter fw = FastaWriter(output, 60)
        list seqidtup = [], seqtup
        str seq_with_tes
    
    crap, chasis = load_chasis(ref_fasta)
    for tmp in tedeftuples: # [(pos1, insid1),(pos2, insid2), ...]
        # seqidtup = [(t[0], str(crap) + "_hg" + str(counter) + ";" + t[1]) for t in tmp]
        seqtup = [(t[0], sc.getTESequence(t[1])) for t in tmp] # [(pos1, TESequence),(pos2, TESequence), ...]
        seq_with_te = insertSequences(chasis, seqtup, outseqf, seqidtup)
        fw.write(str(crap) + "_hg" + str(counter), seq_with_te)
        counter += 1
    print("Done; Wrote {0} genomes to file {1}".format(sub_size, output))
    fw.close()