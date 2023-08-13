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
        list te_seqs = [t[1] for t in tmp] # [te_seq1, te_seq2, ...]
        SequenceContainer sc = SequenceContainer(te_seqs)
        PopGenDefinitionReader pgdr = PopGenDefinitionReader(pgdf, sc)
        list tedeftuples = pgdr.read_transposed()
        str crap, chasis
        int counter = 1 + sub_idx * sub_size
        FastaWriter fw = FastaWriter(output, 60)
        list seqidtup = [], seqtup
        str seq_with_tes
    
    # Load template(chasis) sequence
    crap, chasis = load_chasis(ref_fasta)

    # Insert te sequences into template step by step
    for tmp in tedeftuples: # [(pos1, insid1), (pos2, insid2), ...]
        # Get ids of sequence to be inserted into template
        # seqidtup = [(t[0], str(crap) + "_hg" + str(counter) + ";" + t[1]) for t in tmp]

        # Get instances of sequence to be inserted into template
        seqtup = [(t[0], sc.getTESequence(t[1])) for t in tmp] # [(pos1, TESequence), (pos2, TESequence), ...]

        # Sequence of haplotype genome with inserted sequences
        seq_with_te = insertSequences(chasis, seqtup, outseqf, seqidtup)

        # Write out the haplotype genome to population genome
        fw.write(str(crap) + "_hg" + str(counter), seq_with_te)
        counter += 1
    print("Done; Wrote {0} genomes to file {1}".format(sub_size, output))
    fw.close()