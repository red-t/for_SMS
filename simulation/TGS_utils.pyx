import random
import pysam

cpdef generate_TGS(str pop_gen,
                   float error_rate,
                   str err_frac,
                   int nread,
                   int rlen,
                   float alpha,
                   float loc,
                   float beta,
                   str rldfile,
                   str outfa,
                   int tgs_minl,
                   int tgs_maxl):
    cdef:
        list pgld = get_length_list(pop_gen)
        PacBioMutator mutator = PacBioMutator(error_rate, err_frac)
        RandomReads randr = RandomReads(nread, pgld)
        RLDfactory_gamma rldf = get_rld_factory(rlen, alpha, loc, beta, rldfile)
        FastaWriter fw = FastaWriter(outfa, 60)
        int counter = 0
        int readcount = 1
        str header, seq
        int targetreads, seqlen, i, readlen, firstposition
        str h, read1
    
    with pysam.FastxFile(pop_gen) as fh:
        for entry in fh:
            header = entry.name
            seq = entry.sequence
            targetreads = randr.get_reads(counter)
            counter += 1
            seqlen = len(seq)
            print("Generating {0} reads for haploid genome {1}".format(targetreads, header))
            for i in range(0, targetreads):
                readlen = rldf.nextl()
                if readlen < tgs_minl:
                    readlen = random.randint(tgs_minl, tgs_maxl)
                if readlen > tgs_maxl:
                    readlen = random.randint(tgs_minl, tgs_maxl)

                firstposition = random.randint(0, seqlen - readlen)
                read1 = seq[firstposition:firstposition + readlen]
                if random.random() < 0.5:
                    read1 = rc(read1)
                read1 = mutator.mutateseq(read1)
                
                h = "{0};{1}:{2}".format(readcount, header, firstposition)
                fw.write(h, read1)
                readcount += 1
        
    fw.close()