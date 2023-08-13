import random
import pysam

cpdef generate_PACBIO(str pop_gen,
                      int nread,
                      int minl,
                      int maxl,
                      str outfa,
                      str protocol):
    cdef:
        list pgld = get_length_list(pop_gen)
        TGS_Mutator mutator = TGS_Mutator(protocol)
        RandomReads randr = RandomReads(nread, pgld)
        RLDfactory rldf = RLDfactory(protocol)
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
                # read length
                readlen = rldf.gamma_nextl()
                while readlen < minl or readlen > maxl:
                    readlen = rldf.gamma_nextl()
                
                if seqlen - readlen < minl:
                    readlen = random.randint(minl, seqlen - minl)
                
                # position
                firstposition = random.randint(0, seqlen - readlen)
                read1 = seq[firstposition:firstposition + readlen]

                # reverse complementary
                if random.random() < 0.5:
                    read1 = rc(read1)
                
                # mutation
                read1 = mutator.mutateseq(read1)

                # output
                h = "{0};{1}:{2}".format(readcount, header, firstposition)
                fw.write(h, read1)
                readcount += 1
        
    fw.close()



cpdef generate_ONT(str pop_gen,
                   int nread,
                   int minl,
                   int maxl,
                   str outfa,
                   str protocol):
    cdef:
        list pgld = get_length_list(pop_gen)
        TGS_Mutator mutator = TGS_Mutator(protocol)
        RandomReads randr = RandomReads(nread, pgld)
        RLDfactory rldf = RLDfactory(protocol)
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
                # read length
                readlen = rldf.expon_nextl()
                while readlen < minl or readlen > maxl:
                    readlen = rldf.expon_nextl()
                
                if seqlen - readlen < minl:
                    readlen = random.randint(minl, seqlen - minl)
                
                # position
                firstposition = random.randint(0, seqlen - readlen)
                read1 = seq[firstposition:firstposition + readlen]

                # reverse complementary
                if random.random() < 0.5:
                    read1 = rc(read1)
                
                # mutation
                read1 = mutator.mutateseq(read1)

                # output
                h = "{0};{1}:{2}".format(readcount, header, firstposition)
                fw.write(h, read1)
                readcount += 1
        
    fw.close()