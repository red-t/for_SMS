import random
import pysam

cpdef generate_PE(str pop_gen,
                  float error_rate,
                  int nread,
                  int insertsize,
                  int stddev,
                  int rlen,
                  str fq1,
                  str fq2,
                  float chimera):
    cdef:
        list pgld = get_length_list(pop_gen)
        NGS_Mutator mutator = NGS_Mutator(error_rate)
        RandomReads randr = RandomReads(nread, pgld)
        FastqPairWriter fqwriter = FastqPairWriter(fq1, fq2)
        int counter = 0
        int readcount = 1
        str header, seq
        int targetreads, seqlen, i
        int firstposition, secondposition, innerdistance, outerdistance
        float ischimera
        str h, read1, read2

    with pysam.FastxFile(pop_gen) as fh:
        for entry in fh:
            header = entry.name
            seq = entry.sequence
            targetreads = randr.get_reads(counter)
            print("Generating {0} reads for haploid genome {1}".format(targetreads, header))
            counter += 1
            seqlen = len(seq)
            for i in range(0, targetreads):
                ischimera = random.random()
                if(ischimera < chimera): # read is a chimera
                    firstposition = random.randint(0, seqlen - rlen)
                    secondposition = random.randint(0, seqlen - rlen)
                else:                   # read is normal paire-end
                    innerdistance = int(random.gauss(insertsize, stddev))
                    while(innerdistance < 1): # must be larger than 0
                        innerdistance = int(random.gauss(insertsize, stddev))

                    outerdistance = 2 * rlen + innerdistance
                    firstposition = random.randint(0, seqlen - outerdistance)
                    secondposition = firstposition + rlen + innerdistance
                read1 = seq[firstposition:firstposition + rlen]
                read2 = rc(seq[secondposition:secondposition + rlen])
                read1 = mutator.mutateseq(read1)
                read2 = mutator.mutateseq(read2)
                
                h = "{0};{1}:{2}-{3}".format(readcount, header, firstposition, secondposition)
                fqwriter.write(h, read1, read2)
                readcount += 1
        
    fqwriter.close()