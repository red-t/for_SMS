#!/usr/bin/python
import argparse
from TEMP3.ParallelTemplate import build_cluster_parallel

def parse_args():
    parser = argparse.ArgumentParser(description="Demo of argparse")

    parser.add_argument('-b', '--bam', dest='fpath', type=str,
                                 help='path of input BAM file, mapped by minimap2 -Y', default='')
    parser.add_argument('-o', '--out_path', dest='out_path', type=str,
                                 help='Path of the output', default='./')
    parser.add_argument('-p', '--nprocess', dest='nprocess', type=int,
                                 help='maximum nummber of worker processes.', default=1)
    parser.add_argument('-t', '--nthreads', dest='nthreads', type=int,
                                 help='maximum nummber of extra threads to use in each sub process.', default=1)
    parser.add_argument('-l', '--minl', dest='minl', type=int,
                                 help='minimun segment length, cigar operation with length < minl will not be used to create insert segment.', default=50)
    parser.add_argument('-L', '--maxdist', dest='maxdist', type=int,
                                 help='max merging distance, segments with distance larger than maxdist will not be merged in to the same cluster.', default=50)
    
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    # parse args from command line options
    args = parse_args()

    # build clusters in parallel
    tid_to_clusters = build_cluster_parallel(args.fpath,
                                             args.nprocess,
                                             args.nthreads,
                                             args.minl,
                                             args.maxdist)
    
    # next