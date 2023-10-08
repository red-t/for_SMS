#!/usr/bin/python
import argparse
from TEMP3.ParallelTemplate import build_cluster_parallel

def parse_args():
    parser = argparse.ArgumentParser(description="Demo of argparse")

    parser.add_argument('-b', '--bam', dest='fpath', type=str,
                                 help='path of input BAM file, mapped by minimap2 -Y', default='')
    parser.add_argument('-r', '--rep', dest='rep_path', type=str,
                                 help='path of repeats annotation file', default='')
    parser.add_argument('-g', '--gap', dest='gap_path', type=str,
                                 help='path of gap annotation file', default='')
    parser.add_argument('-T', '--teref', dest='teref', type=str,
                                 help='path of reference transposon (fa/mmi)', default='')
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
    parser.add_argument('-i', '--reftid', dest='reftid', type=int,
                                 help='specify a tid, which will be used to estimate background information', default=0)
    
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    # parse args from command line options
    args = parse_args()

    # build clusters in parallel
    tid_to_clusters = build_cluster_parallel(args.fpath,
                                             args.rep_path,
                                             args.gap_path,
                                             args.teref,
                                             args.nprocess,
                                             args.nthreads,
                                             args.minl,
                                             args.maxdist,
                                             args.reftid)
    
    # next
    print("当前版本构建的 clusters 总数:")
    n = 0
    for i in range(len(tid_to_clusters)):
        n += tid_to_clusters[i][0].shape[0]

    print(n, "\n")

    print("当前版本构建的 segments 总数:")
    m = 0
    for i in range(len(tid_to_clusters)):
        m += tid_to_clusters[i][1].shape[0]
    
    print(m)