#!/usr/bin/python
import argparse
from TEMP3.ParallelTemplate import buildClusterParallel

def parseArgs():
    parser = argparse.ArgumentParser(description="Demo of argparse")

    parser.add_argument('-b', '--bam', dest='genomeBamFilePath', type=str,
                                 help='path of input BAM file, mapped by minimap2 -Y', default='')
    parser.add_argument('-r', '--repeat', dest='repeatPath', type=str,
                                 help='repeat annotation file', default='')
    parser.add_argument('-g', '--gap', dest='gapPath', type=str,
                                 help='gap annotation file', default='')
    parser.add_argument('-T', '--refte', dest='referenceTe', type=str,
                                 help='transposon reference fasta file', default='')
    parser.add_argument('-o', '--outpath', dest='outPath', type=str,
                                 help='output directory', default='./')
    parser.add_argument('-p', '--numprocess', dest='numProcess', type=int,
                                 help='max number of worker processes.', default=1)
    parser.add_argument('-t', '--numthread', dest='numThread', type=int,
                                 help='max number of extra threads to use in each sub process.', default=1)
    parser.add_argument('-l', '--minSegLen', dest='minSegLen', type=int,
                                 help='min segment length, cigar operation with length < minSegLen will not be used to create segment.', default=100)
    parser.add_argument('-L', '--maxdist', dest='maxDistance', type=int,
                                 help='max merging distance, segments with distance larger than maxdist will not be merged in to the same cluster.', default=50)
    
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parseArgs()

    tidToResult = buildClusterParallel(args)

    # Log
    print("当前版本构建的 clusters 总数:")
    numCluster = 0
    for i in range(len(tidToResult)):
        numCluster += tidToResult[i][0].shape[0]

    print(numCluster, "\n")

    print("当前版本构建的 segments 总数:")
    numSeg = 0
    for i in range(len(tidToResult)):
        numSeg += tidToResult[i][1].shape[0]
    
    print(numSeg)