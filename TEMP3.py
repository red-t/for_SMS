#!/usr/bin/python
import os
import argparse
from TEMP3.ParallelModule import runInParallel

def parseArgs():
    parser = argparse.ArgumentParser(description="Demo of argparse")

    parser.add_argument('-b', '--bam', dest='genomeBamFilePath', type=str,
                                 help='path of input BAM file, mapped by minimap2 -Y', default='')
    parser.add_argument('-r', '--repeat', dest='repeatFn', type=str,
                                 help='repeat annotation file', default='')
    parser.add_argument('-g', '--gap', dest='gapFn', type=str,
                                 help='gap annotation file', default='')
    parser.add_argument('-e', '--minEdge', dest='minEdge', type=int,
                                 help='min read depth of a valid edge for wtdbg2', default=0)
    parser.add_argument('-n', '--nodeLen', dest='nodeLen', type=int,
                                 help='node length for wtdbg2, times of 256bp', default=256)
    parser.add_argument('-B', '--blacklist', dest='blackListPath', type=str,
                                 help='blacklist BED file', default='')
    parser.add_argument('-C', '--class', dest='classFn', type=str,
                                 help='transposon class file', default='')
    parser.add_argument('-T', '--refTe', dest='teFn', type=str,
                                 help='transposon reference fasta file', default='')
    parser.add_argument('-R', '--refFa', dest='refFn', type=str,
                                 help='reference genome fasta file', default='')
    parser.add_argument('--germ', dest='germModelPath', type=str,
                                 help='path of germline insertion model', default='')
    parser.add_argument('--soma', dest='somaModelPath', type=str,
                                 help='path of somatic insertion model', default='')
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
    parser.add_argument('-O', '--overhang', dest='minOverhang', type=int,
                                 help='min overhang length', default=200)
    
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    cmdArgs = parseArgs()

    if cmdArgs.germModelPath.endswith('/') == False:
        cmdArgs.germModelPath += '/'
    
    if cmdArgs.somaModelPath.endswith('/') == False:
        cmdArgs.somaModelPath += '/'
    
    if os.path.exists('tmp_build') == False:
        os.mkdir('tmp_build')
    
    if os.path.exists('tmp_assm') == False:
        os.mkdir('tmp_assm')
    
    if os.path.exists('tmp_anno') == False:
        os.mkdir('tmp_anno')

    tidToResult = runInParallel(cmdArgs)

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