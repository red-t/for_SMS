import random
import re

cpdef set getGermOrder(int numGerm, int numSoma):
    cdef set germOrderSet = set()
    cdef int numTotal = numGerm + numSoma - 1
    cdef int order
    
    while len(germOrderSet) < numGerm:
        order = random.randint(0, numTotal)
        if order in germOrderSet:
            continue
        germOrderSet.add(order)
    
    return germOrderSet
#
#
#
cpdef defineHeader(object args, set germOrderSet):
    # Get template information
    cdef int chromLen
    cdef str chrom
    cdef str templateIndex = args.templateFa + ".fai"

    for row in open(templateIndex, "r"):
        row = row.split()
        chrom = row[0]
        chromLen = int(row[1])
        with open("{0}.tmp.pgd.header".format(chrom), "w") as fout:
            fout.write("# Chasis {0}; Length {1} nt\n".format(chrom, row[1]))
    
    # Parse TE information
    cdef int teOrder = 1
    cdef str teIndex = args.teFa + ".fai"
    cdef dict teDict = {}

    for row in open(teIndex, "r"):
        row = row.split()
        if row[0] != "":
            teDict[teOrder] = (row[0], int(row[1]))
            teOrder += 1

    # Define (args.numGerm + args.numSoma) insertions
    cdef list insIdList = []
    cdef int numIns = args.numGerm + args.numSoma
    cdef int insOrder, numTe = len(teDict)
    cdef str insId, insExpress
    cdef dict insIdToExpressionDict = {}
    fout = open("{0}.tmp.pgd.header".format(chrom), "a")

    for insOrder in range(1, numIns+1):
        if (insOrder-1) in germOrderSet:
            insId, insExpress = defineIns(insOrder, teDict, numTe, "g", args)
        else:
            insId, insExpress = defineIns(insOrder, teDict, numTe, "s", args)
        
        insIdList.append(insId)
        insIdToExpressionDict[insId] = insExpress
        fout.write("{0}={1}\n".format(insId, insExpress))
    
    fout.close()
    return insIdToExpressionDict, chrom, chromLen, insIdList
#
#
#
cpdef tuple defineIns(int insOrder, dict teDict, int numTe, str insType, object args):
    cdef int teOrder, teLen, truncateLen
    cdef str teId, truncateId, truncateExpress

    # Generate data for training model
    if args.mode == 1:
        ### Choose TE ###
        teOrder = getTeOrderTrain(insType, numTe, args.species)
        teId = teDict[teOrder][0]
        teLen = teDict[teOrder][1]

        ### Define Truncations ###
        truncateId, truncateExpress, truncateLen = defineTruncateTrain(args.truncProb, teLen)

    else:
        ### Choose TE ###
        teOrder = getTeOrderTest(args.species)
        teId = teDict[teOrder][0]
        teLen = teDict[teOrder][1]

        ### Define Truncations ###
        if args.species == "human":
            truncateId, truncateExpress, truncateLen = defineTruncateTest(teOrder, teLen)
        else:
            truncateId, truncateExpress, truncateLen = defineTruncateTrain(args.truncProb, teLen)

    ### Choose TSD ###
    cdef int tsdLen = random.randint(6,20)

    ### Choose Strand ###
    cdef str strand, strandId
    strand, strandId = getStrand()

    ### Define Nested Insertion ###
    cdef str nestId, nestExpress
    nestId, nestExpress = defineNestIns(teDict, numTe, teLen, truncateLen, insType, args)
    
    ### Output Insertion Definition ###
    cdef str insId, insExpress
    cdef float insDivRate = args.insDivRate

    insId = "{0}~{1}~{2}{3}{4}".format(insOrder, teId, truncateId, strandId, nestId)
    if insDivRate > 0:
        insExpress = "${0}{1}{2}{3}{4}%{5}bp".format(teOrder, truncateExpress, strand, nestExpress, insDivRate, tsdLen)
    else:
        insExpress = "${0}{1}{2}{3}{4}bp".format(teOrder, truncateExpress, strand, nestExpress, tsdLen)

    return insId, insExpress
#
#
#
cdef int getTeOrderTrain(str insType, int numTe, str species):
    if insType == "g":
        return random.randint(1, numTe)
    if species == "human":
        # Alu、HERVK、LINE1、SVA
        return random.choice([1, 2, 3, 4])
    # 3S18、Max_element、blood、HMS_Beagle、I_element、P_element
    return random.choice([15, 44, 54, 72, 100, 106])
#
#
#
cdef int getTeOrderTest(str species):
    if species == "human":
        # Alu、HERVK、LINE1、SVA
        return random.choices([1, 2, 3, 4], [0.83, 0.01, 0.12, 0.04], k=1)[0]
    # 3S18、Max_element、blood、HMS_Beagle、I_element、P_element
    return random.choices([15, 44, 54, 72, 100, 106], [0.06, 0.14, 0.07, 0.07, 0.08, 0.58], k=1)[0]
#
#
#
cdef tuple getStrand():
    if random.random() < 0.5:
        return "-", "R"

    return "+", "F"
#
#
#
cdef tuple defineTruncateTrain(float truncProb, int teLen):
    cdef str truncateId = "U", truncateExpress = ""
    cdef int truncateStart, truncateLen = 0

    if random.random() > truncProb:
        return truncateId, truncateExpress, truncateLen
    
    truncateId = "T"
    truncateLen = int(random.uniform(0.1, 0.3) * teLen)
    truncateStart = random.randint(1, teLen - truncateLen)
    truncateExpress = "[{0}..{1}]".format(truncateStart, truncateStart + truncateLen - 1)
    return truncateId, truncateExpress, truncateLen
#
#
#
cdef tuple defineTruncateTest(int teOrder, int teLen):
    cdef str truncateId = "U", truncateExpress = ""
    cdef int truncateStart, truncateLen = 0

    # Only LINE1
    if teOrder != 3:
        return truncateId, truncateExpress, truncateLen
    
    # Full length = 30%
    cdef float truncProb = random.random()
    if truncProb <= 0.3:
        return truncateId, truncateExpress, truncateLen
    
    truncateId = "T"
    truncateLen = int(random.uniform(0.1, 0.3) * teLen)
    
    # Internal truncation = 20%
    if truncProb <= 0.5:
        truncateStart = random.randint(1, teLen - truncateLen)
        truncateExpress = "[{0}..{1}]".format(truncateStart, truncateStart + truncateLen - 1)
        return truncateId, truncateExpress, truncateLen
    
    # 5' truncation = 50%
    truncateStart = 1
    truncateExpress = "[{0}..{1}]".format(truncateStart, truncateStart + truncateLen - 1)
    return truncateId, truncateExpress, truncateLen
#
#
#
cpdef defineNestIns(dict teDict, int numTe, int parentTeLen, int parentTruncateLen, str insType, object args):
    ### No Nesting ###
    cdef str insId = "", insExpress = ""
    cdef float nestProb = random.random()
    if nestProb > args.nestProb:
        return insId, insExpress

    ### Choose TE ###
    cdef int teOrder, teLen
    cdef str teId
    teOrder = getTeOrderTrain(insType, numTe, args.species)
    teId = teDict[teOrder][0]
    teLen = teDict[teOrder][1]

    ### Define Truncations ###
    cdef int truncateLen
    cdef str truncateId, truncateExpress
    truncateId, truncateExpress, truncateLen = defineTruncateTrain(args.truncProb, teLen)

    ### Choose TSD ###
    cdef int tsdLen = random.randint(6, 20)
    
    ### Choose Strand ###
    cdef str strand, strandId
    strand, strandId = getStrand()
    
    ### Define Nested Insertion ###
    cdef str nestId, nestExpress
    nestId, nestExpress = defineNestIns(teDict, numTe, teLen, truncateLen, insType, args)
    
    ### Output Insertion Definition ###
    cdef int intPos = random.randint(tsdLen, parentTeLen - parentTruncateLen)
    cdef float insDivRate = args.insDivRate

    insId = "({0}~{1}{2}{3})".format(teId, truncateId, strandId, nestId)
    if insDivRate > 0:
        insExpress = "{" + "{0}:${1}{2}{3}{4}{5}%{6}bp".format(intPos, teOrder, truncateExpress, strand, nestExpress, insDivRate, tsdLen) + "}"
    else:
        insExpress = "{" + "{0}:${1}{2}{3}{4}{5}bp".format(intPos, teOrder, truncateExpress, strand, nestExpress, tsdLen) + "}"
    return insId, insExpress
#
#
#
cpdef defineBody(dict insIdToExpressDict, str chrom, int minDist, int maxDist, set germOrderSet, list insIdList, object args):
    cdef int insOrder, genomeOrder
    cdef int minPos, maxPos, insPos
    cdef int numGenome, numEmpty
    cdef int tsdStart, tsdLen
    cdef int numTotal = (args.numGerm + args.numSoma)
    cdef float frequency
    cdef str insId, body, insExpress, strand
    cdef list bodyList

    fout = open("{0}.tmp.pgd.body".format(chrom), "w")
    summary = open("{0}.ins.summary".format(chrom), "w")

    for insOrder in range(numTotal):
        ### Set insertion Position ###
        if insOrder > 0:
            minPos = maxDist * insOrder
            maxPos = maxDist * (insOrder+1)
            if minPos - insPos < minDist:
                minPos = insPos + minDist
            insPos = random.randint(minPos, maxPos)
        else:
            insPos = random.randint(1, maxDist)

        ### Set Insertion Frequency ###
        if insOrder in germOrderSet:
            if args.species == "human":
                frequency = random.random()
                if frequency < 0.1:
                    frequency = random.uniform(0.1, 1)
                else:
                    frequency = random.choice((0.5, 1))
            if args.species == "fly":
                frequency = random.uniform(0.1, 1)
            numGenome = int(frequency * args.numTotalGenome)
        else:
            frequency = float(1) / args.numTotalGenome
            numGenome = 1
        
        numEmpty = args.numTotalGenome - numGenome

        ### Construct bodyList ###
        insId = insIdList[insOrder]
        bodyList = [insId for genomeOrder in range(numGenome)]
        bodyList.extend("*" * numEmpty)
        random.shuffle(bodyList)
        bodyList.insert(0, str(insPos))

        ### Write Out Body ###
        body = " ".join(bodyList)
        fout.write(body + "\n")

        ### Write Out Insertion Information ###
        insExpress = insIdToExpressDict[insId]
        strand = re.search(r"([+-])", insExpress).group(1)
        tsdLen = int(re.search("([\d\.]+)bp$", insExpress).group(1))
        tsdStart = insPos - tsdLen + 1
        summary.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format(chrom, insPos, insPos+1, insId, insExpress, strand, frequency, tsdStart, insPos))
    
    fout.close()
    summary.close()