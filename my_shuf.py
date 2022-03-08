import sys, getopt
import random, linecache
import pickle

######## Help Information ########
def help_info():
    print("my_shuf.py")
    print("\t-p <prefix>\tprefix for the output file.")
    print("\t-M <int>\tinertions number of each subset that use to generate target haplotype.")
    print("\t-R <ref.summary>\ttotal simulated insertions that used for generate subsets.")
    print("\t-H \twhether generate subset of homozygous insertions.")
    print("\t-h \tShow this information")


######## Getting parameters ########
try:
    opts, args = getopt.getopt(argv,"p:M:R:H:h")
except getopt.GetoptError:
    help_info()
    sys.exit(2)

PREFIX=None; N_INS=None; REF_SUMMARY=None; HOMOZYGOUS=False
for opt, arg in opts:
    if opt == '-h':
        help_info()
        sys.exit()
    elif opt == '-p':
        PREFIX = arg
    elif opt == '-M':
        N_INS = int(arg)
    elif opt == '-R':
        REF_SUMMARY = arg
    elif opt == '-H':
        HOMOZYGOUS = True

### Checking ###
print("PREFIX:\t", PREFIX, sep='')
print("N_INS:\t", str(N_INS), sep='')
print("REF_SUMMARY:\t", REF_SUMMARY, sep='')
print("HOMOZYGOUS:\t", HOMOZYGOUS, sep='')
### Checking ###


### GENERATE INSERTIONS SUBSET(S) ###
id = PREFIX.split(".")[0]

if HOMOZYGOUS:
    # GENERATE HOMOZYGOUS SUBSET
    fnrow = len(open(REF_SUMMARY, "r").readlines())
    all_set = set(range(2, fnrow+1))
    homo_set = random.sample(all_set, round(0.5*N_INS))
    with open("{}.homozygous.groundtruth.summary".format(id), "w") as homo_out:
        for tmp_line in homo_set:
            tmpdata = linecache.getline(REF_SUMMARY, tmp_line)
            homo_out.write(tmpdata)
    
    # RANDOMLY GENERATE SUBSET OF INSERTIONS
    all_set = all_set.difference(homo_set)
    rand_set = random.sample(all_set, round(0.5*N_INS))
    with open("{}.groundtruth.summary".format(PREFIX), "w") as rand_out:
        for tmp_line in rand_set:
            tmpdata = linecache.getline(REF_SUMMARY, tmp_line)
            rand_out.write(tmpdata)

    pic_out = open("{}.tmp.pickle".format(id), "wb")
    pickle.dump(all_set, pic_out)
    pic_out.close()

else:
    # RANDOMLY GENERATE SUBSET OF INSERTIONS
    pic_in = open("{}.tmp.pickle".format(id), "rb")
    all_set = pickle.load(pic_in)
    pic_in.close()
    rand_set = random.sample(all_set, round(0.5*N_INS))
    with open("{}.groundtruth.summary".format(PREFIX), "w") as rand_out:
        for tmp_line in rand_set:
            tmpdata = linecache.getline(REF_SUMMARY, tmp_line)
            rand_out.write(tmpdata)