import re

cdef class SequenceContainer:
    """
    Contains all the defined TE sequences
    a.) the ones loaded from the file
    b.) the ones defined in the TE definition file
    (also used for parsing the TE-definition Domain-specific language)
    """
    def __init__(self, startseqlist):
        self.__ssl = startseqlist # sequences loaded from the input file
        self.__sc = {} # sequence container

    def get_count_definitions(self):
        return len(self.__sc.keys())

    cdef addDefinition(self, str definition):
        """
        Add a new TE definition to the container 
        """
        cdef:
            str d, seqid

        if "=" not in definition:
            raise Exception("Invalid defintion; must contain '=' " + definition)
        d, seqid = re.split("=", definition)
        d = d.strip(" ") # remove white spaces
        if d in self.__sc:
            raise Exception("Shortcut for sequence does already exist: " + d)
        self.__sc[d] = self.getTESequence(seqid)

    cpdef TESequence getTESequence(self, str id):
        """
        either direct get sequence by shortcut: e.g. d, e, $1, 1
        or use a definition e.g. $1[100-200]+
        """
        """
        Algorithm for parsing the domain specific language
        a) extract the nested from definition string
        b) extract the deletion-definition from the definition string
        c) based on +- split definition string into left and right half; if no +- use the definition string as left half
        d) get the sequence
        e) introduce deletions
        f) fix the strand (reverse complement sequence if necessary)
        g) introduce mutations and capture tsd 
        f) introduce nested insertions
        """
        ## remove white spaces
        id = id.strip(" ")
        ## direct access; no strand info +-

        # extract NESTED
        cdef str nestdef = ""
        if "{" in id:
            m = re.search(r"^(.*?){(.*)}(.*?)$", id)
            if m is None:
                raise Exception("Invalid TE definition; nested insertion")
            id = m.group(1) + m.group(3)
            nestdef = m.group(2)

        # extract DELETIONS
        cdef str deldef = ""
        if "[" in id:
            m = re.search(r"^(.*?)\[(.*)\](.*?)$", id)
            if m is None:
                raise Exception("Invalid TE definition; deletion")
            id = m.group(1) + m.group(3)
            deldef = m.group(2)
            
        # define left side, right side etc
        cdef:
            str left = id
            str right = ""
            str strand = ""
            TESequence teseq

        # if strand is provided divide into 
        m = re.search(r"([+-])", id)
        if m is not None:
            strand = m.group(1)
            left, right = re.split("[+-]", id)
        left = left.strip(" ") # make the hash access space secure
        
        # get the BASE sequence
        teseq = self.__teseqfromHash(left)

        # make deep copy before modifying the sequence
        teseq = TESequence(teseq.sequence, teseq.id, teseq.tsd)
        
        # process deletions
        teseq = self.__process_deletions(teseq, deldef)

        # fix the strand
        teseq = self.__fixstrand(teseq, strand)

        # process right (mutations and tsd)
        teseq = self.__process_right(teseq, right)

        # process nested insertions (possibly recursion)
        teseq = self.__process_nested(teseq, nestdef)

        return teseq

    cdef TESequence __process_nested(self, TESequence teseq, str nestdef):
        """
        the challenge: this code must be recursion safe (ie no splitting or processing of inner TE definitions, inside of {}
        """
        #### NESTED #####
        if nestdef == "":
            return teseq

        # split at "," in a recursion safe manner (not considering nested)
        cdef:
            list tonest = self.__recursionsafe_split(nestdef)
            list pairs = []
            str t, p, seqid
            int pos
            TESequence seq
        for t in tonest:
            m = re.search(r"^([^:]*):(.*)$", t)  ### RECURSION SAVE; SPLIT only at first occurence; tested: check
            if m is None:
                raise Exception("Invalid definition of nested insertion "+t)
            p = m.group(1)
            seqid = m.group(2)
            pos = self.__getPosition(p, teseq)
            seq = self.getTESequence(seqid) # recursion :)
            pairs.append((pos, seq))

        teseq.sequence = insertSequences(teseq.sequence, pairs)
        return teseq
    
    cdef list __recursionsafe_split(self, str nestdef):
        """
        spliting by comma (,) only allowed at the top level of nested insertions
        2:f{3,5},3:g{h:2,g:5}
        should be: 
        2:f{3,5}
        g{h:2,g:5}
        """
        # TESTED and works
        # FIND POSITIONS
        cdef:
            int nestlev = 0 # nestedlevel
            int dellev = 0
            int i
            str c
            positions = []
        for i, c in enumerate(nestdef):
            if c == "{":
                nestlev += 1
            elif c == "}":
                nestlev -= 1
            elif c == "[":
                dellev += 1
            elif c == "]":
                dellev -= 1
            elif c == ",":
                if nestlev==0 and dellev==0:
                    positions.append(i)
                    
        if nestlev != 0:
            raise Exception("Invalid definition string number of '{' does not match number of '}'")

        positions= reversed(positions)

        # SPLIT
        cdef:
            list toret = []
            str pre = nestdef
            int p
            str post
        for p in positions:
            post = pre[p+1:]
            pre = pre[:p]
            toret.append(post)

        toret.append(pre)
        return toret

    cdef TESequence __process_right(self, TESequence teseq, str right):
        """
        right:
        only TSD and mutations
        """
        if right == "":
            return teseq

        ## Separate nested insertions from the rest % and bp
        # TSD
        m = re.search(r"([\d\.]+)bp", right)
        if m is not None:
            teseq.tsd = int(m.group(1))

        # sequence divergence in percent	
        cdef:
            float div = 0.0
            ExhaustiveSeqMutator mutator
        m = re.search(r"([\d\.]+)%", right)
        if m is not None:
            div = float(m.group(1))
            mutator = ExhaustiveSeqMutator(float(div/100.0))
            teseq.sequence = mutator.mutateseq(teseq.sequence)  ## SET THE NEW SEQUENCE DIVERGENCE ONLY WHEN NEEDED
        return teseq

    cdef TESequence __fixstrand(self, TESequence teseq, str strand):
        if strand == "":
            return teseq
        if strand == "+":
            return teseq
        elif strand == "-":
            teseq.sequence = rc(teseq.sequence)
            return teseq
        else:
            raise Exception("invalid strand " + str(strand))

    cdef TESequence __process_deletions(self, TESequence teseq, str deldef):
        # if no deletions return unchanged
        if deldef == "":
            return teseq

        cdef:
            list deletions = [deldef]
            set delset = set([])
            str d, st, ed
            int start, end, i

        if "," in deldef:
            deletions = re.split(",", deldef)

        for d in deletions:
            st, ed = re.split("\.\.", d)
            start = self.__getPosition(st, teseq)
            end = self.__getPosition(ed, teseq)
            for i in range(start, end+1):
                delset.add(i)
        
        return self.__commence_deleting(teseq, delset)

    cdef TESequence __commence_deleting(self, TESequence teseq, set delset):
        """
        delset = a set of sites that need to be removed
        also handles overlapping deletions
        """
        cdef:
            list seq = list(teseq.sequence)
            list todel = sorted(list(delset), key=lambda i:-i)
            int i
            str newseq
            TESequence toret

        for i in todel:
            del seq[i-1]

        newseq = "".join(seq)
        toret = TESequence(newseq, teseq.id, teseq.tsd)
        return toret
    
    cdef TESequence __teseqfromHash(self, str seqid):
        """
        get the sequence for a sequence id
        - if it starts with $ use the number from the file
        - if not get it from the sequence hash
        """
        cdef:
            int id
            str seq
            TESequence tes
        
        if seqid.startswith("$"):
            # when using the TE from a file, a new TESequence object needs to be generated
            id = int(seqid[1:]) - 1
            seq = self.__ssl[id]
            tes = TESequence(seq, id, 0) # sequence, id,  tsd
            return tes
        else:
            if seqid not in self.__sc:
                raise Exception("Unknown sequence " + seqid)
            return self.__sc[seqid]
    
    cdef int __getPosition(self, str pos, TESequence teseq):
        """
        translate the position provided by the user into an int
        #accomodate ^ $ | beginning, end, middle
        """
        if pos == "^":
            return 1
        elif pos == "$":
            return len(teseq.sequence)
        elif pos == "|":
            return int(len(teseq.sequence)/2.0)
        else:
            return int(pos)