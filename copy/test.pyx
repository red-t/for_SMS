cdef class testIns:
    cdef public:
        int start, end
        list l
    
    def __init__(self, s, e):
        if s is not None:
            self.start = s
        
        if e is not None:
            self.end = e
        
        self.l = []
    
    cdef add_seg(self, int n):
        self.l.append(n)


cpdef test_p(int s, int e, int n):
    cdef:
        testIns a
    
    a = testIns(s, e)
    a.add_seg(n)

    return a

cpdef test_a():
    cdef bytes a = 'hello world\n'.encode()
    cdef const char* b = a
    return b
        

# from libc.stdio cimport printf, fprintf, FILE, fopen, fclose
# # cdef FILE* fout = fopen("t.out", "w")
# cdef bytes m = test_a()
# cdef const char* n = m
# printf(n)
# # fprintf(fout, b)
# # fclose(fout)

from libc.stdio cimport FILE, fopen, fclose, fprintf, printf
to_p = ['hello world {}\n'.format(x) for x in range(1000)]

cpdef write_file_c(list to_p):
    cdef:
        str a
        bytes b, fout_name = "print_c.out".encode("UTF-8")
        const char* fout_path = fout_name
        const char* c
        FILE* fout = fopen(fout_path, "w")
    
    for a in to_p:
        b = a.encode('UTF-8')
        c = b
        fprintf(fout, c)
    
    fclose(fout)

cpdef write_file_p(list to_p):
    cdef:
        object fout
        str line, fout_name = "print_p.out"
    
    with open(fout_name, "w") as fout:
        for line in to_p:
            fout.write(line)