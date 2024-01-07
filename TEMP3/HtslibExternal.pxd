from libc.stdint cimport int8_t, int16_t, int32_t, int64_t
from libc.stdint cimport uint8_t, uint16_t, uint32_t, uint64_t
from libc.stdlib cimport malloc, calloc, realloc, free
from libc.string cimport memcpy, memcmp, strncpy, strlen, strdup
from libc.stdio cimport FILE, printf
from cpython cimport PyBytes_Check
from posix.types cimport off_t


########################
### htslib/kstring.h ###
########################
cdef extern from "htslib/kstring.h" nogil:
    ctypedef struct kstring_t:
        size_t l, m
        char *s
#     int kputc(int c, kstring_t *s)
#     int kputw(int c, kstring_t *s)
#     int kputl(long c, kstring_t *s)
#     int ksprintf(kstring_t *s, const char *fmt, ...)

######################
### htslib/hfile.h ###
######################
cdef extern from "htslib/hfile.h" nogil:
    ctypedef struct hFILE

#     # @abstract  Open the named file or URL as a stream
#     # @return    An hFILE pointer, or NULL (with errno set) if an error occurred.
#     hFILE *hopen(const char *filename, const char *mode, ...)

#     # @abstract  Associate a stream with an existing open file descriptor
#     # @return    An hFILE pointer, or NULL (with errno set) if an error occurred.
#     # @notes     For socket descriptors (on Windows), mode should contain 's'.
#     hFILE *hdopen(int fd, const char *mode)

#     # @abstract  Report whether the file name or URL denotes remote storage
#     # @return    0 if local, 1 if remote.
#     # @notes     "Remote" means involving e.g. explicit network access, with the
#     #   implication that callers may wish to cache such files' contents locally.
#     int hisremote(const char *filename)

#     # @abstract  Flush (for output streams) and close the stream
#     # @return    0 if successful, or EOF (with errno set) if an error occurred.
#     int hclose(hFILE *fp)

#     # @abstract  Close the stream, without flushing or propagating errors
#     # @notes     For use while cleaning up after an error only.  Preserves errno.
#     void hclose_abruptly(hFILE *fp)

#     # @abstract  Return the stream's error indicator
#     # @return    Non-zero (in fact, an errno value) if an error has occurred.
#     # @notes     This would be called herror() and return true/false to parallel
#     #   ferror(3), but a networking-related herror(3) function already exists.  */
#     int herrno(hFILE *fp)

#     # @abstract  Clear the stream's error indicator
#     void hclearerr(hFILE *fp)

#     # @abstract  Reposition the read/write stream offset
#     # @return    The resulting offset within the stream (as per lseek(2)),
#     #   or negative if an error occurred.
#     off_t hseek(hFILE *fp, off_t offset, int whence)

#     # @abstract  Report the current stream offset
#     # @return    The offset within the stream, starting from zero.
#     off_t htell(hFILE *fp)

#     # @abstract  Read one character from the stream
#     # @return    The character read, or EOF on end-of-file or error
#     int hgetc(hFILE *fp)

#     # Read from the stream until the delimiter, up to a maximum length
#     #    @param buffer  The buffer into which bytes will be written
#     #    @param size    The size of the buffer
#     #    @param delim   The delimiter (interpreted as an `unsigned char`)
#     #    @param fp      The file stream
#     #    @return  The number of bytes read, or negative on error.
#     #    @since   1.4
#     #
#     # Bytes will be read into the buffer up to and including a delimiter, until
#     # EOF is reached, or _size-1_ bytes have been written, whichever comes first.
#     # The string will then be terminated with a NUL byte (`\0`).
#     ssize_t hgetdelim(char *buffer, size_t size, int delim, hFILE *fp)

#     # Read a line from the stream, up to a maximum length
#     #    @param buffer  The buffer into which bytes will be written
#     #    @param size    The size of the buffer
#     #    @param fp      The file stream
#     #    @return  The number of bytes read, or negative on error.
#     #    @since   1.4
#     #
#     # Specialization of hgetdelim() for a `\n` delimiter.
#     ssize_t hgetln(char *buffer, size_t size, hFILE *fp)

#     # Read a line from the stream, up to a maximum length
#     #    @param buffer  The buffer into which bytes will be written
#     #    @param size    The size of the buffer (must be > 1 to be useful)
#     #    @param fp      The file stream
#     #    @return  _buffer_ on success, or `NULL` if an error occurred.
#     #    @since   1.4
#     #
#     # This function can be used as a replacement for `fgets(3)`, or together with
#     # kstring's `kgetline()` to read arbitrarily-long lines into a _kstring_t_.
#     char *hgets(char *buffer, int size, hFILE *fp)

#     # @abstract  Peek at characters to be read without removing them from buffers
#     # @param fp      The file stream
#     # @param buffer  The buffer to which the peeked bytes will be written
#     # @param nbytes  The number of bytes to peek at; limited by the size of the
#     #   internal buffer, which could be as small as 4K.
#     # @return    The number of bytes peeked, which may be less than nbytes if EOF
#     #   is encountered; or negative, if there was an I/O error.
#     # @notes  The characters peeked at remain in the stream's internal buffer,
#     #   and will be returned by later hread() etc calls.
#     ssize_t hpeek(hFILE *fp, void *buffer, size_t nbytes)

#     # @abstract  Read a block of characters from the file
#     # @return    The number of bytes read, or negative if an error occurred.
#     # @notes     The full nbytes requested will be returned, except as limited
#     #   by EOF or I/O errors.
#     ssize_t hread(hFILE *fp, void *buffer, size_t nbytes)

#     # @abstract  Write a character to the stream
#     # @return    The character written, or EOF if an error occurred.
#     int hputc(int c, hFILE *fp)

#     # @abstract  Write a string to the stream
#     # @return    0 if successful, or EOF if an error occurred.
#     int hputs(const char *text, hFILE *fp)

#     # @abstract  Write a block of characters to the file
#     # @return    Either nbytes, or negative if an error occurred.
#     # @notes     In the absence of I/O errors, the full nbytes will be written.
#     ssize_t hwrite(hFILE *fp, const void *buffer, size_t nbytes)

#     # @abstract  For writing streams, flush buffered output to the underlying stream
#     # @return    0 if successful, or EOF if an error occurred.
#     int hflush(hFILE *fp)

#####################
### htslib/bgzf.h ###
#####################
cdef extern from "htslib/bgzf.h" nogil:
    ctypedef struct bgzf_mtaux_t
    ctypedef struct bgzidx_t
    ctypedef struct z_stream

    ctypedef struct BGZF:
        unsigned           errcode
        unsigned           is_write
        int           is_be
        int           compress_level
        int           is_compressed
        int           is_gzip
        int           cache_size
        int64_t       block_address
        int64_t       uncompressed_address
        void         *uncompressed_block
        void         *compressed_block
        void         *cache
        hFILE        *fp
        bgzf_mtaux_t *mt
        bgzidx_t     *idx
        int           idx_build_otf
        z_stream     *gz_stream

    # #*****************
    # #  Basic routines *
    # # *****************/

    # #  Open an existing file descriptor for reading or writing.
    # #
    # #  @param fd    file descriptor
    # #  @param mode  mode matching /[rwag][u0-9]+/: 'r' for reading, 'w' for
    # #               writing, 'a' for appending, 'g' for gzip rather than BGZF
    # #               compression (with 'w' only), and digit specifies the zlib
    # #               compression level.
    # #               Note that there is a distinction between 'u' and '0': the
    # #               first yields plain uncompressed output whereas the latter
    # #               outputs uncompressed data wrapped in the zlib format.
    # #  @return      BGZF file handler; 0 on error

    # BGZF* bgzf_dopen(int fd, const char *mode)
    # BGZF* bgzf_fdopen(int fd, const char *mode) # for backward compatibility

    # #  Open the specified file for reading or writing.
    # BGZF* bgzf_open(const char* path, const char *mode)

    # #  Open an existing hFILE stream for reading or writing.
    # BGZF* bgzf_hopen(hFILE *fp, const char *mode)

    # #  Close the BGZF and free all associated resources.
    # #
    # #  @param fp    BGZF file handler
    # #  @return      0 on success and -1 on error
    # int bgzf_close(BGZF *fp)

    # #  Read up to _length_ bytes from the file storing into _data_.
    # #
    # #  @param fp     BGZF file handler
    # #  @param data   data array to read into
    # #  @param length size of data to read
    # #  @return       number of bytes actually read; 0 on end-of-file and -1 on error
    # ssize_t bgzf_read(BGZF *fp, void *data, size_t length)

    # #  Write _length_ bytes from _data_ to the file.  If no I/O errors occur,
    # #  the complete _length_ bytes will be written (or queued for writing).
    # #
    # #  @param fp     BGZF file handler
    # #  @param data   data array to write
    # #  @param length size of data to write
    # #  @return       number of bytes written (i.e., _length_); negative on error
    # ssize_t bgzf_write(BGZF *fp, const void *data, size_t length)

    # #  Read up to _length_ bytes directly from the underlying stream without
    # #  decompressing.  Bypasses BGZF blocking, so must be used with care in
    # #  specialised circumstances only.
    # #
    # #  @param fp     BGZF file handler
    # #  @param data   data array to read into
    # #  @param length number of raw bytes to read
    # #  @return       number of bytes actually read; 0 on end-of-file and -1 on error
    # ssize_t bgzf_raw_read(BGZF *fp, void *data, size_t length)

    # #  Write _length_ bytes directly to the underlying stream without
    # #  compressing.  Bypasses BGZF blocking, so must be used with care
    # #  in specialised circumstances only.
    # #
    # #  @param fp     BGZF file handler
    # #  @param data   data array to write
    # #  @param length number of raw bytes to write
    # #  @return       number of bytes actually written; -1 on error
    # ssize_t bgzf_raw_write(BGZF *fp, const void *data, size_t length)

    # #  Write the data in the buffer to the file.
    # int bgzf_flush(BGZF *fp)

    int SEEK_SET

    #  Return a virtual file pointer to the current location in the file.
    #  No interpretation of the value should be made, other than a subsequent
    #  call to bgzf_seek can be used to position the file at the same point.
    #  Return value is non-negative on success.
    int64_t bgzf_tell(BGZF *fp)

    #  Set the file to read from the location specified by _pos_.
    #
    #  @param fp     BGZF file handler
    #  @param pos    virtual file offset returned by bgzf_tell()
    #  @param whence must be SEEK_SET
    #  @return       0 on success and -1 on error
    # /
    int64_t bgzf_seek(BGZF *fp, int64_t pos, int whence)

    # #  Check if the BGZF end-of-file (EOF) marker is present
    # #
    # #  @param fp    BGZF file handler opened for reading
    # #  @return      1 if the EOF marker is present and correct
    # #               2 if it can't be checked, e.g., because fp isn't seekable
    # #               0 if the EOF marker is absent
    # #               -1 (with errno set) on error
    # int bgzf_check_EOF(BGZF *fp)

    # #  Check if a file is in the BGZF format
    # #
    # #  @param fn    file name
    # #  @return      1 if _fn_ is BGZF; 0 if not or on I/O error
    # int bgzf_is_bgzf(const char *fn)

    # #*********************
    # #  Advanced routines *
    # #*********************

    # #  Set the cache size. Only effective when compiled with -DBGZF_CACHE.
    # #
    # #  @param fp    BGZF file handler
    # #  @param size  size of cache in bytes; 0 to disable caching (default)
    # void bgzf_set_cache_size(BGZF *fp, int size)

    # #  Flush the file if the remaining buffer size is smaller than _size_
    # #  @return      0 if flushing succeeded or was not needed; negative on error
    # int bgzf_flush_try(BGZF *fp, ssize_t size)

    # #  Read one byte from a BGZF file. It is faster than bgzf_read()
    # #  @param fp     BGZF file handler
    # #  @return       byte read; -1 on end-of-file or error
    # int bgzf_getc(BGZF *fp)

    # #  Read one line from a BGZF file. It is faster than bgzf_getc()
    # #
    # #  @param fp     BGZF file handler
    # #  @param delim  delimiter
    # #  @param str    string to write to; must be initialized
    # #  @return       length of the string; 0 on end-of-file; negative on error
    # int bgzf_getline(BGZF *fp, int delim, kstring_t *str)

    # #  Read the next BGZF block.
    # int bgzf_read_block(BGZF *fp)

    # #  Enable multi-threading (only effective on writing and when the
    # #  library was compiled with -DBGZF_MT)
    # #
    # #  @param fp          BGZF file handler; must be opened for writing
    # #  @param n_threads   #threads used for writing
    # #  @param n_sub_blks  #blocks processed by each thread; a value 64-256 is recommended
    # int bgzf_mt(BGZF *fp, int n_threads, int n_sub_blks)


    # # Compress a single BGZF block.
    # #
    # # @param dst    output buffer (must have size >= BGZF_MAX_BLOCK_SIZE)
    # # @param dlen   size of output buffer; updated on return to the number
    # #               of bytes actually written to dst
    # # @param src    buffer to be compressed
    # # @param slen   size of data to compress (must be <= BGZF_BLOCK_SIZE)
    # # @param level  compression level
    # # @return       0 on success and negative on error
    # #
    # int bgzf_compress(void *dst, size_t *dlen, const void *src, size_t slen, int level)

    # #*******************
    # #  bgzidx routines *
    # #   BGZF at the uncompressed offset
    # #
    # #   @param fp           BGZF file handler; must be opened for reading
    # #   @param uoffset      file offset in the uncompressed data
    # #   @param where        SEEK_SET supported atm
    # #
    # #   Returns 0 on success and -1 on error.
    # int bgzf_useek(BGZF *fp, long uoffset, int where)

    # #   Position in uncompressed BGZF
    # #
    # #   @param fp           BGZF file handler; must be opened for reading
    # #
    # #   Returns the current offset on success and -1 on error.
    # long bgzf_utell(BGZF *fp)

    # #  Tell BGZF to build index while compressing.
    # #
    # #  @param fp          BGZF file handler; can be opened for reading or writing.
    # #
    # #  Returns 0 on success and -1 on error.
    # int bgzf_index_build_init(BGZF *fp)

    # #  Load BGZF index
    # #
    # #  @param fp          BGZF file handler
    # #  @param bname       base name
    # #  @param suffix      suffix to add to bname (can be NULL)
    # #
    # #  Returns 0 on success and -1 on error.
    # int bgzf_index_load(BGZF *fp, const char *bname, const char *suffix)

    # #  Save BGZF index
    # #
    # #  @param fp          BGZF file handler
    # #  @param bname       base name
    # #  @param suffix      suffix to add to bname (can be NULL)
    # #
    # #  Returns 0 on success and -1 on error.
    # int bgzf_index_dump(BGZF *fp, const char *bname, const char *suffix)

######################
### htslib/hts.h ###
######################
cdef extern from "htslib/hts.h" nogil:
    # uint32_t kroundup32(uint32_t x)

    ctypedef struct cram_fd

    union FilePointerUnion:
        BGZF    *bgzf
        cram_fd *cram
        hFILE   *hfile
        void    *voidp

    enum htsFormatCategory:
        unknown_category
        sequence_data    # Sequence data -- SAM, BAM, CRAM, etc
        variant_data     # Variant calling data -- VCF, BCF, etc
        index_file       # Index file associated with some data file
        region_list      # Coordinate intervals or regions -- BED, etc
        category_maximum

    enum htsExactFormat:
        unknown_format
        binary_format
        text_format
        sam, bam, bai, cram, crai, vcf, bcf, csi, gzi, tbi, bed
        format_maximum

    enum htsCompression:
        no_compression, gzip, bgzf, custom
        compression_maximum

    # cdef enum hts_fmt_option:
    #     CRAM_OPT_DECODE_MD,
    #     CRAM_OPT_PREFIX,
    #     CRAM_OPT_VERBOSITY,
    #     CRAM_OPT_SEQS_PER_SLICE,
    #     CRAM_OPT_SLICES_PER_CONTAINER,
    #     CRAM_OPT_RANGE,
    #     CRAM_OPT_VERSION,
    #     CRAM_OPT_EMBED_REF,
    #     CRAM_OPT_IGNORE_MD5,
    #     CRAM_OPT_REFERENCE,
    #     CRAM_OPT_MULTI_SEQ_PER_SLICE,
    #     CRAM_OPT_NO_REF,
    #     CRAM_OPT_USE_BZIP2,
    #     CRAM_OPT_SHARED_REF,
    #     CRAM_OPT_NTHREADS,
    #     CRAM_OPT_THREAD_POOL,
    #     CRAM_OPT_USE_LZMA,
    #     CRAM_OPT_USE_RANS,
    #     CRAM_OPT_REQUIRED_FIELDS,
    #     HTS_OPT_COMPRESSION_LEVEL,
    #     HTS_OPT_NTHREADS,

    ctypedef struct htsVersion:
        short major, minor

    ctypedef struct htsFormat:
        htsFormatCategory category
        htsExactFormat    format
        htsVersion        version
        htsCompression    compression
        short             compression_level
        void              *specific  

    ctypedef struct htsFile:
        uint8_t  is_bin
        uint8_t  is_write
        uint8_t  is_be
        uint8_t  is_cram
        int64_t lineno
        kstring_t line
        char *fn
        char *fn_aux
        FilePointerUnion fp
        htsFormat format

    # int hts_verbose

    # cdef union hts_opt_val_union:
    #     int i
    #     char *s

    # ctypedef struct hts_opt:
    #     char *arg
    #     hts_fmt_option opt
    #     hts_opt_val_union val
    #     void *next

    # # @abstract Parses arg and appends it to the option list.
    # # @return   0 on success and -1 on failure
    # int hts_opt_add(hts_opt **opts, const char *c_arg)

    # # @abstract Applies an hts_opt option list to a given htsFile.
    # # @return   0 on success and -1 on failure
    # int hts_opt_apply(htsFile *fp, hts_opt *opts)

    # # @abstract Frees an hts_opt list.
    # void hts_opt_free(hts_opt *opts)

    # int hts_parse_format(htsFormat *opt, const char *str)

    # # @abstract Table for converting a nucleotide character to 4-bit encoding.
    # # The input character may be either an IUPAC ambiguity code, '=' for 0, or
    # # '0'/'1'/'2'/'3' for a result of 1/2/4/8.  The result is encoded as 1/2/4/8
    # # for A/C/G/T or combinations of these bits for ambiguous bases.
    # const unsigned char *seq_nt16_table

    # @abstract Table for converting a 4-bit encoded nucleotide to an IUPAC
    # ambiguity code letter (or '=' when given 0).
    const char *seq_nt16_str

    # # @abstract Table for converting a 4-bit encoded nucleotide to about 2 bits.
    # # Returns 0/1/2/3 for 1/2/4/8 (i.e., A/C/G/T), or 4 otherwise (0 or ambiguous).
    # const int *seq_nt16_int

    # # @abstract  Get the htslib version number
    # # @return    For released versions, a string like "N.N[.N]"; or git describe
    # # output if using a library built within a Git repository.
    # const char *hts_version()

    # # @abstract    Determine format by peeking at the start of a file
    # # @param fp    File opened for reading, positioned at the beginning
    # # @param fmt   Format structure that will be filled out on return
    # # @return      0 for success, or negative if an error occurred.
    # int hts_detect_format(hFILE *fp, htsFormat *fmt)

    # # @abstract    Get a human-readable description of the file format
    # # @return      Description string, to be freed by the caller after use.
    # char *hts_format_description(const htsFormat *format)

    # # @abstract       Open a SAM/BAM/CRAM/VCF/BCF/etc file
    # # @param fn       The file name or "-" for stdin/stdout
    # # @param mode     Mode matching / [rwa][bceguxz0-9]* /
    # # @discussion
    # #     With 'r' opens for reading; any further format mode letters are ignored
    # #     as the format is detected by checking the first few bytes or BGZF blocks
    # #     of the file.  With 'w' or 'a' opens for writing or appending, with format
    # #     specifier letters:
    # #       b  binary format (BAM, BCF, etc) rather than text (SAM, VCF, etc)
    # #       c  CRAM format
    # #       g  gzip compressed
    # #       u  uncompressed
    # #       z  bgzf compressed
    # #       [0-9]  zlib compression level
    # #     and with non-format option letters (for any of 'r'/'w'/'a'):
    # #       e  close the file on exec(2) (opens with O_CLOEXEC, where supported)
    # #       x  create the file exclusively (opens with O_EXCL, where supported)
    # #     Note that there is a distinction between 'u' and '0': the first yields
    # #     plain uncompressed output whereas the latter outputs uncompressed data
    # #     wrapped in the zlib format.
    # # @example
    # #     [rw]b  .. compressed BCF, BAM, FAI
    # #     [rw]bu .. uncompressed BCF
    # #     [rw]z  .. compressed VCF
    # #     [rw]   .. uncompressed VCF
    # htsFile *hts_open(const char *fn, const char *mode)

    # # @abstract       Open a SAM/BAM/CRAM/VCF/BCF/etc file
    # # @param fn       The file name or "-" for stdin/stdout
    # # @param mode     Open mode, as per hts_open()
    # # @param fmt      Optional format specific parameters
    # # @discussion
    # #     See hts_open() for description of fn and mode.
    # #     // TODO Update documentation for s/opts/fmt/
    # #     Opts contains a format string (sam, bam, cram, vcf, bcf) which will,
    # #     if defined, override mode.  Opts also contains a linked list of hts_opt
    # #     structures to apply to the open file handle.  These can contain things
    # #     like pointers to the reference or information on compression levels,
    # #     block sizes, etc.
    # htsFile *hts_open_format(const char *fn, const char *mode, const htsFormat *fmt)

    # # @abstract       Open an existing stream as a SAM/BAM/CRAM/VCF/BCF/etc file
    # # @param fp       The already-open file handle
    # # @param fn       The file name or "-" for stdin/stdout
    # # @param mode     Open mode, as per hts_open()
    # htsFile *hts_hopen(hFILE *fp, const char *fn, const char *mode)

    # # @abstract  Close a file handle, flushing buffered data for output streams
    # # @param fp  The file handle to be closed
    # # @return    0 for success, or negative if an error occurred.
    # int hts_close(htsFile *fp)

    # # @abstract  Returns the file's format information
    # # @param fp  The file handle
    # # @return    Read-only pointer to the file's htsFormat.
    # const htsFormat *hts_get_format(htsFile *fp)

    # # @ abstract      Returns a string containing the file format extension.
    # # @ param format  Format structure containing the file type.
    # # @ return        A string ("sam", "bam", etc) or "?" for unknown formats.
    # const char *hts_format_file_extension(const htsFormat *format)

    # # @abstract  Sets a specified CRAM option on the open file handle.
    # # @param fp  The file handle open the open file.
    # # @param opt The CRAM_OPT_* option.
    # # @param ... Optional arguments, dependent on the option used.
    # # @return    0 for success, or negative if an error occurred.
    # int hts_set_opt(htsFile *fp, hts_fmt_option opt, ...)

    # int hts_getline(htsFile *fp, int delimiter, kstring_t *str)
    # char **hts_readlines(const char *fn, int *_n)

    # #   @abstract       Parse comma-separated list or read list from a file
    # #   @param list     File name or comma-separated list
    # #   @param is_file
    # #   @param _n       Size of the output array (number of items read)
    # #   @return         NULL on failure or pointer to newly allocated array of
    # #                   strings
    # char **hts_readlist(const char *fn, int is_file, int *_n)

    # @abstract  Create extra threads to aid compress/decompression for this file
    # @param fp  The file handle
    # @param n   The number of worker threads to create
    # @return    0 for success, or negative if an error occurred.
    # @notes     THIS THREADING API IS LIKELY TO CHANGE IN FUTURE.
    int hts_set_threads(htsFile *fp, int n)

    # # @abstract  Set .fai filename for a file opened for reading
    # # @return    0 for success, negative on failure
    # # @discussion
    # #     Called before *_hdr_read(), this provides the name of a .fai file
    # #     used to provide a reference list if the htsFile contains no @SQ headers.
    # int hts_set_fai_filename(htsFile *fp, const char *fn_aux)

    # int8_t HTS_IDX_NOCOOR
    # int8_t HTS_IDX_START
    # int8_t HTS_IDX_REST
    # int8_t HTS_IDX_NONE

    # int8_t HTS_FMT_CSI
    # int8_t HTS_FMT_BAI
    # int8_t HTS_FMT_TBI
    # int8_t HTS_FMT_CRAI

    # BGZF *hts_get_bgzfp(htsFile *fp)

    ctypedef struct hts_idx_t

    # ctypedef struct hts_pair64_t:
    #     uint64_t u, v

    # ctypedef int hts_readrec_func(BGZF *fp, void *data, void *r, int *tid, int *beg, int *end)

    # ctypedef struct hts_bins_t:
    #     int n, m
    #     int *a

    ctypedef struct hts_pair64_max_t:
        uint64_t u, v
        uint64_t max

    ctypedef struct hts_itr_t:
        int n_off
        hts_pair64_max_t *off
        uint64_t curr_off

    # hts_idx_t *hts_idx_init(int n, int fmt, uint64_t offset0, int min_shift, int n_lvls)
    void hts_idx_destroy(hts_idx_t *idx)
    # int hts_idx_push(hts_idx_t *idx, int tid, int beg, int end, uint64_t offset, int is_mapped)
    # void hts_idx_finish(hts_idx_t *idx, uint64_t final_offset)

    # #### Save an index to a file
    # #    @param idx  Index to be written
    # #    @param fn   Input BAM/BCF/etc filename, to which .bai/.csi/etc will be added
    # #    @param fmt  One of the HTS_FMT_* index formats
    # #    @return  0 if successful, or negative if an error occurred.
    # int hts_idx_save(const hts_idx_t *idx, const char *fn, int fmt)

    # #### Save an index to a specific file
    # #    @param idx    Index to be written
    # #    @param fn     Input BAM/BCF/etc filename
    # #    @param fnidx  Output filename, or NULL to add .bai/.csi/etc to @a fn
    # #    @param fmt    One of the HTS_FMT_* index formats
    # #    @return  0 if successful, or negative if an error occurred.
    # int hts_idx_save_as(const hts_idx_t *idx, const char *fn, const char *fnidx, int fmt)

    # #### Load an index file
    # #    @param fn   BAM/BCF/etc filename, to which .bai/.csi/etc will be added or
    # #                the extension substituted, to search for an existing index file
    # #    @param fmt  One of the HTS_FMT_* index formats
    # #    @return  The index, or NULL if an error occurred.
    # hts_idx_t *hts_idx_load(const char *fn, int fmt)

    # #### Load a specific index file
    # #    @param fn     Input BAM/BCF/etc filename
    # #    @param fnidx  The input index filename
    # #    @return  The index, or NULL if an error occurred.
    # hts_idx_t *hts_idx_load2(const char *fn, const char *fnidx)

    # uint8_t *hts_idx_get_meta(hts_idx_t *idx, uint32_t *l_meta)
    # void hts_idx_set_meta(hts_idx_t *idx, int l_meta, uint8_t *meta, int is_copy)

    # int hts_idx_get_stat(const hts_idx_t* idx, int tid,
    #                      uint64_t* mapped, uint64_t* unmapped)

    # uint64_t hts_idx_get_n_no_coor(const hts_idx_t* idx)

    # int HTS_PARSE_THOUSANDS_SEP  # Ignore ',' separators within numbers

    # # Parse a numeric string
    # #    The number may be expressed in scientific notation, and optionally may
    # #    contain commas in the integer part (before any decimal point or E notation).
    # #    @param str     String to be parsed
    # #    @param strend  If non-NULL, set on return to point to the first character
    # #                   in @a str after those forming the parsed number
    # #    @param flags   Or'ed-together combination of HTS_PARSE_* flags
    # #    @return  Converted value of the parsed number.
    # #
    # #    When @a strend is NULL, a warning will be printed (if hts_verbose is 2
    # #    or more) if there are any trailing characters after the number.
    # long long hts_parse_decimal(const char *str, char **strend, int flags)

    # # Parse a "CHR:START-END"-style region string
    # #    @param str  String to be parsed
    # #    @param beg  Set on return to the 0-based start of the region
    # #    @param end  Set on return to the 1-based end of the region
    # #    @return  Pointer to the colon or '\0' after the reference sequence name,
    # #             or NULL if @a str could not be parsed.
    # const char *hts_parse_reg(const char *str, int *beg, int *end)

    # hts_itr_t *hts_itr_query(const hts_idx_t *idx, int tid, int beg, int end, hts_readrec_func *readrec)
    # void hts_itr_destroy(hts_itr_t *iter)

    # ctypedef int (*hts_name2id_f)(void*, const char*)
    # ctypedef const char *(*hts_id2name_f)(void*, int)
    # ctypedef hts_itr_t *hts_itr_query_func(
    #     const hts_idx_t *idx,
    #     int tid,
    #     int beg,
    #     int end,
    #     hts_readrec_func *readrec)

    # hts_itr_t *hts_itr_querys(
    #     const hts_idx_t *idx,
    #     const char *reg,
    #     hts_name2id_f getid,
    #     void *hdr,
    #     hts_itr_query_func *itr_query,
    #     hts_readrec_func *readrec)

    int hts_itr_next(BGZF *fp, hts_itr_t *iter, void *r, void *data)
    # const char **hts_idx_seqnames(const hts_idx_t *idx, int *n, hts_id2name_f getid, void *hdr)  # free only the array, not the values

    # # hts_file_type() - Convenience function to determine file type
    # # @fname: the file name
    # #
    # # Returns one of the FT_* defines.
    # #
    # # DEPRECATED:  This function has been replaced by hts_detect_format().
    # # It and these FT_* macros will be removed in a future HTSlib release.
    # int FT_UNKN
    # int FT_GZ
    # int FT_VCF
    # int FT_VCF_GZ
    # int FT_BCF
    # int FT_BCF_GZ
    # int FT_STDIN

    # int hts_file_type(const char *fname)

    # # /***************************
    # #  * Revised MAQ error model *
    # #  ***************************/

    # ctypedef struct errmod_t

    # errmod_t *errmod_init(double depcorr)
    # void errmod_destroy(errmod_t *em)

    # # /*
    # #     n: number of bases
    # #     m: maximum base
    # #     bases[i]: qual:6, strand:1, base:4
    # #     q[i*m+j]: phred-scaled likelihood of (i,j)
    # #  */
    # int errmod_cal(const errmod_t *em, int n, int m, uint16_t *bases, float *Probabilistic)

    # # /*****************************************
    # #  * q banded glocal alignment *
    # #  *****************************************/

    # ctypedef struct probaln_par_t:
    #     float d, e
    #     int bw

    # int probaln_glocal(const uint8_t *ref,
    #                    int l_ref,
    #                    const uint8_t *query,
    #                    int l_query, const uint8_t *iqual,
    #                    const probaln_par_t *c,
    #                    int *state, uint8_t *q)

    # # /**********************
    # #  * MD5 implementation *
    # #  **********************/

    # ctypedef struct hts_md5_context

    # # /*! @abstract   Initialises an MD5 context.
    # #  *  @discussion
    # #  *    The expected use is to allocate an hts_md5_context using
    # #  *    hts_md5_init().  This pointer is then passed into one or more calls
    # #  *    of hts_md5_update() to compute successive internal portions of the
    # #  *    MD5 sum, which can then be externalised as a full 16-byte MD5sum
    # #  *    calculation by calling hts_md5_final().  This can then be turned
    # #  *    into ASCII via hts_md5_hex().
    # #  *
    # #  *    To dealloate any resources created by hts_md5_init() call the
    # #  *    hts_md5_destroy() function.
    # #  *
    # #  *  @return     hts_md5_context pointer on success, NULL otherwise.
    # #  */
    # hts_md5_context *hts_md5_init()

    # # /*! @abstract Updates the context with the MD5 of the data. */
    # void hts_md5_update(hts_md5_context *ctx, const void *data, unsigned long size)

    # # /*! @abstract Computes the final 128-bit MD5 hash from the given context */
    # void hts_md5_final(unsigned char *digest, hts_md5_context *ctx)

    # # /*! @abstract Resets an md5_context to the initial state, as returned
    # #  *            by hts_md5_init().
    # #  */
    # void hts_md5_reset(hts_md5_context *ctx)

    # # /*! @abstract Converts a 128-bit MD5 hash into a 33-byte nul-termninated
    # #  *            hex string.
    # #  */
    # void hts_md5_hex(char *hex, const unsigned char *digest)

    # # /*! @abstract Deallocates any memory allocated by hts_md5_init. */
    # void hts_md5_destroy(hts_md5_context *ctx)

    # int hts_reg2bin(int64_t beg, int64_t end, int min_shift, int n_lvls)
    # int hts_bin_bot(int bin, int n_lvls)

    # # * Endianness *
    # int ed_is_big()
    # uint16_t ed_swap_2(uint16_t v)
    # void *ed_swap_2p(void *x)
    # uint32_t ed_swap_4(uint32_t v)
    # void *ed_swap_4p(void *x)
    # uint64_t ed_swap_8(uint64_t v)
    # void *ed_swap_8p(void *x)

####################
### htslib/sam.h ###
####################
cdef extern from "htslib/sam.h" nogil:
    #**********************
    #*** SAM/BAM header ***
    #**********************

    ctypedef struct sam_hrecs_t

    # @abstract Structure for the alignment header.
    # @field n_targets   number of reference sequences
    # @field l_text      length of the plain text in the header
    # @field target_len  lengths of the reference sequences
    # @field target_name names of the reference sequences
    # @field text        plain text
    # @field sdict       header dictionary

    ctypedef struct sam_hdr_t:
         int32_t n_targets, ignore_sam_err
         uint32_t l_text
         uint32_t *target_len
         uint8_t *cigar_tab
         char **target_name
         char *text
         void *sdict
         sam_hrecs_t *hrecs
         uint32_t ref_count
    
    ctypedef sam_hdr_t bam_hdr_t

    #****************************
    #*** CIGAR related macros ***
    #****************************

    int BAM_CMATCH
    int BAM_CINS
    int BAM_CDEL
    int BAM_CREF_SKIP
    int BAM_CSOFT_CLIP
    int BAM_CHARD_CLIP
    int BAM_CPAD
    int BAM_CEQUAL
    int BAM_CDIFF
    int BAM_CBACK

    char    *BAM_CIGAR_STR
    int      BAM_CIGAR_SHIFT
    uint32_t BAM_CIGAR_MASK
    uint32_t BAM_CIGAR_TYPE

    char bam_cigar_op(uint32_t c)
    uint32_t bam_cigar_oplen(uint32_t c)
    char bam_cigar_opchr(uint32_t)
    uint32_t bam_cigar_gen(char, uint32_t)
    int bam_cigar_type(char o)

    # @abstract the read is paired in sequencing, no matter whether it is mapped in a pair
    int BAM_FPAIRED
    # @abstract the read is mapped in a proper pair
    int BAM_FPROPER_PAIR
    # @abstract the read itself is unmapped; conflictive with BAM_FPROPER_PAIR
    int BAM_FUNMAP
    # @abstract the mate is unmapped
    int BAM_FMUNMAP
    # @abstract the read is mapped to the reverse strand
    int BAM_FREVERSE
    # @abstract the mate is mapped to the reverse strand
    int BAM_FMREVERSE
    # @abstract this is read1
    int BAM_FREAD1
    # @abstract this is read2
    int BAM_FREAD2
    # @abstract not primary alignment
    int BAM_FSECONDARY
    # @abstract QC failure
    int BAM_FQCFAIL
    # @abstract optical or PCR duplicate
    int BAM_FDUP
    # @abstract supplementary alignment
    int BAM_FSUPPLEMENTARY

    #*************************
    #*** Alignment records ***
    #*************************

    # @abstract Structure for core alignment information.
    # @field  pos     0-based leftmost coordinate
    # @field  tid     chromosome ID, defined by sam_hdr_t
    # @field  bin     bin calculated by bam_reg2bin()
    # @field  qual    mapping quality
    # @field  l_extranul length of extra NULs between qname & cigar (for alignment)
    # @field  flag    bitwise flag
    # @field  l_qname length of the query name
    # @field  n_cigar number of CIGAR operations
    # @field  l_qseq  length of the query sequence (read)
    # @field  mtid    chromosome ID of next read in template, defined by sam_hdr_t
    # @field  mpos    0-based leftmost coordinate of next read in template
    # @field  isize   observed template length ("insert size")
    ctypedef struct bam1_core_t:
        int32_t     pos
        int32_t     tid
        uint16_t    bin
        uint8_t     qual
        uint8_t     l_extranul
        uint16_t    flag
        uint8_t     l_qname
        uint32_t    n_cigar
        int32_t     l_qseq
        int32_t     mtid
        int32_t     mpos
        int32_t     isize

    # @abstract Structure for one alignment.
    # @field  core       core information about the alignment
    # @field  id
    # @field  data       all variable-length data, concatenated; structure: qname-cigar-seq-qual-aux
    # @field  l_data     current length of bam1_t::data
    # @field  m_data     maximum length of bam1_t::data
    # @field  mempolicy  memory handling policy, see bam_set_mempolicy()
    #
    # @discussion Notes:
    # 1. The data blob should be accessed using bam_get_qname, bam_get_cigar,
    #    bam_get_seq, bam_get_qual and bam_get_aux macros.  These returns pointers
    #    to the start of each type of data.
    # 2. qname is terminated by one to four NULs, so that the following
    #    cigar data is 32-bit aligned; core.l_qname includes these trailing NULs,
    #    while core.l_extranul counts the excess NULs (so 0 <= l_extranul <= 3).
    # 3. Cigar data is encoded 4 bytes per CIGAR operation.
    #    See the bam_cigar_* macros for manipulation.
    # 4. seq is nibble-encoded according to bam_nt16_table.
    #    See the bam_seqi macro for retrieving individual bases.
    # 5. Per base qualities are stored in the Phred scale with no +33 offset.
    #    Ie as per the BAM specification and not the SAM ASCII printable method.
    ctypedef struct bam1_t:
        bam1_core_t core
        uint64_t    id
        uint8_t     *data
        int         l_data
        uint32_t    m_data
        uint32_t    mempolicy

    # @abstract  Get whether the query is on the reverse strand
    # @param  b  pointer to an alignment
    # @return    boolean true if query is on the reverse strand
    int bam_is_rev(bam1_t *b)

    # @abstract  Get the name of the query
    # @param  b  pointer to an alignment
    # @return    pointer to the name string, null terminated
    char *bam_get_qname(bam1_t *b)

    # @abstract  Get the CIGAR array
    # @param  b  pointer to an alignment
    # @return    pointer to the CIGAR array
    #
    # @discussion In the CIGAR array, each element is a 32-bit integer. The
    # lower 4 bits gives a CIGAR operation and the higher 28 bits keep the
    # length of a CIGAR.
    uint32_t *bam_get_cigar(bam1_t *b)

    # @abstract  Get query sequence
    # @param  b  pointer to an alignment
    # @return    pointer to sequence
    #
    # @discussion Each base is encoded in 4 bits: 1 for A, 2 for C, 4 for G,
    # 8 for T and 15 for N. Two bases are packed in one byte with the base
    # at the higher 4 bits having smaller coordinate on the read. It is
    # recommended to use bam_seqi() macro to get the base.
    uint8_t *bam_get_seq(bam1_t *b)

    # @abstract  Get query quality
    # @param  b  pointer to an alignment
    # @return    pointer to quality string
    uint8_t *bam_get_qual(bam1_t *b)

    # @abstract  Get a base on read
    # @param  s  Query sequence returned by bam1_seq()
    # @param  i  The i-th position, 0-based
    # @return    4-bit integer representing the base.
    char bam_seqi(char *s, int i)

    #**************************
    #*** Exported functions ***
    #**************************

    #***************
    #*** BAM I/O ***
    #***************

    ## /* Header */ ##
    # sam_hdr_t *sam_hdr_init(void)
    # sam_hdr_t *bam_hdr_read(BGZF *fp)

    void sam_hdr_destroy(sam_hdr_t *h)
    sam_hdr_t* sam_hdr_dup(const sam_hdr_t *h0)
    ctypedef htsFile samFile
    sam_hdr_t *sam_hdr_read(samFile *fp)
    int sam_hdr_write(samFile *fp, const sam_hdr_t *h)

    # int sam_hdr_nref(const sam_hdr_t *h)
    # int sam_hdr_name2tid(sam_hdr_t *h, const char *ref)

    const char *sam_hdr_tid2name(const sam_hdr_t *h, int tid)
    int sam_hdr_tid2len(const sam_hdr_t *h, int tid)

    ## /* Alignment */ ##
    bam1_t *bam_init1()
    void bam_destroy1(bam1_t *b)
    int bam_read1(BGZF *fp, bam1_t *b)

    # bam1_t *bam_copy1(bam1_t *bdst, const bam1_t *bsrc)
    # bam1_t *bam_dup1(const bam1_t *bsrc)
    # int bam_cigar2qlen(int n_cigar, const uint32_t *cigar)
    # int bam_cigar2rlen(int n_cigar, const uint32_t *cigar)
    # int bam_endpos(const bam1_t *b)
    # int bam_set_qname(bam1_t *b, const char *qname)

    #*************************
    #*** BAM/CRAM indexing ***
    #*************************

    # Load a specific BAM (.csi or .bai) or CRAM (.crai) index file
    # @param fp     File handle of the data file whose index is being opened
    # @param fn     BAM/CRAM/etc data file filename
    # @param fnidx  Index filename, or NULL to search alongside @a fn
    # @return  The index, or NULL if an error occurred.
    hts_idx_t *sam_index_load2(htsFile *fp, const char *fn, const char *fnidx)
    
    # Free a SAM iterator
    # @param iter   Iterator to free
    void sam_itr_destroy(hts_itr_t *iter)

    # Create a BAM/CRAM iterator
    hts_itr_t *sam_itr_queryi(const hts_idx_t *idx, int tid, int beg, int end)

    # # Get the next read from a BAM/CRAM multi-iterator
    # int sam_itr_multi_next(htsFile *htsfp, hts_itr_t *itr, bam1_t *r)

    #***************
    #*** SAM I/O ***
    #***************

    htsFile *sam_open(const char *fn, const char *mode)
    int sam_close(htsFile *fp)

    # sam_read1 - Read a record from a file
    # @param fp   Pointer to the source file
    # @param h    Pointer to the header previously read (fully or partially)
    # @param b    Pointer to the record placeholder
    # @return >= 0 on successfully reading a new record, -1 on end of stream, < -1 on error
    int sam_read1(htsFile *fp, sam_hdr_t *h, bam1_t *b)

    # sam_write1 - Write a record to a file
    # @param fp    Pointer to the destination file
    # @param h     Pointer to the header structure previously read
    # @param b     Pointer to the record to be written
    # @return >= 0 on successfully writing the record, -1 on error
    int sam_write1(samFile *fp, const sam_hdr_t *h, const bam1_t *b)

    #*************************************
    #*** Manipulating auxiliary fields ***
    #*************************************

    # Return a pointer to an aux record
    # @param b   Pointer to the bam record
    # @param tag Desired aux tag
    # @return Pointer to the tag data, or NULL if tag is not present or on error
    # If the tag is not present, this function returns NULL and sets errno to
    # ENOENT.  If the bam record's aux data is corrupt (either a tag has an
    # invalid type, or the last record is incomplete) then errno is set to
    # EINVAL and NULL is returned.
    uint8_t *bam_aux_get(const bam1_t *b, const char tag[2])

    # Get an integer aux value
    # @param s Pointer to the tag data, as returned by bam_aux_get()
    # @return The value, or 0 if the tag was not an integer type
    # If the tag is not an integer type, errno is set to EINVAL.  This function
    # will not return the value of floating-point tags.
    int64_t  bam_aux2i(const uint8_t *s)

    # Get an integer aux value
    # @param s Pointer to the tag data, as returned by bam_aux_get()
    # @return The value, or 0 if the tag was not an integer type
    # If the tag is not an numeric type, errno is set to EINVAL.  The value of
    # integer flags will be returned cast to a double.
    double   bam_aux2f(const uint8_t *s)