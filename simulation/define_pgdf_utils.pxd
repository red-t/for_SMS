cpdef set get_germline_pos(int ngermline, int nsomatic)

cpdef define_header(str ref_fa, str te_fa, int germline_count, int somatic_count, set germline_pos, float d_rate)

cpdef define_body(dict id2dsl, int popsize, int ntotal, str contig_id, int mindist, int maxdist, set germline_pos, list ins_ids)