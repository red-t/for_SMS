#include "AIList.h"

ailist_t *ailist_init(void)
{
	ailist_t *ail = malloc(1*sizeof(ailist_t));
	ail->ctg = malloc(1*sizeof(ctg_t)); // allocate memory for 1 contig
    ctg_t *p = &ail->ctg[0]; // ctg_t pointer to the first contig
    p->nr = 0; // initialize nr(number of regions) of the first contig to 0
    p->mr = 64;
	p->glist = malloc(p->mr*sizeof(gdata_t)); // allocate memory for mr regions
	return ail;
}


void ailist_destroy(ailist_t *ail)
{
	if (ail == 0) return;
    free(ail->ctg[0].glist); // firstly free the list of regions
	free(ail->ctg[0].maxE); // then free the list of max end
	free(ail->ctg); // then free the contig
	free(ail); // finally free the AIList
}


void ailist_add(ailist_t *ail, int32_t s, int32_t e)
{
	if(s > e) return;
	ctg_t *q = &ail->ctg[0];
	if(q->nr == q->mr) EXPAND(q->glist, q->mr);	// if the region list is full, expand mr more
	gdata_t *p = &q->glist[q->nr++]; // pointer to the nr+1-th region
	p->start = s;
	p->end   = e;
    // TO DO: should add a new feature value, representing the TE family
	return;
}


void readBED(ailist_t *ail, const char* fn, const char* chr)
{
    char buf[1024];
    char *s1, *s2, *s3;
    FILE* fd = fopen(fn, "r");
    while(fgets(buf, 1024, fd)){
        s1 = strtok(buf, "\t");
        s2 = strtok(NULL, "\t");
        s3 = strtok(NULL, "\t");
        // TO DO: need one more strtok to get the information of the TE family,
        // then change it to value convert it to integer value through a dictionary
		if(s1){
            // if the contig of the region is the same as the
            // given chr, then add the region into region list
            int32_t ret = strcmp(s1, chr);
            if(ret == 0) ailist_add(ail, atol(s2), atol(s3));
        }
	}
	fclose(fd);
	return;
}


void ailist_construct(ailist_t *ail, int cLen)
{   //New continueous memory?
    int cLen1=cLen/2, j1, minL = MAX(64, cLen);
    cLen += cLen1;
    int lenT, len, iter, j, k, k0, t;
	//1. Decomposition
	ctg_t   *p  = &ail->ctg[0];
	gdata_t *L1 = p->glist;                         //L1: to be rebuilt
	int32_t  nr = p->nr;
    if(nr<=minL){
        p->nc = 1, p->lenC[0] = nr, p->idxC[0] = 0;
    }
    else{
        gdata_t *L0 = malloc(nr*sizeof(gdata_t)); 	//L0: serve as input list
        gdata_t *L2 = malloc(nr*sizeof(gdata_t));   //L2: extracted list
        memcpy(L0, L1, nr*sizeof(gdata_t));
        iter = 0;	k = 0;	k0 = 0;
        lenT = nr;
        while(iter<MAXC && lenT>minL){
            len = 0;
            for(t=0; t<lenT-cLen; t++){
                int32_t tt = L0[t].end;
                j=1;    j1=1;
                while(j<cLen && j1<cLen1){
                    if(L0[j+t].end>=tt) j1++;
                    j++;
                }
                if(j1<cLen1) memcpy(&L2[len++], &L0[t], sizeof(gdata_t));
                else memcpy(&L1[k++], &L0[t], sizeof(gdata_t));
            }
            memcpy(&L1[k], &L0[lenT-cLen], cLen*sizeof(gdata_t));
            k += cLen, lenT = len;
            p->idxC[iter] = k0;
            p->lenC[iter] = k-k0;
            k0 = k, iter++;
            if(lenT<=minL || iter==MAXC-2){			//exit: add L2 to the end
                if(lenT>0){
                    memcpy(&L1[k], L2, lenT*sizeof(gdata_t));
                    p->idxC[iter] = k;
                    p->lenC[iter] = lenT;
                    iter++;
                }
                p->nc = iter;
            }
            else memcpy(L0, L2, lenT*sizeof(gdata_t));
        }
        free(L2),free(L0);
    }
    //2. Augmentation
    p->maxE = malloc(nr*sizeof(int32_t));
    for(j=0; j<p->nc; j++){
        k0 = p->idxC[j];
        k  = k0 + p->lenC[j];
        int32_t tt = L1[k0].end;
        p->maxE[k0] = tt;
        for(t=k0+1; t<k; t++){
            if(L1[t].end > tt) tt = L1[t].end;
            p->maxE[t] = tt;
        }
    }
}


int32_t bSearch(gdata_t* As, int32_t idxS, int32_t idxE, int32_t qe)
{   //find tE: index of the first item satisfying .start<qe from right
    int tL=idxS, tR=idxE-1, tM, tE=-1;
    // if start of the rightmost region less than qe,
    // all the other regions on the left is eligible
    if(As[tR].start < qe)
        return tR;
    // if start of the leftmost region larger than qe,
    // all the other regions on the right is not eligible
    else if(As[tL].start >= qe)
        return -1;
    while(tL<tR-1){ // binary search between tL & tR
        tM = (tL+tR)/2; 
        if(As[tM].start >= qe)
            tR = tM-1;
        else
            tL = tM;
    }
    if(As[tR].start < qe)
        tE = tR;
    else if(As[tL].start < qe)
        tE = tL;       
    return tE; 
}


void query_dist_p(ailist_t *ail, int32_t rpos, int32_t flank, int32_t *n, int32_t *d)
{   
    int32_t nr = 0; // number of regions overlapped with query
    int32_t qs = (rpos<flank) ? 0 : rpos-flank; // query start of the extended rpos
    int32_t qe = rpos + flank; // query end of the extended rpos
    int32_t k, md;
    int32_t mdist = 0x7fffffff;
    ctg_t   *p = &ail->ctg[0]; // point to the first contig
    
    // k-th component catains p->lenC[k] regions
    // search in components one by one
    for(k=0; k<p->nc; k++){
        int32_t cs = p->idxC[k];        // component start index
        int32_t ce = cs + p->lenC[k];   // component end index
        int32_t t;
        if(p->lenC[k]>15){
            t = bSearch(p->glist, cs, ce, qe); 	// get index of the first item satisfying .start<qe from right
            while(t>=cs && p->maxE[t]>qs){      // search from t-th region to right
                if(p->glist[t].end>qs){         // .start<qe & .end>qs, overlap
                    md    = MIN(abs(rpos - p->glist[t].start), \
                                abs(rpos - p->glist[t].end));
                    // TO DO: should get value of the overlapped
                    // and update it if mdist is updated
                	mdist = MIN(mdist, md);     // calculate the shortest distance from the border
                    nr++;
                }
                t--;
            }
        }
        else{
            for(t=cs; t<ce; t++){
                if(p->glist[t].start<qe && p->glist[t].end>qs){
                	md    = MIN(abs(rpos - p->glist[t].start), \
                                abs(rpos - p->glist[t].end));
                    // TO DO: should get value of the overlapped
                    // and update it if mdist is updated
                	mdist = MIN(mdist, md);
                    nr++;
                }
			}
        }
    }
    *d = MIN(mdist, *d);
    *n = *n + nr;
    // return nr;
}


void query_dist_c(ailist_t *ail, int32_t st, int32_t ed, int32_t flank, int32_t *n, int32_t *d)
{
    int32_t nr = 0;
    int32_t qs = (st<flank) ? 0 : st-flank; // query start of the extended st
    int32_t qe = ed + flank; // query end of the extended ed
    int32_t k, ldist, rdist;
    int32_t mdist = 0x7fffffff;
    ctg_t   *p = &ail->ctg[0];

    for(k=0; k<p->nc; k++){
        int32_t cs = p->idxC[k];
        int32_t ce = cs + p->lenC[k];
        int32_t t;
        if(p->lenC[k]>15){
            t = bSearch(p->glist, cs, ce, qe);
            while(t>=cs && p->maxE[t]>qs){
                if(p->glist[t].end>qs){
                    ldist = MIN(abs((int32_t)st - (int32_t)p->glist[t].start), \
                                abs((int32_t)st - (int32_t)p->glist[t].end));
                    rdist = MIN(abs((int32_t)ed - (int32_t)p->glist[t].start), \
                                abs((int32_t)ed - (int32_t)p->glist[t].end));
                    // TO DO: should get value of the overlapped
                    // and update it if mdist is updated
                    mdist = MIN(mdist, MIN(ldist, rdist));
                    nr++;
                }
                t--;
            }
        }
        else{
            for(t=cs; t<ce; t++){
                if(p->glist[t].start<qe && p->glist[t].end>qs){
                	ldist = MIN(abs((int32_t)st - (int32_t)p->glist[t].start), \
                                abs((int32_t)st - (int32_t)p->glist[t].end));
                    rdist = MIN(abs((int32_t)ed - (int32_t)p->glist[t].start), \
                                abs((int32_t)ed - (int32_t)p->glist[t].end));
                    // TO DO: should get value of the overlapped
                    // and update it if mdist is updated
                    mdist = MIN(mdist, MIN(ldist, rdist));
                    nr++;
                }
			}
        }
    }
    *d = MIN(mdist, *d);
    *n = *n + nr;
}


void aln_loc_flag(ailist_t *rep_ail, ailist_t *gap_ail, seg_dtype_struct segs[])
{
    int32_t n1 = 0, n2 = 0;
    int32_t d1 = 0x7fffffff, d2 = 0x7fffffff;
    uint8_t flag1, flag2;

    // intersecting with repeats
    query_dist_p(rep_ail, segs[0].refst, 50, &n1, &d1);
    query_dist_p(rep_ail, segs[0].refed, 50, &n2, &d2);

    // intersecting with gaps
    query_dist_p(gap_ail, segs[0].refst, 50, &n1, &d1);
    query_dist_p(gap_ail, segs[0].refed, 50, &n2, &d2);

    // flag of refst
    if (n1>0){
        if (d1<50) {
            // at boundary
            flag1 = 1;
        } else {
            // inside repeats/gap
            flag1 = 2;
        }
    } else {
        // at normal
        flag1 = 0;
    }

    // flag of refed
    if (n2>0){
        if (d2<50) {
            // at boundary
            flag2 = 1;
        } else {
            // inside repeats/gap
            flag2 = 2;
        }
    } else {
        // at normal
        flag2 = 0;
    }

    // flag of the alignment
    if (flag1==1 && flag2==1){ // both sides at repeat/gap boundary
        segs[0].loc_flag = 4;
    } else if (flag1==2 && flag2==2) { // both sides inside repeat/gap
        segs[0].loc_flag = 2;
    } else if (flag1==1 || flag2==1) { // one side at repeat/gap boundary
        if (flag1==0 || flag2==0) { // the other side at normal region
            segs[0].loc_flag = 8;
        } else { // the other side inside repeat/gap
            segs[0].loc_flag = 16;
        }
    } else { // at least one side at normal region
        segs[0].loc_flag = 1;
    }
}


void clt_loc_flag(ailist_t *rep_ail, ailist_t *gap_ail, cluster_dtype_struct clts[])
{
    int32_t n = 0, d = 0x7fffffff;

    // intersecting with repeats
    query_dist_c(rep_ail, clts[0].st, clts[0].ed, 50, &n, &d);

    // intersecting with gaps
    query_dist_c(gap_ail, clts[0].st, clts[0].ed, 50, &n, &d);

    // flag of the cluster
    if (n>0){
        if (d<50) {
            // at boundary
            clts[0].cloc_flag = 2;
        } else {
            // inside repeats/gap
            clts[0].cloc_flag = 4;
        }
    } else {
        // at normal
        clts[0].cloc_flag = 1;
    }
}