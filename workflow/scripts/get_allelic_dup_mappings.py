'''
process sam file reads. output relation between sample core duplicons.
'''
import pysam
import numpy as np
import pandas as pd 


bam_path= 
core_dup_path=

bam = pysam.AlignmentFile(bam_path, "r")
#get references and sort them
aln_refs = list(bam.references)
aln_refs.sort()

def get_q_map(seg , ref_coord):
    '''given segment and a reference coordinate, return the query coordinate that maps to that location.'''
    ####_______________ Parsing Cigar String
    M=0 #M  BAM_CMATCH      0
    I=1 #I  BAM_CINS        1
    D=2 #D  BAM_CDEL        2
    N=3 #N  BAM_CREF_SKIP   3
    S=4 #S  BAM_CSOFT_CLIP  4
    H=5 #H  BAM_CHARD_CLIP  5
    P=6 #P  BAM_CPAD        6
    E=7 #=  BAM_CEQUAL      7
    X=8 #X  BAM_CDIFF       8
    B=9 #B  BAM_CBACK       9
    NM=10 #NM       NM tag  10
    conRef  =       [M, D, N, E, X] # these ones "consume" the reference
    conQuery=       [M, I, S, E, X] # these ones "consume" the query
    conAln  =       [M, I, D, N, S, E, X] # these ones "consume" the alignments
    ####________________

    r_loc = seg.reference_start
    q_step = 0  #how far query coordinate will be from query start
    for opt, l in seg.cigartuples:
        if(opt in conRef):
            l = min(l, ref_coord - r_loc)
            r_loc += l
        if(opt in conQuery):
            q_step += l
        if(r_loc == ref_coord):
            q_coord = seg.qend - q_step if seg.is_reverse else seg.qstart + q_step
            return q_coord
    raise Exception(f"Got to end of cigar tuple without reaching ref_coordinate {seg.reference} : {ref_coord}")

def get_dup(dup_bed_df , coord_start , coord_end):
    '''identify bed duplicon coordinates corresponding to a given start and end, or return -1
    @output: subset of dup_bed_df for core dups corresponding to the region of interest
    '''
    coords = pd.Interval( coord_start , coord_end )
    bed_intervals = [pd.Interval(x,y, 'both') for x,y in dup_bed_df.loc[:, ["start", "stop"]].values ]
    overlaps =  [coords.overlaps(x) for x in bed_intervals]
    if any(overlaps):
        overlaps = np.where( [coords.overlaps(x) for x in bed_intervals] )[0].tolist()
        return dup_bed_df.loc[overlaps , :]
    return -1

def get_score(seg):
    return [tag[1] for tag in seg.get_tags() if tag[0] == 'AS'][0]

def process_segment(seg, ref_dup_bed_df , q_dup_bed_df ):
    '''for each segment, identify reference query coreDup annotations.'''
    nested_dups = get_dup(ref_dup_bed_df , seg.reference_start, seg.reference_end)
    for i, row in nested_dups.iterrows():
        q_1 = get_q_map(seg, row.start)
        q_2 = get_q_map(seg, row.stop)
        q_start = q_1 if not seg.is_reverse else q_2 
        q_end = q_2 if not seg.is_reverse else q_1
        q_dup = get_dup(q_dup_bed_df, q_start , q_end)
        if(q_dup == -1):
            return -1
        assert len(q_dup) == 1 , f"Issue: getting {q_dup} query duplicon alignments to ref dup: {row} "
        #get scores.
        score = get_score(seg)
        return pd.Series( data = [seg.reference_name , row['start'], row['stop'] , row['name'], 
                                  seg.qname, q_dup['start'], q_dup['stop'] , q_dup['name'] , score ] ,
                                  index = ['ref', 'ref_start', 'ref_stop', 'ref_loc_name', 'q', 'q_start', 'q_stop', 'q_dup_name', score ] )


for cur_ref in aln_refs:
    sub_bam = bam.fetch(reference = cur_ref) # returns iterator of alignment file
    dup_aln_stats = map(process_read, sub_bam)
    aggro_dup_results.append(aln_stats)





