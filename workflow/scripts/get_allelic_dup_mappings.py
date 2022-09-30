'''
process sam file reads. output relation between sample core duplicons.
'''
import pysam
import numpy as np
import pandas as pd 
import os
import sys
#import functions and classes from repo process_cigar
sys.path.append('/net/eichler/vol26/home/guitarfx/software/github_clones/alignment_operations')
from cigar_opt_class import CigarOperation
from cigar_opt_class import cigar_dict
import process_cigar
#import functions and classes from repo process_cigar
from multiprocessing import cpu_count
from multiprocessing.pool import Pool
import time

pd.options.mode.chained_assignment = None  # default='warn'


bam_path =  snakemake.input.bam     #bam file being processed
locus_bed_dir = snakemake.params.locus_bed_dir   #directory holding all dup_beds for samples in bam.

bam = pysam.AlignmentFile(bam_path, "r")
#get references and sort them
aln_refs = list(bam.references)
aln_refs.sort()

def get_q_map(seg , ref_coord):
    '''given segment and a reference coordinate, return the query coordinate that maps to that location.'''
    return process_cigar.get_qCoord_from_rCoord( seg , ref_coord )
    

def get_dup(dup_bed_df , coord_start , coord_end, is_seg = False):
    '''identify bed duplicon coordinates corresponding to a given start and end, or return None
    @output: subset of dup_bed_df for core dups corresponding to the region of interest
    @input: is_seg: if coord start and end are the boundaries of an alignment segment, then need to adjust dup bed to not go beyond alignment.
    '''
    assert coord_start <= coord_end , f"coordinates wrong: start > end. coord_start: {coord_start} ; coord_end : {coord_end}"
    coords = pd.Interval( coord_start , coord_end )
    bed_intervals = [pd.Interval(x,y, 'both') for x,y in dup_bed_df.loc[:, ["start", "stop"]].values ]
    overlaps =  [coords.overlaps(x) for x in bed_intervals]
    if any(overlaps):
        sub_bed = dup_bed_df.loc[overlaps , :]
        if(is_seg):
            sub_bed['start'] = sub_bed['start'].map(lambda x: max(x, coord_start ) )
            sub_bed['stop'] = sub_bed['stop'].map(lambda x: min(x, coord_end ) )
        return sub_bed
    return None

def get_dup_flank_coords(dup_bed_df, cur_dup_i, seg ):
    '''Extend and return reference coordinates containing and spanning to adjacent duplicons or end of alignment.
       @input : dup_bed_df : bedframe generated from script: scripts/get_cluster_annotations.py
       @input : cur_dup_i : index of dup row in dup_bed_df is being looked at.
       @input : seg : pysam AlignmentSegment being processed.
       @output: coords : tuple of start and end coordinates alignemnts (reference coordinates)
    '''
    dup_bed_df = dup_bed_df.sort_values(by = 'start') #make sure it's sorted
    coord_start = seg.reference_start if cur_dup_i == 0 else dup_bed_df.loc[cur_dup_i-1, 'stop']
    coord_start = max(seg.reference_start, coord_start)
    coord_end = seg.reference_end if cur_dup_i == dup_bed_df.shape[0] - 1 else dup_bed_df.loc[cur_dup_i+1, 'start']
    coord_end = min(seg.reference_end, coord_end )
    return (coord_start , coord_end)


def get_score(seg, ref_start = None, ref_stop = None):
    '''process and return alignment score from a pysam AlignedSegment. Return integer score.'''
    if( ref_start is not None and ref_stop is not None ):
        subset_cigar_tuples = process_cigar.subset_cigar(seg ,ref_start , ref_stop)
        return process_cigar.get_cigarTuple_alignment_score( subset_cigar_tuples )
    else:
        return process_cigar.get_cigarTuple_alignment_score(seg.cigartuples)


def get_dup_bed_df(samp, hap , file_path_sep = "_"):
    '''find and return name of locus bed for a sample haplotype. Or return None if not present'''
    bed_match = [x for x in os.listdir(locus_bed_dir) if f"{samp}{file_path_sep}{hap}" in x]
    if bed_match:
        assert len(bed_match) == 1 , f"Didn't find only one bed file to match sample hap. found: {len(bed_match)}"
        bed_df = pd.read_csv( f"{locus_bed_dir}/{bed_match[0]}", sep = "\t", header = None )
        bed_df = bed_df.rename( columns = dict( zip( range(0,6) , ['sample', 'start', 'stop', 'name', '.', 'strand']  ) ) )  #change column namee 
        return bed_df
    return None

def process_segment(seg ):
    '''for each segment, identify reference query coreDup annotations.'''
    ref_samp, ref_hap = seg.reference_name.split("__")[0] , seg.reference_name.split("__")[1]
    ref_dup_bed_df = get_dup_bed_df( ref_samp,  ref_hap, file_path_sep = "_" )
    ref_dup_bed_df = ref_dup_bed_df.sort_values(by = "start").reset_index(drop= True)
    assert ref_dup_bed_df is not None , f"Couldn't find reference locus bed: {ref_samp} , {ref_hap} in {locus_bed_dir }. Check file_path_sep char."
    nested_dups = get_dup(ref_dup_bed_df , seg.reference_start, seg.reference_end, is_seg = True)
    if(nested_dups is None):
        return None
    aln_stats = pd.DataFrame(columns =  ['ref', 'ref_start', 'ref_stop', 'ref_loc_name', 'r_score_start', 'r_score_stop' , 'q', 'q_start', 'q_stop', 'q_dup_name', 'q_score_start', 'q_score_stop' ,'score', 'strand' ]) #add coordinates of flank.
    for i, row in nested_dups.iterrows():
        q_samp , q_hap = seg.qname.split("__")[0] , seg.qname.split("__")[1]
        q_dup_bed_df = get_dup_bed_df(q_samp, q_hap, file_path_sep = "_")
        assert q_dup_bed_df is not None , f"Couldn't find query locus bed: {q_samp} , {q_hap} in {locus_bed_dir }. Check file_path_sep char."

        q_1 , q_2 = get_q_map(seg, row.start) , get_q_map(seg, row.stop)
        q_start = q_1 if not seg.is_reverse else q_2
        q_end = q_2 if not seg.is_reverse else q_1

        q_dup = get_dup(q_dup_bed_df, q_start , q_end, is_seg = False) ######!!!!!
        if(q_dup is None):
            continue
        assert len(q_dup) == 1 , f"Issue: getting {q_dup} query duplicon alignments to ref dup: {row} "
        q_dup = q_dup.iloc[0] #turn into pandas Series
        #get scores.
        flank_tuple = get_dup_flank_coords( ref_dup_bed_df , i , seg ) #get coordinates of flanks of dup.
        score = get_score(seg, flank_tuple[0] , flank_tuple[1])
        r_flank_start , r_flank_stop = flank_tuple[0] , flank_tuple[1]
        q_flank_1, q_flank_2 = get_q_map(seg, r_flank_start ) , get_q_map(seg, r_flank_stop)
        q_flank_start = q_flank_2 if seg.is_reverse else q_flank_1
        q_flank_stop = q_flank_1 if seg.is_reverse else q_flank_2
        #get strand
        strand = '-' if seg.is_reverse else '+'
        assert score <= flank_tuple[1] - flank_tuple[0] , f"Issue: problem with segment: getting higher score than size of locus: score {score}  flank_tuple : {flank_tuple}"
        aln_summary = pd.Series( data = [seg.reference_name , row['start'], row['stop'] , row['name'], r_flank_start, r_flank_stop ,
                                    seg.qname, q_dup['start'], q_dup['stop'] , q_dup['name'] , q_flank_start, q_flank_stop , score  , strand] ,
                                    index = list(aln_stats.columns.values))
        aln_stats = pd.concat([aln_stats, aln_summary.to_frame().T])   
        #aln_stats = aln_stats.append(aln_summary, ignore_index=True)
    return aln_stats

def try_process_seg(segment):
    '''process segment or return error
       @input: pysam AlignmentSegment object
       @output: result from process_segment() or None
    '''
    try:
        cur_aln_stats = process_segment(segment)
        if(cur_aln_stats is not None):
            return cur_aln_stats
        # if(cur_aln_stats is not None):
        #     aln_stats = aln_stats.append(cur_aln_stats, ignore_index=True)
    except:
        e = sys.exc_info()[0]
        print("ERROR")
        print(e)
    return None

def get_aln_stats(alignments):
    '''logic to map try_process_seg across a set of alignments'''
    x = map( try_process_seg , alignments)
    return pd.concat(x , axis = 0)

def open_pysam(filepath_or_object, mode="r"):
    '''open a pysam object in read mode'''
    return pysam.AlignmentFile(filepath_or_object, mode)

def process_ref(cur_ref , bam_path = bam_path ):
    '''process bam reads for one reference'''
    with open_pysam(bam_path) as bam:
        alignments = list( bam.fetch(reference = cur_ref) )
    return get_aln_stats(alignments)

start = time.time()
pool = Pool(processes = min( snakemake.threads - 1 , cpu_count()-1 ) )
x = pool.map( process_ref, aln_refs , chunksize = 1 ) #lambda x : process_ref(x , bam) , aln_refs[0:3] ) 
pool.close()
pool.join()
mapped_alignments = pd.concat(x , axis = 0)
end = time.time()
print ("Time elapsed:", end - start)

mapped_alignments.to_csv(snakemake.output.locus_mappings , sep='\t' , header = True, index = False , na_rep='NA')

    



