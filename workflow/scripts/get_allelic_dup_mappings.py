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
            sub_bed['start'] = sub_bed['start'].map(lambda x: max(x, seg.reference_start) )
            sub_bed['stop'] = sub_bed['stop'].map(lambda x: min(x, seg.reference_end) )
        return sub_bed
    return None

def get_flank_dup(dup_bed_df, cur_dup_i, seg ):
    '''Extend and return reference coordinates containing and spanning to adjacent duplicons or end of alignment.
       @input : dup_bed_df : bedframe generated from script: scripts/get_cluster_annotations.py
       @input : cur_dup_i : index of dup row in dup_bed_df is being looked at.
       @input : seg : pysam AlignmentSegment being processed.
       @output: coords : tuple of start and end coordinates alignemnts (reference coordinates)
    '''
    coord_start = seg.reference_start if cur_dup_i == 0 else dup_bed_df.loc[i-1, 'end']
    coord_end = seg.reference_end if cur_dup_i == dup_bed_df.shape[0] - 1 else dup_bed_df.loc[i-1, 'start']
    return (coord_start , coord_end)


def get_score(seg):
    '''process and return alignment score from a pysam AlignedSegment. Return integer score.'''
    return [tag[1] for tag in seg.get_tags() if tag[0] == 'AS'][0]


def get_dup_bed_df(samp, hap , file_path_sep = "_"):
    '''find and return name of locus bed for a sample haplotype. Or return None if not present'''
    bed_match = [x for x in os.listdir(locus_bed_dir) if f"{samp}{file_path_sep}{hap}" in x]
    if bed_match:
        assert len(bed_match) == 1 , f"Didn't find only one bed file to match sample hap. found: {len(bed_match)}"
        bed_df = pd.read_csv( f"{locus_bed_dir}/{bed_match[0]}", sep = "\t", header = None, names = ['sample', 'start', 'stop', 'name'] )
        return bed_df
    return None


def process_segment(seg ):
    '''for each segment, identify reference query coreDup annotations.'''
    ref_samp, ref_hap = seg.reference_name.split("__")[0] , seg.reference_name.split("__")[1]
    ref_dup_bed_df = get_dup_bed_df( ref_samp,  ref_hap, file_path_sep = "_" )
    assert ref_dup_bed_df is not None , f"Couldn't find reference locus bed: {ref_samp} , {ref_hap} in {locus_bed_dir }. Check file_path_sep char."
    nested_dups = get_dup(ref_dup_bed_df , seg.reference_start, seg.reference_end, is_seg = True)
    if(nested_dups is None):
        return None
    aln_stats = pd.DataFrame(columns =  ['ref', 'ref_start', 'ref_stop', 'ref_loc_name', 'q', 'q_start', 'q_stop', 'q_dup_name', 'score' ])
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
        score = get_score(seg)
        aln_summary = pd.Series( data = [seg.reference_name , row['start'], row['stop'] , row['name'], 
                                    seg.qname, q_dup['start'], q_dup['stop'] , q_dup['name'] , score ] ,
                                    index = list(aln_stats.columns.values))
        aln_stats = aln_stats.append(aln_summary, ignore_index=True)
    return aln_stats

#def main():
aggro_aln_stats = pd.DataFrame(columns =  ['ref', 'ref_start', 'ref_stop', 'ref_loc_name', 'q', 'q_start', 'q_stop', 'q_dup_name', 'score' ])
for cur_ref in aln_refs:
    sub_bam = bam.fetch(reference = cur_ref) # returns iterator of alignment file
    alignments = list(sub_bam)
    #dup_aln_stats = map(process_segment, sub_bam) 
    aln_stats = pd.DataFrame(columns =  ['ref', 'ref_start', 'ref_stop', 'ref_loc_name', 'q', 'q_start', 'q_stop', 'q_dup_name', 'score' ])
    for seg in alignments:
        print(f"processing: {seg.reference_name}:{seg.reference_start}-{seg.reference_end} ;; reverse: {seg.is_reverse} ;; {seg.qname}:{seg.qstart}-{seg.qend}")
        try:
            cur_aln_stats = process_segment(seg)
            if(aln_stats is not None):
                aln_stats = aln_stats.append(cur_aln_stats, ignore_index=True)
        except:
            e = sys.exc_info()[0]
            print("ERROR")
            print(e)
    aggro_aln_stats = aggro_aln_stats.append(aln_stats)

aggro_aln_stats.to_csv(snakemake.output.locus_mappings , sep='\t' , header = True, index = False , na_rep='NA')

#for main
# if(aggro_aln_stats.shape[0] > 0):
#     return aggro_aln_stats
# return -1
    



