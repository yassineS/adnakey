import os
import sys
import math

from cosmos.lib.ezflow.dag  import DAG, add_, split_, sequence_, map_, reduce_, reduce_split_, apply_
from cosmos.lib.ezflow.tool import INPUT

from genomekey.tools        import pipes
from genomekey              import settings

def CteamPipeline(input_bams):

    bam_seq = None
    bam_dup = []  # to check duplicate input files
    
    for b in input_bams:
        # extract genome_id from file, add as a tag
        genome_id = os.path.basename(b).partition('.')[0]

        if genome_id in bam_dup:
            print '\n\nERROR: \"%s\" was already included in the input file list.\n' % b
            sys.exit()
        else:
            bam_dup.append(genome_id)
        
        # TEMPORARILLY, use genome_id as RG_ID, too
        s = sequence_( add_([INPUT(b, tags={'rg':genome_id})], stage_name="Load Input"))

        # append to sequence
        if bam_seq is None:   bam_seq = s
        else:                 bam_seq = sequence_(bam_seq, s, combine=True)

    nInput = len(input_bams)
    nNodes = settings.settings['nNode']
    nSplit = min (256, 16 * max(nNodes /nInput,1)) # will use floor, min 16, up to 256 splits

    settings.settings['nSplit'] = nSplit

    chrom  = ('chrom', range(1,23) + ['X', 'Y', 'MT'])
    split  = ('split', range(1,nSplit+1)) 

    return sequence_(
    bam_seq,
    map_(             pipes.CteamSortSplitBam),      # sort bam by readname (== shuffling)
    split_([split],   pipes.CteamTrimReadGroup),     # 
    map_(             pipes.CteamBwaAln),            # bwa aln
    map_(             pipes.CteamBwaSampe),          # bwa sampe
    reduce_(['rg'],   pipes.CteamSplitByChromosome), # merge split files and (re)split by chromosome
    split_([chrom],   pipes.CteamRmDup_BuildIndex),  # samtools rmdup + index
    map_(             pipes.CteamRealignTarget),     # gatk indel realign target creator
    map_(             pipes.CteamIndelRealigner),    # gatk indel realigner
    map_(             pipes.CteamUnifiedGenotyper)   # gatk unifiedGenotyper

        # #map_(pipes.CteamVariantFiltration)          # gatk variantFilter
    )
