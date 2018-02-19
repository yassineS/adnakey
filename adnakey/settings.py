import os

#from cosmos.config import settings as cosmos_settings

#import pwd
#userName=pwd.getpwuid(os.getuid()).pw_name

def _get_specification_GE(jobAttempt):

    cpu_req  = jobAttempt.task.cpu_requirement
    #mem_req = jobAttempt.task.memory_requirement

    #return '-l spock_mem={0}M,num_proc={1}'.format(mem_req,cpu_req)
    return '-l num_proc={}'.format(cpu_req)

opj = os.path.join

tools_path = '/WGA/tools/'     #'/groups/cbi/WGA/tools/cteam'     in Orchestra
ref_path   = '/WGA/reference/' #'/groups/cbi/WGA/reference/cteam'


settings = {
    'java'                  : 'java -d64 -XX:+AggressiveOpts -XX:+UseLargePages',
    
    # Output directory can be S3
    'outputDir'             : 's3://lpm.reichkey.test',
    'scratch'               : '/mnt/tmp',
    
    'samtools'              : opj(tools_path, 'samtools.v0.1.18_r580'),
    'bwa'                   : opj(tools_path, 'bwa.v0.5.10'),
    'gatk'                  : opj(tools_path, 'gatk.v2.5-2.jar'),
    
    'removeTagMapping'      : opj(tools_path, 'removeTagsMapping'),  
    'addRgCteam'            : opj(tools_path, 'addRG_CTEAM'),
    'mergeTrimReadsBam'     : opj(tools_path, 'mergeTrimReadsBAM'),
    
    'ref1'                  : opj(ref_path,   'hs37d5.fa.gz'),
    'ref2'                  : opj(ref_path,   'hs37d5.fa'),
    'dbsnp_vcf'             : opj(ref_path,   'dbsnp_137.b37.vcf'),
    'mask'                  : opj(ref_path,   'hs37m_mask35_50.flt.bed'),
    
    'drmaa_spec'            : _get_specification_GE,

    # should be in the same directory with this file
    'no_et_key'             : opj(os.path.dirname(os.path.realpath(__file__)),'no_et.key'),

    # Cluster specification
    'nNode'               : 1,   # total number of computing nodes in the cluster
    'nCpuPerNode'         : 32,  # number of cores in a node
    'memPerNode'          : 200  # memory in GB per node
    }
