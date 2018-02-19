from cosmos.lib.ezflow.tool import Tool

def _list2input(l, opt):
    return " ".join(map(lambda x: str(x), l))

# Init: Make a temp directory and set traps
cmd_init = r"""
        tmpDir="$( mktemp -d --tmpdir={s[scratch]})" && mkdir -p $tmpDir/out
        trap 'exit 1'              ERR SIGHUP SIGINT SIGTERM SIGQUIT
        trap '/bin/rm -rf $tmpDir' EXIT
        id={p[rg]} && dest={s[outputDir]}
        """

# Download input from S3 bucket
cmd_input_s3 = r"""
        [[ "$input" == s3://* ]] && printf "%s downloading from S3\n" "$(date "+%T %D")" && \
        aws s3 cp "$input" "$tmpDir"/ --quiet                                    && \
        [[ $? -eq 0 ]]           && printf "%s downloading done\n"    "$(date "+%T %D")" && input=$tmpDir/"$( basename $input )"        
        [[ ! -f "$input"      ]] && printf "ERROR: couldn't find %s\n" "$input"  && exit 1
        """

# Move output to S3 or outdir
cmd_out = r"""
        printf "%s ########## Move output files ##########\n" "$(date "+%T %D")"
        for f in $tmpDir/out/*;   # no quotes for proper expansion
        do
           out=$outDir/$( basename "$f" )
           [[ "$dest" == s3://* ]] && aws s3 cp "$f" "$out" --quiet && printf "%s copied %s to %s\n" "$(date "+%T %D")" "$(basename $f)" "$outDir"
           [[ "$dest" == /*     ]] &&        mv "$f" "$out"
           echo "$out" >> $OUT.txt
        done        
        """

# Attach another SSD EBS if space is not enough
cmd_ebs_init = r"""
        printf "%s ########## Attach EBS Drive if needed ##########\n" "$(date "+%T %D")"
        inputSize=$( ls -l $input                           | cut -d' ' -f5 )
         diskSize=$( df | awk '/\/mnt/' | sed 's/[ ]\+/ /g' | cut -d' ' -f4 )
          volSize=$(( 3 * inputSize / 1000000000 )) && [[ $volSize -eq 0 ]] && volSize=1   # min 1 GB (for testing)

        #if [[ $diskSize -lt 10*$inputSize ]] ; then
          printf "%s InputSize= %s, DiskSize= %s => Attach\n" "$(date "+%T %D")" "$inputSize" "$diskSize"

            iid=$( GET http://169.254.169.254/latest/meta-data/instance-id )
           zone=$( GET http://169.254.169.254/latest/meta-data/placement/availability-zone )
          volId=$( aws ec2 create-volume --size $volSize --volume-type gp2 --availability-zone $zone --output text | cut -f8 )
          while $(sleep 5); do
              status=$( aws ec2 describe-volumes --volume-id $volId --filter Name=status,Values=available ); 
              [[ -n "$status" ]] && break; echo "Volume not available yet"
          done

          aws ec2 attach-volume --volume-id $volId --instance-id $iid --device /dev/sdf
          while $(sleep 5); do
              status=$( aws ec2 describe-volumes --volume-id $volId --filter Name=attachment.status,Values=attached ); 
              [[ -n "$status" ]] && break; echo "Volume not attached yet"
          done

          sudo mkfs.ext4 /dev/xvdf -q && sudo mount /dev/xvdf /ebs && sudo chown ubuntu:ubuntu /ebs
          [[ $? -eq 0 ]] && printf "%s new EBS volume %s was attached to /ebs\n" "$(date "+%T %D")" "$volId"
        #fi
        """

# Should be always paired with cmd_ebs_init
cmd_ebs_remove = r"""
        printf "%s ########## Remove EBS Drive ##########\n" "$(date "+%T %D")"
        check=$( mount | grep "/ebs" )
        [[ -n "$check" ]] && sudo umount /ebs
        [[ -n "$volId" ]] && aws ec2 detach-volume --volume-id $volId --force
        while $(sleep 5); do
          status=$( aws ec2 describe-volumes --volume-id $volId --filter Name=status,Values=available ); 
          [[ -n "$status" ]] && break; echo "Volume not detached yet"
        done
        [[ -n "$volId" ]] && aws ec2 delete-volume --volume-id $volId
        """

class CteamSortSplitBam(Tool):
    """
    DIR=/files/Genetics/1000Genomes/cteam/RAW_DATA/hard_drives/disk#/
    IN=$DIR/$ID/Assembly/genome/bam/$ID.bam

    OUTDIR=$DIR/sorted_bams/$ID && mkdir -p $OUTDIR
    OUT=$OUTDIR/$ID.sort   # output prefix

    $SAMTOOLS sort -n $IN $OUT
    """
    
    name    = "SortSplitBAM" # gridengine DRM does NOT allow names starting with numbers
    cpu_req = 30
    mem_req = 200 * 1024     # 200 GB

    inputs  = ['bam']
    outputs = ['txt']        # txt file that has the location of actual data file

    def cmd(self,i,s,p):
        # -n: Sort by read names rather than by chromosomal coordinates.
        # -o: Output the final alignment to the standard output.
        # -m INT: Approximately the maximum required memory per thread. [500000000] 
        # -l: compression level - works?

        # sort -> convert to sam -> cut tags -> split -> convert to bam
        cmd_init_local = r"""
        input={i[bam][0]}
        outDir=$dest/sorted_split_bams && [[ "$dest" != s3://* ]] && mkdir -p $outDir
        printf "%s input=  %s\n\t\t\ttmpDir= %s\n\t\t\toutDir= %s\n" "$(date "+%T %D")" "$input" "$tmpDir" "$outDir"
        """

        cmd_main = r"""
        printf "%s ########## SortSplitBAM Main ##########\n" "$(date "+%T %D")"
        # Build index
        [[ ! -f "$input".idx  ]] && printf "%s building index\n"      "$(date "+%T %D")" && \
        {s[samtools]} index "$input"                                             && \
        [[ $? -eq 0 ]]           && printf "%s building index done\n" "$(date "+%T %D")"

        # Calculate number of reads in each split file
        nRead=$( {s[samtools]} idxstats "$input" | cut -f 3,4 | sed 's/\t/+/;s/^/.+/' | bc | tail -1 )
        lines=$((  (nRead/({s[nSplit]} * 1000)) * 1000 + 1000 ))
        printf "Number of reads in input: %s, Number of reads in each of %s split files: %s\n" "$nRead" "{s[nSplit]}" "$lines"

        # Write output to a tmp directory first, move them to the actual output dir
        tmpOut=$tmpDir/out/$id".split" # prefix

        # -@ > 16 is not recommend as samtools will make >2 temporary disk files per theread.
        # Also, '-l 0' option is not recommended unless disk space > 5x input file.
        # Last, '-S -1' (fast compression) can be used instead of '-Su' (uncompressed)       
        {s[samtools]} sort -o -n -@ 16 -m 12G -l 0 "$input" /ebs/tmpSort |\
        {s[samtools]} view - | cut -f1-11                                |\
        split -l "$lines" --numeric-suffixes --filter='{s[samtools]} view -Su -T {s[ref1]} - > $FILE.bam 2> /dev/null' - "$tmpOut"
        [[ $? -eq 0 ]] && printf "%s splitting input files done.\n" "$(date "+%T %D")"
        """
        
        return (cmd_init + cmd_init_local + cmd_input_s3 + cmd_ebs_init + cmd_main + cmd_out + cmd_ebs_remove)

class CteamTrimReadGroup(Tool):
    """
    IN=$DIR/sorted_bams/$ID/$ID.sort.bam

    OUTDIR=$DIR/merged/$ID && mkdir -p $OUTDIR
    OUT=$OUTDIR/$ID.sort.bam

    $REMOVE_TAG_MAPPING $IN  /dev/stdout |
    $ADD_RG_CTEAM /dev/stdin /dev/stdout $RG_ID |
    $MERGE_TRIM_READS_BAM --keepOrig -o $OUT --log $OUT.trim.log /dev/stdin
    """

    name = "TrimReadGroup"
    cpu_req = 2
    mem_req = 4 * 1024 # 4GB

    inputs  = ['txt']
    outputs = ['txt']
    
    def cmd(self,i,s,p):
        cmd_init_local = r"""
        input=$( head -n {p[split]} {i[txt][0]} | tail -n 1 )
        outDir=$dest/trim && [[ "$dest" != s3://* ]] && mkdir -p $outDir
        printf "%s input= %s\n\t\t\ttmpDir= %s\n\t\t\tsplitId= %s\n" "$(date "+%T %D")" "$input" "$tmpDir" "{p[split]}"
        """
        cmd_main = r"""
        printf "%s ########## TrimReadGroup Main ##########\n" "$(date "+%T %D")"
        tmpOut=$tmpDir/out/$id".{p[split]}".trim.bam

        {s[removeTagMapping]}  $input /dev/stdout |\
        {s[mergeTrimReadsBam]} --keepOrig -o $tmpOut --log $tmpOut".log" /dev/stdin
        [[ $? -eq 0 ]]          && printf "%s trimming %s done.\n" "$(date "+%T %D")" "{p[split]}"
        """

        return (cmd_init + cmd_init_local + cmd_input_s3 + cmd_main + cmd_out)

class CteamBwaAln(Tool):
    """
    IN=$DIR/merged/$ID/$ID.sort.bam
    
    OUTDIR=$DIR/mapped_hg19_v3/$ID && mkdir -p $OUTDIR
    OUT=$OUTDIR/$ID  # output prefix

    ($BWA aln -q 15 -t 2 $REF1 -b1 $IN > $OUT.sai1 ) >& $OUT.sai1.oe 
    ($BWA aln -q 15 -t 2 $REF1 -b2 $IN > $OUT.sai2 ) >& $OUT.sai2.oe 
    """

    name = "BWA_ALN"
    cpu_req = 2
    mem_req = 4 * 1024 # 4 GB
    
    inputs  = ['txt']
    outputs = ['txt']

    def cmd(self,i,s,p):
        # Usage:   bwa aln [options] <prefix> <in.fq>
        
        # Options: -n NUM    max #diff (int) or missing prob under 0.02 err rate (float) [0.04]
        #          -o INT    maximum number or fraction of gap opens [1]
        #          -e INT    maximum number of gap extensions, -1 for disabling long gaps [-1]
        #          -i INT    do not put an indel within INT bp towards the ends [5]
        #          -d INT    maximum occurrences for extending a long deletion [10]
        #          -l INT    seed length [32]
        #          -k INT    maximum differences in the seed [2]
        #          -m INT    maximum entries in the queue [2000000]
        #          -t INT    number of threads [1]
        #          -M INT    mismatch penalty [3]
        #          -O INT    gap open penalty [11]
        #          -E INT    gap extension penalty [4]
        #          -R INT    stop searching when there are >INT equally best hits [30]
        #          -q INT    quality threshold for read trimming down to 35bp [0]
        #          -f FILE   file to write output to instead of stdout
        #          -B INT    length of barcode
        #          -c        input sequences are in the color space
        #          -L        log-scaled gap penalty for long deletions
        #          -N        non-iterative mode: search for all n-difference hits (slooow)
        #          -I        the input is in the Illumina 1.3+ FASTQ-like format
        #          -b        the input read file is in the BAM format
        #          -0        use single-end reads only (effective with -b)
        #          -1        use the 1st read in a pair (effective with -b)
        #          -2        use the 2nd read in a pair (effective with -b)
        #          -Y        filter Casava-filtered sequences

        cmd_init_local = r"""
        input=$(head -n 1 {i[txt][0]})
        outDir=$dest/mapped_hg19_v3 && [[ "$dest" != s3://* ]] && mkdir -p $outDir;
        printf "%s input= %s\n\t\t\ttmpDir= %s\n\t\t\tsplitId= %s\n" "$(date "+%T %D")" "$input" "$tmpDir" "{p[split]}"
        """
        cmd_main = r"""
        printf "%s ########## BwaAln Main ##########\n" "$(date "+%T %D")"
        tmpOut=$tmpDir/out/$id."{p[split]}"   # prefix
        {s[bwa]} aln -q 15 -t 2 {s[ref1]} -b1 "$input" > "$tmpOut.sai1"
        {s[bwa]} aln -q 15 -t 2 {s[ref1]} -b2 "$input" > "$tmpOut.sai2"
        
        # Need to copy ORIGINAL input to the next stage of bwa sampe
        echo "$(head -n 1 {i[txt][0]})" > $OUT.txt
        """

        return (cmd_init + cmd_init_local + cmd_input_s3 + cmd_main + cmd_out)

class CteamBwaSampe(Tool):
    """
    IN=$DIR/merged/$ID/$ID.sort.bam
    
    OUTDIR=$DIR/mapped_hg19_v3/$ID && mkdir -p $OUTDIR
    OUT=$OUTDIR/$ID  # output prefix

    ($BWA sampe  $REF $OUT.sai1 $OUT.sai2 $IN $IN | $SAMTOOLS view -Su - | $SAMTOOLS sort - $OUT.sort) >& $OUT.sais2aln.oe;
    """

    name = "BWA_SAMPE"
    cpu_req = 2
    mem_req = 12 * 1024 # 12 GB
    
    inputs  = ['txt']
    outputs = ['txt']

    def cmd(self,i,s,p):
        # Usage:   bwa sampe [options] <prefix> <in1.sai> <in2.sai> <in1.fq> <in2.fq>

        # Options: -a INT   maximum insert size [500]
        #          -o INT   maximum occurrences for one end [100000]
        #          -n INT   maximum hits to output for paired reads [3]
        #          -N INT   maximum hits to output for discordant pairs [10]
        #          -c FLOAT prior of chimeric rate (lower bound) [1.0e-05]
        #          -f FILE  sam file to output results to [stdout]
        #          -r STR   read group header line such as `@RG\tID:foo\tSM:bar' [null]
        #          -P       preload index into memory (for base-space reads only)
        #          -s       disable Smith-Waterman for the unmapped mate
        #          -A       disable insert size estimate (force -s)

        # Notes: 1. For SOLiD reads, <in1.fq> corresponds R3 reads and <in2.fq> to F3.
        #        2. For reads shorter than 30bp, applying a smaller -o is recommended to
        #           to get a sensible speed at the cost of pairing accuracy.

        cmd_init_local = r"""
        input=$(head -n 1 {i[txt][0]}); 
         sai1=$(head -n 2 {i[txt][0]} | tail -n 1);
         sai2=$(head -n 3 {i[txt][0]} | tail -n 1);
        outDir="$dest/bwa" && [[ "$dest" != s3://* ]] && mkdir -p $outDir;
        printf "%s input= %s\n\t\t\tsai1= %s\n\t\t\tsai2=%s\n\t\t\toutDir= %s\n" "$(date "+%T %D")" "$input" "$sai1" "$sai2" "$outDir"
        """
        cmd_main = r"""
        printf "%s ########## BwaSampe Main ##########\n" "$(date "+%T %D")"
        [[ "$sai1" == s3://* ]] && aws s3 cp "$sai1" "$tmpDir" --quiet && sai1=$tmpDir/$(basename $sai1)
        [[ "$sai2" == s3://* ]] && aws s3 cp "$sai2" "$tmpDir" --quiet && sai2=$tmpDir/$(basename $sai2)

        tmpOut=$tmpDir/out/$id."{p[split]}".sampe.bam;
        
        RG="@RG\\tID:$id\\tSM:{p[rg]}\\tPL:ILLUMINA" && echo "$RG";
        
        {s[bwa]} sampe -P -r "$RG" {s[ref1]} $sai1 $sai2 $input $input |\
        {s[samtools]} view    -Su   -                                  |\
        {s[samtools]} sort -o -@ 2 -m 3G - $tmpDir/tmpSort > $tmpOut;
        
        [[ -s "$tmpOut" ]] && {s[samtools]} index "$tmpOut";        
        """

        return (cmd_init + cmd_init_local + cmd_input_s3 + cmd_main + cmd_out)

class CteamSplitByChromosome(Tool):
    """
      As Samtool rmdup can handle one file only and it does rmdup per chromosome, 
      merge split bam files here and (re)split by chromosome, using gatk printreads.
    """

    name = "SplitByChrom"
    cpu_req = 30
    mem_req = 10 * 1024

    inputs  = ['txt']
    outputs = ['txt']

    def cmd(self,i,s,p):
        # GATK: bam file should be indexed
        # GATK: can't handle gziiped reference
        cmd_init_local = r"""
        txt=$(cat {input_files});
        for f in $txt; do
           if [[ "$f" == s3://* ]] ; then
             aws s3 cp "$f" $tmpDir;
             [[ "$f" == *.bam ]] && echo "$tmpDir/$( basename $f )" >> $tmpDir/input.list;
           else
             [[ "$f" == *.bam ]] && echo "$f"                       >> $tmpDir/input.list;
           fi
        done
        outDir="$dest/PerChrom" && [[ "$dest" != s3://* ]] && mkdir -p $outDir;
        """
        
        cmd_main = r"""
        printf "%s ########## SplitByChrom Main ##########\n" "$(date "+%T %D")"
               
        for chr in $(seq 1 22) X Y MT;
        do 
           tmpOut=$tmpDir/out/$id.$chr.bam;

           {s[java]} -Djava.io.tmpdir=$tmpDir -Xmx8G \
           -jar {s[gatk]} -R {s[ref2]} -T PrintReads -et NO_ET -K {s[no_et_key]} \
           -L "$chr"                    \
           -o "$tmpOut"                 \
           -nct 2                       \
           -I $tmpDir/input.list > $tmpDir/gatk.$chr.out &

           echo "$!" >> $tmpDir/pids;
        done

        for job in $(cat $tmpDir/pids); do wait $job; done
        """

        return (cmd_init + cmd_init_local + cmd_main + cmd_out),{'input_files': _list2input(i['txt']," ")}

class CteamRmDup_BuildIndex(Tool):
    """
     IN=$DIR/mapped_hg19_v3/$ID/$ID.sort.bam
    OUT=$DIR/mapped_hg19_v3/$ID/$ID.dedup.bam

    $SAMTOOLS rmdup $IN $OUT
    $SAMTOOLS index     $OUT
    """

    name = "RmDup-Index"
    cpu_req = 1
    mem_req = 2 * 1024 # 2 GB

    inputs  = ['txt']
    outputs = ['txt']

    def cmd(self,i,s,p):
        cmd_init_local = r"""
        input=$( cat {i[txt][0]} | grep "\.{p[chrom]}\.bam" )
        index=$( cat {i[txt][0]} | grep "\.{p[chrom]}\.bai" )
        outDir="$dest/rmdup" && [[ "$dest" != s3://* ]] && mkdir -p $outDir;
        printf "%s input= %s\n\t\t\tchrom= %s\n\t\t\toutDir= %s\n" "$(date "+%T %D")" "$input" "{p[chrom]}" "$outDir"
        """
        cmd_main = r"""
        printf "%s ########## RemoveDuplicate Main ##########\n" "$(date "+%T %D")"
        [[ $index == s3://* ]] && aws s3 cp $index $tmpDir --quiet

        tmpOut=$tmpDir/out/$id.{p[chrom]}.rmdup.bam;
        
        [[ -s $input  ]] && {s[samtools]} rmdup $input $tmpOut;        
        [[ -s $tmpOut ]] && {s[samtools]} index        $tmpOut;
        """
        return (cmd_init + cmd_init_local + cmd_input_s3 + cmd_main + cmd_out)

class CteamRealignTarget(Tool):
    """
        IN=$DIR/mapped_hg19_v3/$ID/$ID.dedup.bam
    OUTDIR=$DIR/realign/$ID && mkdir -p $OUTDIR
    
    java -Xmx4g -jar $GATK -T RealignerTargetCreator -R $REF2 -L $chr -I $IN -o $OUTDIR/$chr.intervals     
    """
    name = "RealignTarget"
    cpu_req = 4
    mem_req = 8 * 1024

    inputs  = ['txt']
    outputs = ['txt']

    def cmd(self,i,s,p):
        cmd_init_local = r"""
        input=$( head -n 1 {i[txt][0]} )
        index=$input".bai"
        outDir="$dest/realign" && [[ "$dest" != s3://* ]] && mkdir -p $outDir;
        printf "%s input= %s\n\t\t\tchrom= %s\n\t\t\toutDir= %s\n" "$(date "+%T %D")" "$input" "{p[chrom]}" "$outDir"
        """
        cmd_main = r"""
        printf "%s ########## RealignTarget Main ##########\n" "$(date "+%T %D")"
        [[ $index == s3://* ]] && aws s3 cp $index $tmpDir --quiet

        tmpOut=$tmpDir/out/$id.{p[chrom]}.intervals;
        
        {s[java]} -Djava.io.tmpdir=$tmpDir -Xmx8G \
        -jar {s[gatk]} -R {s[ref2]} -T RealignerTargetCreator  -et NO_ET -K {s[no_et_key]} \
        -L {p[chrom]} \
        -I $input     \
        -nt 4         \
        -o $tmpOut
        
        echo $( head -n 1 {i[txt][0]} ) >  $OUT.txt;  # need to pass ORIGINAL input to the next stage
        """
        return (cmd_init + cmd_init_local + cmd_input_s3 + cmd_main + cmd_out)

class CteamIndelRealigner(Tool):
    """
    IN=$DIR/mapped_hg19_v3/$ID/$ID.dedup.bam
    
    OUTDIR=$DIR/realign/$ID && mkdir -p $OUTDIR
    OUT=$OUTDIR/$chr.realign.bam

    java -Xmx4g -jar $GATK -T IndelRealigner -R $REF2 -L $chr -I $IN -targetIntervals $OUTDIR/$chr.intervals -o $OUT
    """
    name = "IndelRealign"
    cpu_req = 1
    mem_req = 6 * 1024 

    inputs  = ['txt']
    outputs = ['txt']

    def cmd(self,i,s,p):
        cmd_init_local = r"""
        input=$(   head -n 1 {i[txt][0]})   # first  line = input
        interval=$(tail -n 1 {i[txt][0]})   # second line = interval
        index=$input".bai"
        outDir="$dest/realign" && [[ "$dest" != s3://* ]] && mkdir -p $outDir; 
        printf "%s input= %s\n\t\t\tchrom= %s\n\t\t\toutDir= %s\n" "$(date "+%T %D")" "$input" "{p[chrom]}" "$outDir"
        """

        cmd_main = r"""
        printf "%s ########## IndelRealigner Main ##########\n" "$(date "+%T %D")"
        [[ $interval == s3://* ]] && aws s3 cp $interval $tmpDir --quiet && interval=$tmpDir/$( basename $interval )
        [[ $index    == s3://* ]] && aws s3 cp $index    $tmpDir --quiet

        tmpOut=$tmpDir/out/$id.{p[chrom]}.realign.bam;
        
        {s[java]} -Djava.io.tmpdir=$tmpDir -Xmx6G \
        -jar {s[gatk]} -R {s[ref2]} -et NO_ET -K {s[no_et_key]} \
        -T IndelRealigner \
        -L {p[chrom]} \
        -I $input     \
        -o $tmpOut    \
        -targetIntervals $interval
        """
        return (cmd_init + cmd_init_local + cmd_input_s3 + cmd_main + cmd_out)

class CteamUnifiedGenotyper(Tool):
    """
     IN=$DIR/realign/$ID/$chr.realign.bam
    OUT=$DIR/realign/$ID/$chr.realign.vcf

    java -Xmx4g -jar $GATK -T UnifiedGenotyper -R $REF2 -L $chr -I $IN -D $DBSNP -o $OUT
    -dcov 600 -glm SNP -out_mode EMIT_ALL_SITES -stand_call_conf 5.0 -stand_emit_conf 5.0
    -inputPrior 0.0010 -inputPrior 0.4995 -A GCContent -A BaseCounts 
    """
    name = "UnifiedGenotyper"
    cpu_req = 4
    mem_req = 20 * 1024

    inputs  = ['txt']
    outputs = ['txt']

    def cmd(self,i,s,p):
        cmd_init_local = r"""
        input=$( cat {i[txt][0]} | grep bam );
        index=$( cat {i[txt][0]} | grep bai );
        outDir="$dest/ug" && [[ "$dest" != s3://* ]] && mkdir -p $outDir; 
        printf "%s input= %s\n\t\t\tchrom= %s\n\t\t\toutDir= %s\n" "$(date "+%T %D")" "$input" "{p[chrom]}" "$outDir"
        """
        cmd_main = r"""
        printf "%s ########## UnifiedGenotyper Main ##########\n" "$(date "+%T %D")"
        [[ $index    == s3://* ]] && aws s3 cp $index    $tmpDir --quiet

        tmpOut=$tmpDir/out/$id.{p[chrom]}.vcf.gz
        
        {s[java]} -Djava.io.tmpdir=$tmpDir -Xmx20g \
        -jar {s[gatk]} -R {s[ref2]} -et NO_ET -K {s[no_et_key]} \
        -T UnifiedGenotyper \
        -L {p[chrom]} \
        -I $input \
        -D {s[dbsnp_vcf]} \
        -dcov 600 \
        -glm SNP  \
        -out_mode EMIT_ALL_SITES \
        -stand_call_conf 5.0 -stand_emit_conf 5.0 \
        -inputPrior 0.0010 -inputPrior 0.4995 \
        -A GCContent -A BaseCounts \
        -nt 2 -nct 2 \
        -o $tmpOut
        """
        return (cmd_init + cmd_init_local + cmd_input_s3 + cmd_main + cmd_out)

class CteamVariantFiltration(Tool):
    """
     IN=$DIR/realign/$ID/$chr.realign.vcf
    OUT=$DIR/realign/$ID/$chr.realign.filt.vcf

    (java -Xmx8g -jar $GATK -T VariantFiltration -R $REF2 -L $chr -V $IN -mask:BED $MASK -o $OUT
     -filter \"QUAL\<20.0\" -filter \"QD\<2.0\" -filter \"MQ\<20.0\" -filter \"DP\<5\" -filter \"DP\>50\"
     -filterName LowQual -filterName LowQD -filterName LowMQ -filterName LowDP -filterName HighDP -maskName NotInMask) >& $ID.realign.$chr.filt.vcf.oe
    """
    name = "VariantFiltration"
    cpu_req = 1
    mem_req = 8 * 1024

    inputs  = ['txt']
    outputs = ['txt']

    def cmd(self,i,s,p):
        cmd_input = r"""
        input=$(head -n 1 {i[txt][0]}) && echo "input=$input";
        """
        cmd_main = r"""
        printf "%s ########## VariantFiltration Main ##########\n" "$(date "+%T %D")"
        outDir="{s[outputDir]}/realign/filtered" && mkdir -p $outDir;
        output=$outDir/$id."{p[chrom]}".filt.vcf;
        
        {s[java]} -Djava.io.tmpdir=$tmpDir -Xmx8g \
        -jar {s[gatk]} -R {s[ref2]} -T VariantFiltration  -et NO_ET -K {s[no_et_key]} \
        -L {p[chrom]} \
        -V $input \
        -mask:BED {s[mask]} -maskName NotInMask \
        -filter \"QUAL\<20.0\" -filter \"QD\<2.0\" -filter \"MQ\<20.0\" -filter \"DP\<5\" -filter \"DP\>50\" \
        -filterName LowQual -filterName LowQD -filterName LowMQ -filterName LowDP -filterName HighDP \
        -o $output;
        
        echo $output > $OUT.txt;
        /bin/rm -rf $tmpDir;
        # echo "Remove input file for disk space" && /bin/rm -rf $input && touch $input;
        """
        return (cmd_init + cmd_input + cmd_input_s3 + cmd_main)

################################################################################
#    NOT USED IN PIPELINE
################################################################################
# class OLD_CteamLinkMerge(Tool):
#     """
#     ln -s $OUTDIR/$ID.sort.rg.dedup_RG.bam  $OUTDIR3/$RG.bam
#     print $ID.hdr \"\@\RG\tID:$RG\tSM:$RG\tPL:ILLUMINA\n\"
#     cat nop.header.txt | $SAMTOOLS view -Sbo $ID.nop.bam -

#     IN=$OUTDIR/$SAMP_ID
#     OUT=$IN.rg.final.bam

#     ($SAMTOOLS merge -rh $ID.hdr $OUT $IN.bam $ID.nop.bam ) >& $IN.merge.oe
#     """

#     name = "LinkMerge"
#     cpu_req = 2
#     mem_req = 1 * 1024 # 1 GB
#     time_req = 1 * 60

#     inputs  = ['bam']
#     ouptuts = ['final.bam']

#     def cmd(self,i,s,p):
#         return r"""
#         set -e -o pipefail;
        
#         {s[samtools]} merge -rh {p[hdr]} $OUT.final.bam {i['bam'][0]} $OUT.nop.bam
#         """

# class OLD_CteamSortBam(Tool):
#     """
#     DIR=/files/Genetics/1000Genomes/cteam/RAW_DATA/hard_drives/disk#/
#     IN=$DIR/$ID/Assembly/genome/bam/$ID.bam

#     OUTDIR=$DIR/sorted_bams/$ID && mkdir -p $OUTDIR
#     OUT=$OUTDIR/$ID.sort   # output prefix

#     $SAMTOOLS sort -n $IN $OUT
#     """
    
#     name = "SortBAM"  # gridengine DRM does NOT allow names starting with numbers
#     cpu_req  = 30
#     mem_req  = 200 * 1024  # 200 GB
#     time_req = 2 * 60

#     inputs  = ['bam']
#     outputs = ['txt']  # txt file that has the location of actual data file

#     def cmd(self,i,s,p):
#         # -n: Sort by read names rather than by chromosomal coordinates.
#         # -o: Output the final alignment to the standard output.
#         # -m INT: Approximately the maximum required memory per thread. [500000000] 
#         # -l: compression level - works?
#         # r3.8xlarge has 32 vCPUs and 240 GB mem: max 6G per thread?
#         return r"""
#         id={p[rg]};
        
#         outDir={s[outputDir]}/sorted_bams/$id && mkdir -p $outDir;
#         out=$outDir/$id.sort; # prefix;
        
#         {s[samtools]} sort -o -n -@ 16 -m 12G -l 0 {i[bam][0]} {s[scratch]}/$id.sort > $out.bam;
        
#         echo $out.bam > $OUT.txt;
#         """

# class OLD_CteamTrimReadGroup(Tool):
#     """
#     IN=$DIR/sorted_bams/$ID/$ID.sort.bam

#     OUTDIR=$DIR/merged/$ID && mkdir -p $OUTDIR
#     OUT=$OUTDIR/$ID.sort.bam

#     $REMOVE_TAG_MAPPING $IN /dev/stdout |
#     $ADD_RG_CTEAM /dev/stdin /dev/stdout $RG_ID |
#     $MERGE_TRIM_READS_BAM --keepOrig -o $OUT --log $OUT.trim.log /dev/stdin
#     """

#     name = "2-TrimReadGroup"
#     cpu_req = 2
#     mem_req = 4 * 1024 # 4GB
#     time_req = 2 * 60

#     inputs  = ['txt']
#     outputs = ['txt']
    
#     def cmd(self,i,s,p):
#         return r"""
#         input=$(head -n 1 {i[txt][0]});
        
#         outDir={s[outputDir]}/merged/$id && mkdir -p $outDir;
#         output=$outDir/$id.sort.bam;  # actual bam file name, not prefix;
#         outlog=$outDir/$id.trim.log;
        
#         {s[removeTagMapping]}  $input /dev/stdout             |
#         {s[addRgCteam]}        /dev/stdin /dev/stdout {p[rg]} |
#         {s[mergeTrimReadsBam]}  --keepOrig -o $output --log $outlog /dev/stdin;
        
#         echo $output > $OUT.txt;               
#         """

