import argparse
import os
import sys
import os.path
import subprocess
import shutil
import csv
os.system("module load bamtools/2.3.0")
os.system("module load augustus")
os.system('module load bowtie2')
os.system('module load tophat/2.0.10')
os.system('module load samtools')

__version__ = 'v0.0.0'
__date__ = '2017-04-04'
__updated__ = '2017-04-04'

class CLIError(Exception):
    '''Generic exception to raise and log different fatal errors.'''
    def __init__(self, msg):
        super(CLIError).__init__(type(self))
        self.msg = 'E: %s' % msg
    def __str__(self):
        return self.msg
    def __unicode__(self):
        return self.msg

def test_dir(option, opt_str, value, parser):
    if os.path.exists(value):
        setattr(parser.values, option.dest, value)
    else:
        print("directory of reads cannot be found")
        sys.exit()

def map_reads_to_exon_junctions(genome, genome_base, R1_numbered, R2_numbered):
    exex = '%s_exex1' % genome_base
    tophat2_out1 = 'tophat2_toMaskedGenome'
    bowtie_out = 'bowtie.sam'
    bowtie_view = 'bowtie.F.sam'
    map = 'map.psl'
    global_bowtie = 'bowtie.global.sam'
    input_bam = 'accepted_hits.bam'
    output_bam = 'accepted_hits.noN.bam'
    header = 'header.txt'
    sam_and_header = 'bowtie.global.h.sam'
    bam_and_header = 'bowtie.global.h.bam'
    both_bams = 'both.samtools.bam'
    both_sams = 'both.2.samtools.s'
    #subprocess.Popen(['cp', R1_numbered, tophat2_out1]).wait()
    #subprocess.Popen(['cp', R2_numbered, tophat2_out1]).wait()
    os.chdir(tophat2_out1)
    with open(bowtie_view, 'w') as handle, open(global_bowtie, 'w') as gbow, open(output_bam, 'w') as obam, open(sam_and_header, 'w') as shead, open(bam_and_header, 'w') as bhead, open(both_bams, 'w') as bbams, open(both_sams, 'w') as bsams:
        subprocess.Popen(['bowtie2', exex, '-1', R1_numbered, '-2', R2_numbered, '-S', bowtie_out]).wait()
        subprocess.Popen(['samtools', 'view', '-S', '-F', '4', bowtie_out], stdout=handle).wait()
        subprocess.Popen(['/packages/augustus/3.0.1/scripts/samMap.pl', bowtie_view, map, '100'], stdout=gbow).wait()
        subprocess.Popen(['bamtools', 'filter', '-in', input_bam, '-out', output_bam, '-script', '/packages/augustus/3.0.1/auxprogs/auxBamFilters/operation_N_filter.txt']).wait()
        subprocess.Popen(['cat', header, global_bowtie], stdout=shead).wait()
        subprocess.Popen(['samtools', 'view', '-bs', '-o', bam_and_header, sam_and_header]).wait()
        subprocess.Popen(['samtools', 'merge', both_bams, bam_and_header, output_bam]).wait()
        subprocess.Popen(['samtools', 'sort', '-n', both_bams, both_sams]).wait()
        #subprocess.Popen(['/packages/augustus/3.0.1/auxprogs/filterBam/bin/filterBam', '--uniq', '--paired', '--in',
        
def make_exon_exon_boundaries(aug1_predictions, genome_base, genome_masked_name):
    intron_locations = 'intron_locations.txt'
    tophat2_out1 = 'tophat2_toMaskedGenome'
    aug_prelim = 'aug.prelim.gff'
    aug1_introns = 'aug1.introns.gff'
    introns = 'introns.lst'
    exex = 'exex.fa'
    map = 'map.psl'
    exex_reference = '%s_exex1' % genome_base
    subprocess.Popen(['cp', genome_masked_name, tophat2_out1]).wait()
    os.chdir(tophat2_out1)
    with open(aug_prelim, 'w') as aug1_prelim, open(aug1_introns, 'w') as aug_introns, open(intron_locations, 'w') as intron_file:
        s1 = subprocess.Popen(['cat', aug1_predictions], stdout=subprocess.PIPE, stderr=sys.stderr)
        #s2 = subprocess.Popen(['tee', aug_prelim], stdin=s1.stdout, stdout=subprocess.PIPE).wait()
        s3 = subprocess.Popen(['grep', '-P', "\tintron\t"], stdin=s1.stdout, stdout=aug_introns).wait()
        subprocess.Popen(['cp', aug1_introns, aug_prelim])
        with open(aug1_introns, 'r') as handle:
            for inline in handle:
                infields = inline.strip().split('\t')
                contig_name = str(infields[0])
                intron_start = str(infields[3])
                intron_end = str(infields[4])
                dash = '-'
                intron_locs = (contig_name + ':' + intron_start + dash + intron_end + '\n')
                intron_file.write(intron_locs)
        with open(introns, 'w') as final, open(exex, 'w') as ex, open(map, 'w') as mp:
            subprocess.Popen(['sort', '-u', intron_locations], stdout=final).wait()
            subprocess.Popen(['/scratch/bin/augustus.2.7/scripts/intron2exex.pl', '--introns=%s' % introns, '--seq=%s' % genome_masked_name, '--exex=%s' % exex, '--map=%s' % map, '--flank=100']).wait()
            subprocess.Popen(['bowtie2-build', exex, exex_reference]).wait()

def run_augustus_with_hints(aug_species, genome, genome_base, aug1_predictions):
    tophat2_out1 = 'tophat2_toMaskedGenome'
    extrinsic_file = 'extrinsic.cfg'
    hints = 'hints.gff'
    aug1_predictions = '%s_aug1.out' % genome_base
    subprocess.Popen(['cp', genome, tophat2_out1]).wait()
    os.chdir(tophat2_out1)
    with open(aug1_predictions, 'w') as handle:
        subprocess.Popen(['cp', '/packages/augustus/3.0.1/config/extrinsic/extrinsic.M.RM.E.W.cfg', extrinsic_file]).wait()
        subprocess.Popen(['augustus', '--species=%s' % aug_species, '--extrinsicCfgFile=%s' % extrinsic_file, '--alternatives-from-evidence=true', '--hintsfile=%s' % hints, '--allow_hinted_splicesites=atac', '--introns=on', '--genemodel=complete', genome], stdout=handle).wait()

def make_augustus_hints(bowtie_index, R1_numbered, R2_numbered, header):
    initial_tophat_run = 'accepted_hits.bam'
    sorted_bam = 'accepted_hits.s'
    sorted_fbam = 'accepted_hits.sf'
    new_name_sorted_bam = 'accepted_hits.s.bam'
    new_sorted_fbam = 'accepted_hits.sf.bam'
    tophat2_out1 = 'tophat2_toMaskedGenome'
    header = 'header.txt'
    hints_input = 'both.ssf'
    hints_new_name = 'both.ssf.bam'
    hints = 'hints.gff'
    subprocess.Popen(['mkdir', tophat2_out1]).wait()
    #with open(initial_tophat_run, 'w') as run1:
    subprocess.Popen(['tophat2', '-o', tophat2_out1, '-i', '35', '-I', '1000', '-p', '8', bowtie_index, R1_numbered, R2_numbered, '-r', '50', '--mate-std-dev', '100']).wait()
    os.chdir(tophat2_out1)
    with open(header, 'w') as head:
        sort_bam = subprocess.Popen(['samtools', 'sort', '-n', initial_tophat_run, sorted_bam]).wait()
        subprocess.Popen(['/scratch/bin/augustus.2.7/auxprogs/filterBam/bin/filterBam', '--uniq', '--paired', '--in', new_name_sorted_bam, '--out', new_sorted_fbam]).wait()
        view_bam = subprocess.Popen(['samtools', 'view', '-H', new_sorted_fbam], stdout=head).wait()
        subprocess.Popen(['samtools', 'sort', new_sorted_fbam, hints_input]).wait()
        subprocess.Popen(['/scratch/bin/augustus.2.7/auxprogs/bam2hints/bam2hints', '--intronsonly', '--in', hints_new_name, '--out', hints]).wait()
def building_reference(genome, masked_genome, genome_masked_name, bowtie_index):
    shutil.move(masked_genome, genome_masked_name)
    subprocess.Popen(['bowtie2-build', genome_masked_name, bowtie_index]).wait()

def mark_repeats(genome, genome_base, masked_genome):

    lfreq_table = "%s_lfreq.txt" % genome_base
    repeat_out = "%s.out" % 'repeat_scout'
    repeat_filter1 = 'repeat_scout_filter1.fasta'
    repeat_filter2 = 'repeat_scout_filter2_10.fasta'
    genome_out = '%s.out' % genome
    fungal_repbase = 'fungal_repbase.fa'
    all_repeat = 'all_repeats.fa'
    masked_genome_dir = 'masked_genome'
    all_genomedot = '%s.out' % genome
    removal_of_ori_repeat_masker = '%s.ori.out' % genome
    removal_of_more_repeat_files_to_rerun = '%s.tbl' % genome
    more_removal = '%s.cat.gz' % genome
    if genome.endswith('.fasta'):
        with open(lfreq_table, 'w') as lfreq_table_file, open(repeat_out, 'w') as repeat_final, open(repeat_filter1, 'w') as filter1, open (repeat_filter2, 'w') as filter2, open (fungal_repbase, 'w') as fungal_rep, open (all_repeat, 'w') as all_repeats:
            subprocess.Popen(['/scratch/bin/RepeatScout-1/build_lmer_table', '-sequence', genome, '-freq', lfreq_table], stdout=sys.stdout, stderr=sys.stderr).wait()
            print('create lfreq table: DONE')
            subprocess.Popen(['/scratch/bin/RepeatScout-1/RepeatScout', '-sequence', genome, '-output', repeat_out, '-freq', lfreq_table]).wait()
            print('run repeat scout: DONE')
            #filter-stage from repeat scout, filters any sequence deemed to be more than %50 low complexity
            ps = subprocess.Popen(['cat', repeat_out], stdout=subprocess.PIPE, stderr=sys.stderr)
            subprocess.Popen('/scratch/bin/RepeatScout-1/filter-stage-1.prl', stdin=ps.stdout, stdout=filter1).wait()
            print('first round of repeat filter is complete.')
            subprocess.Popen(['/packages/RepeatMasker/4.0.3/RepeatMasker', '-pa', '8', '-lib', repeat_filter1, genome, '-e', 'ncbi']).wait()
            print('first repeat masker is complete.')
            filt2 = subprocess.Popen(['cat', repeat_out], stdout=subprocess.PIPE, stderr=sys.stderr)
            subprocess.Popen((['/scratch/bin/RepeatScout-1/filter-stage-2.prl', '--cat', genome_out, '--thres', '10']),  stdin=filt2.stdout, stdout=filter2).wait()
            print('second round of repeat scout is complete.')
            subprocess.Popen((['/packages/RepeatMasker/4.0.3/util/queryRepeatDatabase.pl', '-species', 'fungi']), stdout=fungal_rep).wait()
            print('query repeat database is complete.')
            subprocess.Popen((['cat', repeat_filter2, fungal_repbase]), stdout=all_repeats).wait()
            subprocess.Popen((['rm', all_genomedot, removal_of_ori_repeat_masker, removal_of_more_repeat_files_to_rerun, more_removal, masked_genome]))
            subprocess.Popen(['/packages/RepeatMasker/4.0.3/RepeatMasker', '-pa', '8', '-lib', all_repeat, genome, '-e', 'ncbi']).wait()
            print('final repeat masker is complete.')
            subprocess.Popen((['mkdir', masked_genome_dir])).wait()
            subprocess.Popen((['mv', repeat_filter2, masked_genome_dir]))
            subprocess.Popen((['mv', all_repeat, masked_genome_dir]))
            subprocess.Popen((['mv', repeat_filter1, masked_genome_dir]))
            subprocess.Popen((['mv', repeat_out, masked_genome_dir]))
            subprocess.Popen((['mv', fungal_repbase, masked_genome_dir]))
            subprocess.Popen((['mv', lfreq_table, masked_genome_dir]))
            #shutil.move(masked_genome, masked_genome_dir)

def decompress(read1, read2, unzip_R1, unzip_R2):
    # TODO: add other compression formats
    os.system('gunzip -c %s > %s' % (read1, unzip_R1))
    print(read1 + ' unzipped')
    os.system('gunzip -c %s > %s' % (read2, unzip_R2))
    print(read2 + ' unzipped')
   
def parse_input_reads(read1, read2, unzip_R1, unzip_R2, R1_numbered, R2_numbered):
    #check if the read file extension is .gz
    extension1 = os.path.splitext(read1)[1][1:]
    extension2 = os.path.splitext(read2)[1][1:]
    if extension1 and extension2 == 'gz':
        decompress(read1, read2, unzip_R1, unzip_R2)
    else:
        unzip_R1 = read1
        unzip_R2 = read2

    def prep_reads_with_numbers(unzip_R1, unzip_R2, R1_numbered, R2_numbered):
        with open(unzip_R1, 'r') as f, open(R1_numbered, 'w') as out1:
            for inline in f:
                split_line = inline.strip().split()
                if len(split_line) > 1:
                    read_name = split_line.pop(0)
                    read_name += '-1\n'
                    out1.write(read_name)
                else:
                    out1.write(inline)
        with open(unzip_R2, 'r') as f, open(R2_numbered, 'w') as out2:
            for inline in f:
                split_line = inline.strip().split()
                if len(split_line) > 1:
                    read_name = split_line.pop(0)
                    read_name += '-2\n'
                    out2.write(read_name)
                else:
                    out2.write(inline)
    prep_reads_with_numbers(unzip_R1, unzip_R2, R1_numbered, R2_numbered)
  
def main(argv=None):
    '''Command line options.'''

    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    program_name = os.path.basename(sys.argv[0])
    program_version = 'v%s' % __version__
    program_build_date = str(__updated__)
    program_version_message = '%%(prog)s %s (%s)' % (program_version, program_build_date)
    program_license = ''' %s

    Created by TGen North on %s.
    Copyright 2017 TGen North. All rights reserved.

    Available for academic and research use only under a license
    from The Translational Genomics Research Institute (TGen)
    that is free for non-commercial use.

    Distributed on as 'AS IS' basis without warranties
    or conditions of any kind, either express or implied.

    USAGE
    ''' % ('Augustus RNASeq Pipeline Step 1 ', str(__date__))

    parser = argparse.ArgumentParser(description=program_license, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-g', '--genome', required=True, help='DNA assembly [REQUIRED]')
    parser.add_argument('-r1', '--read1', required=True,
                                help='read 1 [REQUIRED]')
    parser.add_argument('-r2', '--read2', required=True, help='read 2 [REQUIRED]')
    parser.add_argument('-s', '--species', required=True, help='species identifier from augustus'
                                                               ' --species=help [REQUIRED)')
    parser.add_argument('-o', '--out-dir', dest='outdir', help='directory to write output files. [optional]')
    parser.add_argument('-v', '--version', action='version', version=__version__)

    args = parser.parse_args()

    genome = args.genome
    read1 = args.read1
    read2 = args.read2
    aug_species = args.species
    genome_base = genome[:-len(".fasta")]
    masked_genome_dir = 'masked_genome'
    genome_masked_name = '%s.masked.fa' % genome_base
    masked_genome = '%s.masked' % genome
    bowtie_index = '%s.masked' % genome_base
    unzip_R1 = read1[:-len(".gz")]
    unzip_R2 = read2[:-len(".gz")]
    R1_base = unzip_R1[:-len(".fasta")]
    R2_base = unzip_R2[:-len(".fasta")]
    R1_numbered = '%s_numbered_1.fastq' % R1_base
    R2_numbered = '%s_numbered_2.fastq' % R2_base
    aug1_predictions = '%s_aug1.out' % genome_base
    header = 'header.txt'
#    parse_input_reads(read1, read2, unzip_R1, unzip_R2, R1_numbered, R2_numbered)
 #   mark_repeats(genome, genome_base, masked_genome)
#   building_reference(genome, masked_genome, genome_masked_name, bowtie_index)
#    make_augustus_hints(bowtie_index, R1_numbered, R2_numbered, header)
#    run_augustus_with_hints(aug_species, genome, genome_base, aug1_predictions)
#    make_exon_exon_boundaries(aug1_predictions, genome_base, genome_masked_name)
    map_reads_to_exon_junctions(genome, genome_base, R1_numbered, R2_numbered)
if __name__ == "__main__":
    main()
