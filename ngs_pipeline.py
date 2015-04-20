__author__ = 'avivro'
import StringIO
import os
import subprocess
import re
from Bio import SeqIO
import gzip
from Bio.Seq import Seq
import csv
import sys
import swalign
from os import path
import shelve
import time
#primers act as the adapters in the seqprep tool, need to see how this acts with umi's
# DEFAULT_FORWARD_PRIMER='NNNNNNNNNCAGCTCTTCGCCTTTACGCATATG' with UMI
DEFAULT_FORWARD_PRIMER='CAGCTCTTCGCCTTTACGCATATG'
# DEFAULT_REVERSE_PRIMER ='ATGAAAAGCTTAGTCATGGCG' on the sense strand
# DEFAULT_REVERSE_PRIMER ='CGCCATGACTAAGCTTTTCAT' not including the part of the restriction
DEFAULT_REVERSE_PRIMER ='GGCGCGCCATGACTAAGCTTTTCAT'
# DEFAULT_FORWARD_PRIMER='GGCGCGCCATGACTAAGCTTTTCATTGTCATGC'
# DEFAULT_REVERSE_PRIMER = 'CATATGCGTAAAGGCGAAGAGCTGCTGTGTAGATCT'

#home directory where everthing happens
DEFAULT_HOME_DIR = '/home/labs/pilpel/avivro/workspace/'
#the variants refrences in the format of <id> <sequence> in tab delimited format
DEFAULT_REFERENCE_FILE = DEFAULT_HOME_DIR + 'data/reference_variant_full_sequences.tab'

#the location of the directory with the fastq files of each sample
DEFAULT_BIN_DIR = DEFAULT_HOME_DIR  + "data/fitseq_sample_data/Project_fitseq_1anc/"
# DEFAULT_BIN_DIR = DEFAULT_HOME_DIR  + 'FitSeq/data/goodman_raw_data/Project_goodman'

#the directory of the results of the pipeline
# DEFAULT_RESULT_DIR = DEFAULT_HOME_DIR  + 'FitSeq/data/ngs_sample_data/ngs_pipeline_result'
DEFAULT_RESULT_DIR = DEFAULT_HOME_DIR  + 'data/fitseq_sample_data/ngs_pipeline_result_1anc_1804/'

#the file with the result frequency matrix
# DEFAULT_RESULT_FILE = DEFAULT_RESULT_DIR + '/Project_avivro_result.csv'
DEFAULT_RESULT_FILE = DEFAULT_RESULT_DIR + 'fitseq_sample_result.csv'

#the files that were discarded at the trim umi count stage
# DEFAULT_DISCARDED_FILE = DEFAULT_RESULT_DIR + '/Project_avivro_discarded.csv'
DEFAULT_DISCARDED_TRIM_FILE = DEFAULT_RESULT_DIR + 'fitseq_sample_discarded_trim'


#the files that were discarded at the variant count stage
# DEFAULT_DISCARDED_FILE = DEFAULT_RESULT_DIR + '/Project_avivro_discarded.csv'
DEFAULT_DISCARDED_VARIANT_FILE = DEFAULT_RESULT_DIR + 'fitseq_sample_discarded_variant'

#summary file to summarize stats. line 1:merged reads 2:trimmed reads 3:number of unique fragments found 4:number of unique umis for all fragments
DEFAULT_SUMMARY_FILE= DEFAULT_RESULT_DIR + 'fitseq_sample_summary.txt'



DEFAULT_UMI_LENGTH = 9
DEFAULT_R_SHIFT_LENGTH = 4
DEFAULT_DESIGN_MAX_LENGTH = 94
DEFAULT_DESIGN_MIN_LENGTH = 91
DEFAULT_MATCH_SCORE= 2
DEFAULT_MISMATCH_SCORE = -1
DEFAULT_SCORING_MATRIX = swalign.NucleotideScoringMatrix(DEFAULT_MATCH_SCORE, DEFAULT_MISMATCH_SCORE)
DEFAULT_SW_ALIGN = swalign.LocalAlignment(DEFAULT_SCORING_MATRIX)


def load_reference_dict(reference_location = DEFAULT_REFERENCE_FILE):
    #method for loading the refence file into a variant counting dictionary
    #recieves the reference file path

    #opens the refence file as a seqio module
    handle = open(reference_location, "rb")
    reference_dict = {}
    for record in SeqIO.parse(handle, "tab") :
        #sequences are fed into the reference dict where the id is the key and the value is a list with the id and the frequency
        reference_dict[str(record.seq.upper().reverse_complement())] = [record.id,0]
    handle.close()
    return reference_dict

def trim_to_restricition_site(seq):
    #method to trim one sequence to the canonical restriction sites and cut it to the variable regions
    #recieves the sequence to trim

    #left restriction site regular expression
    re1 = re.compile('(?<=CATATG)(.+)', flags=re.IGNORECASE)
    #right restriction site regular expression
    re2 = re.compile('(.+?)(?=GGCGCGCC)', flags=re.IGNORECASE)
    #both sides restriction site regular expression
    re3 = re.compile('(?<=CATATG)(.+?)(?=GGCGCGCC)', flags=re.IGNORECASE)

    #a list showing if the right or left sites where trimmed, the middle value is the sequnece
    trim = [None, None, None]

    # Check if both restriction sites are there

    hit1 = re3.search(seq)

    if hit1:
        trim[0] = 1
        trim[1] = str(hit1.group(1))
        trim[2] = 1
        return trim
    elif not hit1:
        # Check if only CATATG is present
        hit2 = re1.search(seq)
        if hit2:
            trim[0] = 1
            trim[1] = str(hit2.group(1))
            trim[2] = None
            return trim
        elif not hit2:
            # Check if only GGCGCGCC is present
                hit3 = re2.search(seq)
                if hit3:
                    trim[0] = None
                    trim[1] = str(hit3.group(1))
                    trim[2] = 1
                    return trim
                else:
                    #if none are present
                    trim[0] = None
                    trim[1] = str(seq)
                    trim[2] = None
                    return trim


def trim_all_merged_sequences(bin, merged_seq_file_location):
    #method to take all the merged sequences and trims them to the variable regions
    #recieves the bin name, and the merged sequences file path

    #open the merged sequences which are outputed as gzips
    merged_seq_file = gzip.open(merged_seq_file_location,'rb')
    #create new trimmed file location
    trimmed_merged_seq_file_location = merged_seq_file_location[:-7] + 'T.fq'
    #open file as csv writer
    trimmed_merged_seq_csv = csv.writer(open(trimmed_merged_seq_file_location, 'wb'))


    for record in SeqIO.parse(merged_seq_file, "fastq"):
        #open as fastq send to trim method and write to csv so we also have the bin and record name for discarded records
        trimmed_seq = trim_to_restricition_site(str(record.seq))
        trimmed_merged_seq_csv.writerow([bin, record.name, trimmed_seq[1]])
    return trimmed_merged_seq_file_location



def run_seqprep(file_set,location,file_number,bin,all_reads,merged_reads,
                adapter_1=DEFAULT_FORWARD_PRIMER, adapter_2=DEFAULT_REVERSE_PRIMER):
        #this method runs the seqprep tool on a set of files from a certain sample and results in a merged and
        #adapter trimmed file
        #recieves the a tuple with the locations of the forward and reverse read input files, the file number, the bin number
        #the forward and revers adapter

        #create output directory path
        output_dir = '/'.join([location,'out/'])
        #if the the output directory doesn't exist than create one
        if not os.path.exists(output_dir):
          os.makedirs(output_dir)

        #names for the forward read output file, the reverse read output file, and the merged read output file
        f_prefix = ''.join([output_dir,'output.',os.path.split(file_set['F'])[1]])
        r_prefix = ''.join([output_dir,'output.',os.path.split(file_set['R'])[1]])
        merged_result_file = output_dir+bin + '.file_' + str(file_number)+'.M.fq.gz'
        summary_file = output_dir+bin + '.file_' + str(file_number)+'.summary.txt'

        #call the seqprep tool using the subprocesses library which integratges with linux os
        # pipe = subprocess.call([
        pipe = subprocess.Popen([
        'SeqPrep',
        '-f', file_set['F'],
        '-r', file_set['R'],
        '-1', f_prefix[:-6]+'.fq.gz',
        '-2', r_prefix[:-6]+'.fq.gz',
        '-3', f_prefix[:-6]+'.disc.fq.gz',
        '-4', r_prefix[:-6]+'.disc.fq.gz',
        '-s', merged_result_file,
        # '-A', adapter_1, '-B', adapter_2,
        '-X', '1', '-g', '-L', '5']
        # )
        , stdout=subprocess.PIPE,stderr=subprocess.PIPE, )
        output = pipe.communicate()[1]
        output_file = StringIO.StringIO(output)
        for line in output_file:
            if line.find('Pairs Processed:')> -1:
                # print 'in pairs processed!!!'
                # print line.replace('Pairs Processed:','').strip()
                all_reads =  all_reads + int(line.replace('Pairs Processed:','').strip())
            if line.find('Pairs Merged:')>-1:
                # print 'in pairs processed!!!'
                # print line.replace('Pairs Merged:','').strip()
                merged_reads = merged_reads + int(line.replace('Pairs Merged:','').strip())
        # '-X', '1', '-g', '-L', '5'])

        subprocess.call(['gzip', '-d','-f', merged_result_file])
        merged_result_file = merged_result_file[:-3]
        return merged_result_file, all_reads, merged_reads


def count_variants(bin,fragment_umi_record_dict,reference_dictionary, discarded_variant_csv, discarded_variant_fasta):
    #a method that receives the trimmed merged sequences of one sequencing result file and counts the variant frequency
    #returns a dictionary with the variant names and a count of them in this file
    # the method saves  all the sequences that don't map to a discarded sequences csv
    #recieves the trimmed merged file path, the reference sequence dictionary and a csv writer object to write
    # the discarded sequences
    # TODO maybe create a method to summarize the discarded and mapped sequences
    # TODO create an output varaint frequencies to file mode



    for fragment,umi_dict in fragment_umi_record_dict.iteritems():
        #parse bin, sequnce name in sequencing output file and sequence from the csv
        umi_record_dict = umi_dict['umi_dict']
        umi_set = set(umi_record_dict.keys())
        #if the reference dictionary has a variant with this sequence, count it
        reverse = str(Seq(fragment).reverse_complement())
        if reference_dictionary.has_key(str(fragment)):
            reference_dictionary[str(fragment)][1] = reference_dictionary[str(fragment)][1] + len(umi_set)
        #if not check reverse complement
        elif reference_dictionary.has_key(reverse):
            ('has  reverse key')
            reference_dictionary[reverse][1] = reference_dictionary[reverse][1] + len(umi_set)
        #else write the bin, name and sequence to discarded csv
        else:
            for umi,record_list in umi_record_dict.items():
                for record in record_list:
                    discarded_variant_csv.writerow([bin,record.name,fragment,umi])
                    SeqIO.write(record,discarded_variant_fasta,'fastq')
            # print '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    #take the sequnce to id:frequency dictionary and transform it to a id:frequency dictionary
    variant_count = {variant:frequency for variant, frequency in reference_dictionary.values()}
    return variant_count


def bin_variant_frequency(bin_name,location,ref_file,discarded_trim_file,discarded_variant_file,summary_file_location):
    #method that goes over the contents of a directory finds F R file pairs and sends them to be merged, trimemed and
    #counted. this method will return a dictionary of variant ids and numbers of counts for a bin
    #receives the anme of the bin, the directory of the sequncing fils. the list of files, the reference dictionary
    #the discarded sequenc csv writer object, and the write to file flag which isn't used
    # TODO output variant frequency to file option

    #create a list for the files that will be merged which is the size of all the R F pairs but starting from 1
    files = []
    files = os.listdir(location)


    reference_dictionary = load_reference_dict(ref_file)
    result_dir = path.dirname(discarded_trim_file)
    if not os.path.exists(result_dir):
                os.makedirs(result_dir)
    discarded_trim_csv  = csv.writer(open(discarded_trim_file + '.csv','wb'))
    discarded_trim_fasta = open(discarded_trim_file + '.fq','wb')
    discarded_variant_csv = csv.writer(open(discarded_variant_file + '.csv','wb'))
    discarded_variant_fasta = open(discarded_variant_file + '.fq','wb')

    mergable_files =[{} for _ in xrange(len(files)/2 + 1)]

    #go over all the files
    for file in files:
        if file == 'out':
            continue
        #split the file in order to parse it
        split_file = file.split('_')
        #get the file number
        file_num = split_file[-1].split('.')[0]

        #create the dictionary for the f and r files if it doesn't exist
        # if mergable_files[int(file_num)] is None:
        #     mergable_files[int(file_num)] = {}
        #parse the read and use it to set the direction value
        file_read = split_file[-2]
        if file_read == 'R1':
            file_direction = 'F'
        if file_read == 'R3':
            file_direction = 'R'
        #put the file directory concatenated to the name in the current file number slot in the F or R slot
        cur_dict =  mergable_files[int(file_num)]
        cur_dict[file_direction] = '/'.join([location,file])

    #create list of variant frequency the length of the number of mergable files plus one (so it starts from 1)
    variant_frequencies = [{} for _ in xrange(len(mergable_files)+1)]

    all_reads = 0
    merged_reads = 0
    trimmed_reads = 0
    # if the file set doesn't have any files than don't loop
    for index,file_set in  enumerate(mergable_files):
        if len(file_set)<1:
            continue
        #take the mergable files list of sublists and run seqprep on it
        merged_result_file,all_reads, merged_reads =  run_seqprep(file_set,location,index,bin_name,all_reads, merged_reads)
        #trim merged file results from seqprep
        # trimmed_merged_result_file = trim_all_merged_sequences(bin_name,merged_result_file)
        #trim merged file results from seqprep
        fragment_umi_record_dict,trimmed_reads =  get_umi_counts(bin_name,merged_result_file,discarded_trim_csv,
                                                                 discarded_trim_fasta,trimmed_reads)
        #take the trimmed frequencines and count the variants in each file, put the resulting dictionary in the index of the file number
        variant_frequencies[index] = count_variants(bin_name,fragment_umi_record_dict,reference_dictionary,
                                                    discarded_variant_csv,discarded_variant_fasta)

    #create empty dictionary to aggregate the varaint counts to one dictionary
    bin_var_freq_dict ={}
    #go over variant count dictionary and merge all counts to one dictionary using the keys from any one of them
    for var_freq_dict in variant_frequencies:
        #if this dictionary exists
        if len(var_freq_dict.items()) > 0:
            for variant,frequency in var_freq_dict.items():
                #if the variant already exists in the result dict, add the value from this dict to it
                if bin_var_freq_dict.has_key(variant):
                    bin_var_freq_dict[variant] = bin_var_freq_dict[variant] + frequency
                #if not create variant in the result dict
                else:
                    bin_var_freq_dict[variant]  = frequency
    #print the sum of all bin counts
    found_variants = 0
    all_freqs = 0
    for variant,freq in bin_var_freq_dict.items():
        if int(freq) > 0:
            found_variants = found_variants +1
        all_freqs = all_freqs + freq
    summary_csv = csv.writer(open(summary_file_location, 'wb'))
    summary_csv.writerow(['',bin_name])
    summary_csv.writerow(['all_reads',all_reads])
    summary_csv.writerow(['merged_reads',merged_reads])
    summary_csv.writerow(['trimmed_reads',trimmed_reads])
    summary_csv.writerow(['found_variants',found_variants])
    summary_csv.writerow(['all_freqs',all_freqs])
    return bin_var_freq_dict



def go_over_bins(bin_dir  = DEFAULT_BIN_DIR, result_dir = DEFAULT_RESULT_DIR, reference_file = DEFAULT_REFERENCE_FILE,
                 result_file = DEFAULT_RESULT_FILE, discarded_trim_file = DEFAULT_DISCARDED_TRIM_FILE,
                 discarded_variant_file = DEFAULT_DISCARDED_VARIANT_FILE, summary_file_location = DEFAULT_SUMMARY_FILE):
    #go over the bin directory, parse bin name from directory names and get the file contents of each
    #create a dictionary of the variant count of each bin and the aggregate them to a matrix of varaint to count per bin
    #receives raw data directory, result directory, result file name, discarded file name
    # can think of creating third dimension for seperating time point and repeat maybe better to do straigh in R
    # or do that by writing a list in one of the each of the csv cells

    #dictionary for end result of table with bin to variant frequencies
    bin_freq_dict = {}

    #if the result directory doesn't exist than create it
    if not os.path.exists(result_dir):
        os.makedirs(result_dir)
    #opne the discarded files and csv writer that will be available for all bins to write to
    discarded_trim_csv  = csv.writer(open(discarded_trim_file + '.csv','wb'))
    discarded_trim_fasta = open(discarded_trim_file + '.fq','wb')
    discarded_variant_csv = csv.writer(open(discarded_variant_file + '.csv','wb'))
    discarded_variant_fasta = open(discarded_variant_file + '.fq','wb')


    reference_dictionary = load_reference_dict(reference_file)

    #counting which point we are at in the walk
    walk_count = -1
    for root, dirs, files in os.walk(bin_dir, followlinks=True):
        #if the directory we're on is one of the working output directories than skip it
        if os.path.split(root)[1] == 'out':
            continue

        walk_count = walk_count + 1
        #if the walk count is zero then we are at the top level and should be parsing the bin names
        if walk_count == 0:
            bins = []
            for dir in dirs:
                #take the directory name and parse out the bin name, put it in the bin name list which we will use later
                bin_name = '_'.join(dir.split('_')[1:])
                bins.append(bin_name)
                # bin_freq_dict[bin_name] = variant_frequency(bin_name, )
        # if there are files in the file list this means we are inside a bin and should be running the variant count
        if len(files) > 0:
            # we will use the bin name as the key for this bin in the bin-variant_freq dictionary
            bin_freq_dict[bins[walk_count-1]] = bin_variant_frequency(bins[walk_count-1],root,files,reference_dictionary,
                                                                      discarded_trim_csv,discarded_trim_fasta,
                                                                      discarded_variant_csv,discarded_variant_fasta,
                                                                      summary_file_location)
    # a variable to sum all reads mapped
    all_sum = 0
    #for each dictionary of each bin sum all values and put into the sum vairable and prin it
    for freq_dict in bin_freq_dict.values():
        all_sum = all_sum + sum(freq_dict.values())

    # create the result file csv and open it as a csv writer
    result_csv = csv.writer(open(result_file, 'wb'))
    result_csv.writerow(['variant']+ bins)
    # for each variant (which should exist in all bins) go over the each freq dictionary find the variants frequency
    # and add it to a list which will then be writen as a row in the csv
    for variant in  bin_freq_dict[bins[0]].keys():
        variant_freq_list = [variant]
        for freq_dict in bin_freq_dict.values():
            all_sum = all_sum + sum(freq_dict.values())
            variant_freq_list.append(freq_dict[variant])
    #     #create a matrix with variants as rows and bins as columns
        result_csv.writerow( variant_freq_list)

def align_and_trim_primers_nosw(seq, f_primer = DEFAULT_FORWARD_PRIMER, r_primer = DEFAULT_REVERSE_PRIMER,
                           sw = DEFAULT_SW_ALIGN , umi_length = DEFAULT_UMI_LENGTH,
                           max_r_shift_length = DEFAULT_R_SHIFT_LENGTH, design_max_length = DEFAULT_DESIGN_MAX_LENGTH,
                           design_min_length = DEFAULT_DESIGN_MIN_LENGTH ):
    if len(seq) < umi_length + len(f_primer) + design_min_length + len(r_primer):
        return ('entire sequence shorther than minimu possible read length','', '',str(seq)),1
    umi = seq[:umi_length]
    r_shift_overhang = ''
    fragment = ''
    read_f_primer_area = seq[umi_length:len(f_primer) + umi_length]
    print f_primer
    print read_f_primer_area
    f_primer_percent_identity = get_percent_identity(f_primer,read_f_primer_area )
    if f_primer_percent_identity < 0.05:
        fragment_start = len(f_primer) + umi_length
        for r_shift_length in xrange(0,max_r_shift_length+1):
            if r_shift_length==0:
                read_r_primer_area = seq[-len(r_primer):]
            else:
                read_r_primer_area = seq[-len(r_primer)-r_shift_length:-r_shift_length]
            r_primer_percent_identity = get_percent_identity(r_primer,read_r_primer_area )
            if r_primer_percent_identity  < 0.05:
                fragment_end = -len(r_primer)-r_shift_length
                fragment = seq[fragment_start:fragment_end]
                if len(fragment) <= design_max_length and len(fragment) >= design_min_length:
                    return (str(umi),str(fragment)),0
                else:
                    # print len(fragment)
                    # print f_primer_percent_identity
                    # print r_primer_percent_identity
                    # print 'fragment not correct length',umi, r_shift_overhang,str(fragment)
                    return ('fragment not correct length',umi, r_shift_overhang,str(fragment)),4
        # print r_primer_percent_identity
        # print 'reverse primer alignment not full',umi, r_shift_overhang,str(seq)
        return ('reverse primer alignment not full',umi, r_shift_overhang,str(seq)),3
    else:
        # print f_primer_percent_identity
        # print 'forward primer alignment not full',umi, r_shift_overhang,str(seq)
        return ('forward primer alignment not full',umi, r_shift_overhang,str(seq)),2




def get_percent_identity(primer, read_primer_area):
    assert len(primer) == len(read_primer_area)
    dif = 0
    for index,base in enumerate(primer):
        if base != read_primer_area[index]:
            dif = dif + 1
    return float(dif)/float(len(primer))


def align_and_trim_primers(seq, f_primer = DEFAULT_FORWARD_PRIMER, r_primer = DEFAULT_REVERSE_PRIMER,
                           sw = DEFAULT_SW_ALIGN , umi_length = DEFAULT_UMI_LENGTH,
                           r_shift_length = DEFAULT_R_SHIFT_LENGTH, design_max_length = DEFAULT_DESIGN_MAX_LENGTH,
                           design_min_length = DEFAULT_DESIGN_MIN_LENGTH ):
    # a method to align a sequence to primers, trim the primers and validate if all the parts are the right length
    # variables: seq: sequence to align, f_primer: forward primer, r_primer: reverse primer,
    # scoring_matrix: Smith Waterman scoring matrix, umi_length: primer design UMI length,
    # r_shift_length: reverse primer design shift length, design_max_length: maximum length of designs,
    # design_min_length: minimum length of designs
    # returns UMI and sequence if the sequence passes all tests, returns fails string, umi, reverse shift overhang and
    # fragment if one of the tests isn't passed
    umi = ''
    r_shift_overhang = ''
    fragment = ''
    f_alignment = sw.align(f_primer,seq)
    if f_alignment.r_end - f_alignment.r_pos== len(f_primer):
        umi = seq[:f_alignment.q_pos]
        if len(umi) == umi_length:
            fragment_start =f_alignment.q_end
            r_alignment = sw.align(r_primer,seq)
            if r_alignment.r_end - r_alignment.r_pos == len(r_primer):
                r_shift_overhang =seq[r_alignment.q_end:]
                if len(r_shift_overhang)<= r_shift_length:
                    fragment_end =r_alignment.q_pos
                    fragment = seq[fragment_start:fragment_end]
                    if len(fragment) <= design_max_length and len(fragment) >= design_min_length:
                        return (str(umi),str(fragment)),0
                    else:
                        # print len(fragment)
                        # print f_alignment.dump()
                        # print r_alignment.dump()
                        # print 'fragment not correct length',umi, r_shift_overhang,str(fragment)
                        return ('fragment not correct length',umi, r_shift_overhang,str(fragment)),5
                else:
                    # print r_alignment.dump()
                    # print 'reverse shift overhang not correct length',umi, r_shift_overhang,str(seq)
                    return ('reverse shift overhang not correct length',umi, r_shift_overhang,str(seq)),4
            else:
                # print f_alignment.dump()
                # print 'reverse primer alignment not full',umi, r_shift_overhang,str(seq)
                return ('reverse primer alignment not full',umi, r_shift_overhang,str(seq)),3
        else:
            # print f_alignment.dump()
            # print 'UMI not correct length',umi, r_shift_overhang,str(seq)
            return ('UMI not correct length',umi, r_shift_overhang,str(seq)),2
    else:
        # print f_alignment.dump()
        # print 'forward primer alignment not full',umi, r_shift_overhang,str(seq)
        return ('forward primer alignment not full',umi, r_shift_overhang,str(seq)),1


def get_umi_counts(bin, merged_seq_file_location,discarded_trim_csv,discarded_trim_fasta,trimmed_reads,fragment_umi_dict_location):
    #todo condisder having an output to json or xml function for this because umi count is undetermined
    # if os.stat(summary_file_location).st_size == 0:
    #     all =0
    #     trim_num = 0
    # else:
    #     summary_file = open(summary_file_location,'rb')
    #     all = int(summary_file.next())
    #     trim_num = int(summary_file.next())
    #     summary_file.close()
    fragment_umi_dict = shelve.open(fragment_umi_dict_location,writeback=True)
    trim_num = 0
    shelve_counter = 0
    for record in SeqIO.parse(open(merged_seq_file_location,'rb'),'fastq'):
        # all = all + 1
        print '+++++++++++stragith attempt+++++++++++'
        print record
        trimmed,fail_stage = align_and_trim_primers_nosw(record.seq)
        if len(trimmed)>2:
            print '------------------reverse attempt------------------'
            rev_trimmed,rev_fail_stage = align_and_trim_primers_nosw(record.seq.reverse_complement())
            if rev_trimmed == 2:
                trimmed = rev_trimmed
            else:
                if fail_stage > rev_fail_stage:
                    fail = trimmed
                else:
                    fail = rev_trimmed

        if len(trimmed)==2:
            print '==================== SUCCSESS ===================='
            trim_num = trim_num + 1
            umi,fragment = trimmed
            if fragment in fragment_umi_dict:
                # fragment_umi_dict[fragment] = fragment_umi_dict[fragment]
                if umi in fragment_umi_dict[fragment]['umi_dict']:
                    fragment_umi_dict[fragment]['umi_dict'][umi].append(record)
                else:
                    fragment_umi_dict[fragment]['umi_dict'][umi] = [record]
                fragment_umi_dict[fragment]['count'] = fragment_umi_dict[fragment]['count'] + 1
                # fragment_umi_dict[fragment] = temp_fragment
            else:
                fragment_umi_dict.setdefault(fragment,{'umi_dict':{umi:[record]},'count':1})
            shelve_counter = shelve_counter + 1
            if shelve_counter%100 ==0:
                fragment_umi_dict.sync()
        else:
            print '~~~~~~~~~~~~~~~~~~~~~~ failure ~~~~~~~~~~~~~~~~~~~~~~'

            SeqIO.write(record,discarded_trim_fasta,'fastq')
            discarded_trim_csv.writerow([bin] + [fail])
    # summary_file = open(summary_file_location,'wb')
    # summary_file.write(str(str(all)) +'\n')
    # summary_file.write(str(str(trim_num)) +'\n')
    # summary_file.close()
    trimmed_reads = trimmed_reads + trim_num

    for fragment, dictionary in fragment_umi_dict.items():
        umi_set = set(dictionary['umi_dict'].keys())
    return fragment_umi_dict, trimmed_reads


def main(argv):
    # opts, args = getopt.getopt(argv,"bin:resd:resf:disf:fp:rp:ref:hd")
    # for opt, arg in opts:
    #     if opt == 'bin':
    #         bin_dir = arg
    #     if opt == 'resd':
    #         result_dir = arg
    #     if opt == 'resf':
    #         result_file = arg
    #     if opt == 'disf':
    #         discarded_file = arg
    #     if opt == 'fp':
    #         forward_primer = arg
    #     if opt == 'rp':
    #         reverse_primer = arg
    #     if opt == 'ref':
    #         reference_file = arg
    #     if opt == 'hd':
    #         home_dir= arg
    # discarded_trim_csv = csv.writer(open('C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\fitseq_sample_data\\discarded_trim.csv','wb'))
    # discarded_trim_fasta = open('C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\fitseq_sample_data\\discarded_trim.fq','wb')
    # discarded_variant_csv = csv.writer(open('C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\fitseq_sample_data\\discarded_variant.csv','wb'))
    # discarded_variant_fasta = open('C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\fitseq_sample_data\\discarded_variant.fq','wb')
    # reference_location = 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\reference_variant_full_sequences.tab'
    # summary_file_location = 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\fitseq_sample_data\\fitseq_sample_summary.txt'
    # summary_file = open( summary_file_location ,'wb')
    # summary_file.close()
    # bin_name = '1_anc'
    # merged_seq_file_location='C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\fitseq_sample_data\\1_anc.file_1.M.fq'
    # fragment_umi_dict_location = 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\fitseq_sample_data\\frag_umi_dict.shelve'
    # before = time.time()
    # fragment_umi_counts = get_umi_counts(bin_name,merged_seq_file_location,discarded_trim_csv, discarded_trim_fasta,0,fragment_umi_dict_location)
    # print time.time()-before
    # print fragment_umi_counts
    # reference_dictionary = load_reference_dict(reference_location = reference_location)
    # variant_frequencies = count_variants(bin_name,fragment_umi_counts,reference_dictionary,
    #                                                 discarded_variant_csv,discarded_variant_fasta)
    # all_freqs = 0
    # found_variants = 0
    # for variant,freq in variant_frequencies.items():
    #     if int(freq) > 0:
    #         found_variants = found_variants +1
    #     all_freqs = all_freqs + freq
    # summary_file = open(summary_file_location,'ab')
    # summary_file.write(str(found_variants)+'\n')
    # summary_file.write(str(all_freqs)+'\n')
    # home_dir = argv[0]
    # ref_file = home_dir + argv[1]
    # bin_dir = home_dir + argv[2]
    # res_dir = home_dir  + argv[3]
    # res_file = home_dir  + argv[4]
    # discarded_trim_file = home_dir  + argv[5]
    # discarded_variant_file= home_dir  + argv[6]
    # summary_file =home_dir  + argv[7]
    # seq ='TTCTTGTTGCAGCTCTTCGCCTTTACGCATATGCGCATCGTTGGTCAGTTCCATCGGTTCAAACATCTAGTATTTCTCCTCTTTAATCGGTGGCTAGCATTATACCTAGGACTGAGCTAGCTGTCAGGGCGCGCCATGACTAAGCTTTTCATTGTC'
    # print align_and_trim_primers(seq)
    # print align_and_trim_primers_nosw(seq)
    # # go_over_bins(bin_dir,res_dir,ref_file,res_file,discarded_trim_file,discarded_variant_file,summary_file,)
    #
    bin_name =  argv[0]
    bin_location = argv[1]
    home_dir = argv[2]
    ref_file = argv[3]
    res_file = argv[4]
    discarded_trim_file = argv[5]
    discarded_variant_file= argv[6]
    summary_file = argv[7]
    design_frequency = bin_variant_frequency(bin_name,bin_location,ref_file,discarded_trim_file,discarded_variant_file,summary_file)
    result_csv = csv.writer(open(res_file,'wb'))
    result_csv.writerow(['',bin_name])
    for design,frequency in design_frequency.items():
        result_csv.writerow([design,frequency])
    print 'run completed'


if __name__ == "__main__":
    #call main giving arguments as so 1-home directory 2-reference file name 3-bin directory 4-result directory
    # 5- result file 6-discarded trim file 7-discarded variant file 7-summary_file

     main(sys.argv[1:])


