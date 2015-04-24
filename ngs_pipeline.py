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

#goodman primers
# DEFAULT_FORWARD_PRIMER='GGCGCGCCATGACTAAGCTTTTCATTGTCATGC'
# DEFAULT_REVERSE_PRIMER = 'CATATGCGTAAAGGCGAAGAGCTGCTGTGTAGATCT'

#home directory where everthing happens
DEFAULT_HOME_DIR = '/home/labs/pilpel/avivro/workspace/'
#the designs refrences in the format of <id> <sequence> in tab delimited format
DEFAULT_REFERENCE_FILE = DEFAULT_HOME_DIR + 'data/reference_design_full_sequences.tab'

#the location of the directory with the fastq files of each sample
DEFAULT_SAMPLE_DIR = DEFAULT_HOME_DIR  + "data/fitseq_sample_data/Project_fitseq_1anc/"
# DEFAULT_SAMPLE_DIR = DEFAULT_HOME_DIR  + 'FitSeq/data/goodman_raw_data/Project_goodman'

#the directory of the results of the pipeline
# DEFAULT_RESULT_DIR = DEFAULT_HOME_DIR  + 'FitSeq/data/ngs_sample_data/ngs_pipeline_result'
DEFAULT_RESULT_DIR = DEFAULT_HOME_DIR  + 'data/fitseq_sample_data/ngs_pipeline_result_1anc_1804/'

#the file with the result frequency matrix
# DEFAULT_RESULT_FILE = DEFAULT_RESULT_DIR + '/Project_avivro_result.csv'
DEFAULT_RESULT_FILE = DEFAULT_RESULT_DIR + 'fitseq_sample_result.csv'

#the files that were discarded at the trim umi count stage
# DEFAULT_DISCARDED_FILE = DEFAULT_RESULT_DIR + '/Project_avivro_discarded.csv'
DEFAULT_DISCARDED_TRIM_FILE = DEFAULT_RESULT_DIR + 'fitseq_sample_discarded_trim'


#the files that were discarded at the design count stage
# DEFAULT_DISCARDED_FILE = DEFAULT_RESULT_DIR + '/Project_avivro_discarded.csv'
DEFAULT_DISCARDED_DESIGN_FILE = DEFAULT_RESULT_DIR + 'fitseq_sample_discarded_design'

#summary file to summarize stats. line 1:merged reads 2:trimmed reads 3:number of unique fragments found 4:number of unique umis for all fragments
DEFAULT_SUMMARY_FILE= DEFAULT_RESULT_DIR + 'fitseq_sample_summary.txt'

#primer trimming constants
DEFAULT_PRIMER_PERCENT_IDENTITY = 0.9
DEFAULT_UMI_LENGTH = 9
DEFAULT_R_SHIFT_LENGTH = 4
DEFAULT_DESIGN_MAX_LENGTH = 94
DEFAULT_DESIGN_MIN_LENGTH = 91


#constants for smith watermant alignment
DEFAULT_MATCH_SCORE= 2
DEFAULT_MISMATCH_SCORE = -1
DEFAULT_SCORING_MATRIX = swalign.NucleotideScoringMatrix(DEFAULT_MATCH_SCORE, DEFAULT_MISMATCH_SCORE)
DEFAULT_SW_ALIGN = swalign.LocalAlignment(DEFAULT_SCORING_MATRIX)


def load_reference_dict(reference_location = DEFAULT_REFERENCE_FILE):
    #method for loading the refence file into a design counting dictionary
    #recieves the reference file path
    #returns a dictionary with the reference design sequences as keys and a list of the id and an empty set as values

    #opens the refence file as a seqio module
    handle = open(reference_location, "rb")
    reference_dict = {}
    for record in SeqIO.parse(handle, "tab") :
        #sequences are fed into the reference dict where the sequence is the key and the value is a list with the id
        # and the umi set
        reference_dict[str(record.seq.upper().reverse_complement())] = [record.id,set()]
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


def trim_all_merged_sequences(sample, merged_seq_file_location):
    #method to take all the merged sequences and trims them to the variable regions
    #recieves the sample name, and the merged sequences file path

    #open the merged sequences which are outputed as gzips
    merged_seq_file = gzip.open(merged_seq_file_location,'rb')
    #create new trimmed file location
    trimmed_merged_seq_file_location = merged_seq_file_location[:-7] + 'T.fq'
    #open file as csv writer
    trimmed_merged_seq_csv = csv.writer(open(trimmed_merged_seq_file_location, 'wb'))


    for record in SeqIO.parse(merged_seq_file, "fastq"):
        #open as fastq send to trim method and write to csv so we also have the sample and record name for discarded records
        trimmed_seq = trim_to_restricition_site(str(record.seq))
        trimmed_merged_seq_csv.writerow([sample, record.name, trimmed_seq[1]])
    return trimmed_merged_seq_file_location



def run_seqprep(file_set,location,file_number,sample,all_reads,merged_reads,
                adapter_1=DEFAULT_FORWARD_PRIMER, adapter_2=DEFAULT_REVERSE_PRIMER):
        #this method runs the seqprep tool on a set of files from a certain sample and results in a merged file
        # important we no longer use seqprep to trim the adapters since we need to trim the primers on our own to get the umis
        #recieves the a tuple with the locations of the forward and reverse read input files, the file number, the sample number
        #the forward and reverse adapter
        #returns the merged result file name as well as the number of input reads and number of reads merged

        #create output directory path
        output_dir = '/'.join([location,'out/'])
        #if the the output directory doesn't exist than create one


        #names for the forward read output file, the reverse read output file, and the merged read output file
        f_prefix = ''.join([output_dir,'output.',os.path.split(file_set['F'])[1]])
        r_prefix = ''.join([output_dir,'output.',os.path.split(file_set['R'])[1]])
        merged_result_file = output_dir+sample + '.file_' + str(file_number)+'.M.fq.gz'
        summary_file = output_dir+sample + '.file_' + str(file_number)+'.summary.txt'

        #call the seqprep tool using the subprocesses library which integratges with linux os
        #paramters are the ones goodman used, could tweak them again to fit not using the primers as adapters
        #use pipe and Popen so that I can get the stderr feed and parse it to get the run stats
        pipe = subprocess.Popen([
        'SeqPrep',
        '-f', file_set['F'],
        '-r', file_set['R'],
        '-1', f_prefix[:-6]+'.fq.gz',
        '-2', r_prefix[:-6]+'.fq.gz',
        '-3', f_prefix[:-6]+'.disc.fq.gz',
        '-4', r_prefix[:-6]+'.disc.fq.gz',
        '-s', merged_result_file,
        '-X', '1', '-g', '-L', '5']
        , stdout=subprocess.PIPE,stderr=subprocess.PIPE, )

        #get the run stats from teh standard error
        #the standard error is the second value in the pipe list
        output = pipe.communicate()[1]
        #turn the string into a file like object
        output_file = StringIO.StringIO(output)
        #loop over the output file object and find the correct lines, get the value out of the relevant lines
        for line in output_file:
            if line.find('Pairs Processed:')> -1:
                # print 'in pairs processed!!!'
                # print line.replace('Pairs Processed:','').strip()
                all_reads =  all_reads + int(line.replace('Pairs Processed:','').strip())
            if line.find('Pairs Merged:')>-1:
                # print 'in pairs processed!!!'
                # print line.replace('Pairs Merged:','').strip()
                merged_reads = merged_reads + int(line.replace('Pairs Merged:','').strip())

        #unzip the result file in the same location since seqprep only returns as a gz file
        subprocess.call(['gzip', '-d','-f', merged_result_file])

        #change the merged result file name
        merged_result_file = merged_result_file[:-3]
        return merged_result_file, all_reads, merged_reads


def count_designs(sample,fragment_umi_record_dict,reference_dictionary, discarded_design_csv, discarded_design_fasta):
    #a method that receives the trimmed merged sequences of one sequencing result file and counts the design frequency
    #returns a dictionary with the design names and a count of them in this file
    # the method saves  all the sequences that don't map to a discarded sequences csv
    #recieves the trimmed merged file path, the reference sequence dictionary and a csv writer object to write
    # the discarded sequences
    # TODO maybe create a method to summarize the discarded and mapped sequences
    # TODO create an output varaint frequencies to file mode



    for fragment,umi_dict in fragment_umi_record_dict.iteritems():
        #parse sample, sequnce name in sequencing output file and sequence from the csv
        umi_record_dict = umi_dict['umi_dict']
        umi_set = set(umi_record_dict.keys())
        #if the reference dictionary has a design with this sequence, count it
        reverse = str(Seq(fragment).reverse_complement())
        if reference_dictionary.has_key(str(fragment)):
            reference_dictionary[str(fragment)][1] = reference_dictionary[str(fragment)][1] + len(umi_set)
        #if not check reverse complement
        elif reference_dictionary.has_key(reverse):
            reference_dictionary[reverse][1] = reference_dictionary[reverse][1] + len(umi_set)
        #else write the sample, name and sequence to discarded csv
        else:
            for umi,record_list in umi_record_dict.iteritems():
                for record in record_list:
                    discarded_design_csv.writerow([sample,record.name,fragment,umi])
                    SeqIO.write(record,discarded_design_fasta,'fastq')
            # print '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    #take the sequnce to id:frequency dictionary and transform it to a id:frequency dictionary
    design_count = {design:frequency for design, frequency in reference_dictionary.values()}
    return design_count


def count_design_umis(design_umi_dictionary,umi_output_location,record_umi_sets = False):
    #method that counts the umis for every design
    # recieves a dictionary of design sequences and their id and umi sets, the umi output location and if the umi output should be recorded
    # returns a dictioinary with the design id and the umi count

    #create design umi count dictionary
    design_count = {}

    #for each design sequence and id umi set pair in dictionary
    for design,value in design_umi_dictionary.iteritems():
        #find the id and the umi set
        design_id = value[0]
        design_umi_set = value[1]
        #if the user wants to record the umi set to a fasta send to the method
        if record_umi_sets:
            write_umi_set(design_id,design_umi_set,umi_output_location)
        # if len(design_umi_set) > 0:
        #     print design_id
        #     print design_umi_set

        #set the design frequencey to be the length of the set
        design_count[design_id] = len(design_umi_set)
    return design_count

def write_umi_set(design_id,design_umi_set,umi_output_location):
    #method to record the umis of a certain design to a fasta file for further analysis
    #receives the design id, the umi set, and the output location

    #if the output directory doesn't exist create it
    if not os.path.exists(umi_output_location):
        os.makedirs(umi_output_location)
    #create the file for the output and ipen it for writing
    umi_file_name = umi_output_location + design_id + '_umi_output.fa'
    umi_file = open(umi_file_name,'wb')
    #write the fasta as a title and sequnce
    for index, umi in enumerate(design_umi_set):
        title = '>' + 'umi_' + str(index) + '_' + umi + "\n"
        umi_file.write(title)
        umi_file.write(umi  + "\n")




def sample_design_frequency(sample_name,location,ref_file,discarded_trim_file,discarded_design_file,summary_file_location,
                          umi_output_location,record_umi_sets):
    #method that goes over the contents of a directory finds F R file pairs and sends them to be merged, trimmed and
    #counted. this method will return a dictionary of design ids and numbers of counts for a sample
    #receives the name of the sample, the directory of the sequencing files, the reference file, the discaded files for
    # trimming and design matching, the summary files, the umi output location and a boolean indicating if the umis should be recorded
    #returns a dictionary with the frequenceis of each design in this sample
    # TODO output design frequency to file option

    #create a list for the files that will be merged which is the size of all the R F pairs but starting from 1
    files = os.listdir(location)

    #create the intermediate output data dir for this stage
    output_dir = '/'.join([location,'out/'])
    if not os.path.exists(output_dir):
          os.makedirs(output_dir)

    # load the reference dictionary into the dictionary we'll use to record the results
    design_umi_dictionary = load_reference_dict(ref_file)

    #create the result dir for this stage
    result_dir = path.dirname(discarded_trim_file)
    if not os.path.exists(result_dir):
                os.makedirs(result_dir)

    #open all the discarded reads files
    discarded_trim_csv  = csv.writer(open(discarded_trim_file + '.csv','wb'))
    discarded_trim_fasta = open(discarded_trim_file + '.fq','wb')
    discarded_design_csv = csv.writer(open(discarded_design_file + '.csv','wb'))
    discarded_design_fasta = open(discarded_design_file + '.fq','wb')

    #create a list of mergable files, I have no idea how this works
    mergable_files =[{} for _ in xrange(len(files)/2 + 1)]

    #go over all the files
    for file in files:
        #if it's an intermediate data directory skip over it
        if file == 'out':
            continue
        #split the file in order to parse it
        split_file = file.split('_')
        #get the file number
        file_num = split_file[-1].split('.')[0]

        #create the dictionary for the f and r files if it doesn't exist
        #parse the read and use it to set the direction value
        file_read = split_file[-2]
        if file_read == 'R1':
            file_direction = 'F'
        if file_read == 'R3':
            file_direction = 'R'
        #put the file directory concatenated to the name in the current file number slot in the F or R slot
        cur_dict =  mergable_files[int(file_num)]
        cur_dict[file_direction] = '/'.join([location,file])

    #create all the run stat variables which will be recorded in the summary
    all_reads = 0
    merged_reads = 0
    trimmed_reads = 0
    design_count_pre_umi = 0

    #go over the files in the mergable file lists
    for index,file_set in  enumerate(mergable_files):
        # if the file set doesn't have any files than don't loop
        if len(file_set)<1:
            continue
        #take the mergable files list of sublists and run seqprep on it
        merged_result_file,all_reads, merged_reads =  run_seqprep(file_set,location,index,sample_name,all_reads, merged_reads)

        #send the merged files to be trimmed annd matched to the design librarym get the resulting frequencies and the run stats
        design_umi_dictionary,trimmed_reads,design_count_pre_umi =  get_umi_counts(sample_name,merged_result_file,discarded_trim_csv,discarded_trim_fasta,
                                        trimmed_reads,design_umi_dictionary,discarded_design_csv,
                                        discarded_design_fasta,design_count_pre_umi)
    #take the design frequency dictionary and change it so that the key is the id and not the sequence, count UMIs
    design_frequencies = count_design_umis(design_umi_dictionary, umi_output_location,record_umi_sets)

    #go over the dictionary and aggregate the run stats
    found_designs = 0
    design_count_post_umi = 0
    # for design,freq in sample_var_freq_dict.iteritems():
    for design,freq in design_frequencies.iteritems():
        if int(freq) > 0:
            found_designs = found_designs +1
        design_count_post_umi = design_count_post_umi + freq
    #write the summary file for this sample
    summary_csv = csv.writer(open(summary_file_location, 'wb'))
    summary_csv.writerow(['',sample_name])
    summary_csv.writerow(['all_reads',all_reads])
    summary_csv.writerow(['merged_reads',merged_reads])
    summary_csv.writerow(['trimmed_reads',trimmed_reads])
    summary_csv.writerow(['found_designs',found_designs])
    summary_csv.writerow(['design_count_pre_umi',design_count_pre_umi])
    summary_csv.writerow(['design_count_post_umi',design_count_post_umi])
    return design_frequencies


def align_and_trim_primers_nosw(seq, f_primer = DEFAULT_FORWARD_PRIMER, r_primer = DEFAULT_REVERSE_PRIMER,
                            umi_length = DEFAULT_UMI_LENGTH, primer_percent_identity = DEFAULT_PRIMER_PERCENT_IDENTITY,
                           max_r_shift_length = DEFAULT_R_SHIFT_LENGTH, design_max_length = DEFAULT_DESIGN_MAX_LENGTH,
                           design_min_length = DEFAULT_DESIGN_MIN_LENGTH ):
    #does an alignment to primers and trim algorithm with no indels and no sequence shift, just mismatch,
    # receives the sequence and can receive the primers, the length of the umi, the identity you want to the primer,
    # the length of the reverse shift and the limits of the design sizes
    # returns a trimmed sequence, and the umi, and a value of the fail stage (0 if succeeded)

    #if the entire read is shorter the the shortest possible sequence to match than fail
    if len(seq) < umi_length + len(f_primer) + design_min_length + len(r_primer):
        # print ('entire sequence shorter than minimum possible read length','', '',str(seq)),1
        return ('entire sequence shorter than minimum possible read length','', '',str(seq)),1
    #get the umi
    umi = str(seq[:umi_length])
    #create the reverse shift string
    r_shift_overhang = ''
    #create the fragment string
    fragment = ''
    #substring the area in the read that should be the forward primer
    read_f_primer_area = seq[umi_length:len(f_primer) + umi_length]

    # compare the forward primer region on the read to the forward primer
    f_primer_percent_identity = get_percent_identity(f_primer,read_f_primer_area )
    #if the primer area is above the level of mismatches than continue write the point the fragment starts and continue
    if f_primer_percent_identity > primer_percent_identity:
        fragment_start = len(f_primer) + umi_length
        #go over all the possible shifts of the reverse primer
        for r_shift_length in xrange(0,max_r_shift_length+1):
            #special case for substring when shift is 0
            if r_shift_length==0:
                read_r_primer_area = seq[-len(r_primer):]
            else:
                #the reverse primer area is where the primer should be on the sequence condidering the reverse shift
                read_r_primer_area = seq[-len(r_primer)-r_shift_length:-r_shift_length]
            #get the percent identity to the reverse primer
            r_primer_percent_identity = get_percent_identity(r_primer,read_r_primer_area )
            #if percent identity is above threshold
            if r_primer_percent_identity  > primer_percent_identity:
                #the framgent ends in the start of the reverse primer area
                fragment_end = -len(r_primer)-r_shift_length
                #substring the fragment from start to finish
                fragment = seq[fragment_start:fragment_end]
                #if the fragment is inside the limit of lengths
                if len(fragment) <= design_max_length and len(fragment) >= design_min_length:
                    return (str(umi),str(fragment)),0
                else:
                    # if the seq isnt within length limits
                    # print len(fragment)
                    # print f_primer_percent_identity
                    # print r_primer_percent_identity
                    # print 'fragment not correct length',umi, r_shift_overhang,str(fragment)
                    return ('fragment not correct length',umi, r_shift_overhang,str(fragment)),4
        # if the reverse primer isnt on the sequence
        # print r_primer_percent_identity
        # print 'reverse primer alignment not full',umi, r_shift_overhang,str(seq)
        return ('reverse primer alignment not full',umi, r_shift_overhang,str(seq)),3
    else:
        # if the forward primer isnt on the sequence
        # print f_primer_percent_identity
        # print 'forward primer alignment not full',umi, r_shift_overhang,str(seq)
        return ('forward primer alignment not full',umi, r_shift_overhang,str(seq)),2



def get_percent_identity(primer, read_primer_area):
    # a method that get's the percent identity between two sequences
    # accepts 2 sequences of the same length
    #returns percent identity
    assert len(primer) == len(read_primer_area)
    dif = 0
    for index,base in enumerate(primer):
        if base != read_primer_area[index]:
            dif = dif + 1
    return 1-float(dif)/float(len(primer))


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


def get_umi_counts(sample, merged_seq_file_location,discarded_trim_csv,discarded_trim_fasta,trimmed_reads,design_umi_dictionary,
                   discarded_design_csv,discarded_design_fasta,design_count_pre_umi):
    #goes over a merged seq file, trims the sequences to get a fragment and a umi, and searched the design library for a match
    #receives the sample name, the merged file location, objects for all the discarded seq files, the reference file
    # dictionary for storing the matched sequence umis, and the count of the number of designs matched before umi filter
    #returns the reference dictionary with the umi sets, the number of trimmed sequences, and the number of matched sequences before umi filter

    #initialize the number of trimmed sequences
    trim_num = 0
    #parse over the merged file using biopython
    for record in SeqIO.parse(open(merged_seq_file_location,'rb'),'fastq'):
        #attempt to trim the sequence on the straigh sequence, get the trimmed string object and the fail stage
        # print '+++++++++++stragith attempt+++++++++++'
        trimmed,fail_stage = align_and_trim_primers_nosw(record.seq)
        #if the trimmed value is more than two elements long then the trim failed
        if len(trimmed)>2:
            #try to trim again on the reverse complement
            # print '------------------reverse attempt------------------'
            rev_trimmed,rev_fail_stage = align_and_trim_primers_nosw(record.seq.reverse_complement())
            #if the trim succeeded switched the value of trimmed with reverse trimmed
            if len(rev_trimmed) == 2:
                trimmed = rev_trimmed
            #if no trim succeeded take the faill at the further along stage for logging
            else:
                if fail_stage > rev_fail_stage:
                    fail = trimmed
                else:
                    fail = rev_trimmed
        #if one of the trims succeeded try to find the design in the library
        if len(trimmed)==2:
            # print '==================== SUCCSESS ===================='
            #count another trimmed sequence
            trim_num = trim_num + 1
            #get the umi and the fragment strings
            umi,fragment = trimmed
            #if the reference dictionary has a design with this sequence add to the umi set try for both regular and reverse complement
            reverse = str(Seq(fragment).reverse_complement())
            if design_umi_dictionary.has_key(str(fragment)):
                design_umi_dictionary[str(fragment)][1].add(umi)
                #add to matched sequence counter
                design_count_pre_umi = design_count_pre_umi + 1
            #if not check reverse complement
            elif design_umi_dictionary.has_key(reverse):
                design_umi_dictionary[reverse][1].add(umi)
                design_count_pre_umi = design_count_pre_umi + 1
        #if the read was trimmed but no match write to the deiscarded match files, the csv value and the fastq value
            else:
                discarded_design_csv.writerow([sample,record.name,fragment,umi])
                SeqIO.write(record,discarded_design_fasta,'fastq')
        #if the read was not trimmed trimmed write to the deiscarded trim files, the csv value and the fastq valuu
        else:
            # print '~~~~~~~~~~~~~~~~~~~~~~ failure ~~~~~~~~~~~~~~~~~~~~~~'

            SeqIO.write(record,discarded_trim_fasta,'fastq')
            discarded_trim_csv.writerow([sample] + [fail])
    trimmed_reads = trimmed_reads + trim_num

    return design_umi_dictionary,trimmed_reads,design_count_pre_umi

def main(argv):
    # # test get_umi_counts
    # discarded_trim_csv = csv.writer(open('C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\fitseq_sample_data\\discarded_trim.csv','wb'))
    # discarded_trim_fasta = open('C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\fitseq_sample_data\\discarded_trim.fq','wb')
    # discarded_design_csv = csv.writer(open('C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\fitseq_sample_data\\discarded_design.csv','wb'))
    # discarded_design_fasta = open('C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\fitseq_sample_data\\discarded_design.fq','wb')
    # reference_location = 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\reference_design_full_sequences.tab'
    # design_umi_dictionary = load_reference_dict(reference_location)
    # summary_file_location = 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\fitseq_sample_data\\fitseq_sample_summary.txt'
    # summary_file = open( summary_file_location ,'wb')
    # summary_file.close()
    # sample_name = '1_anc'
    # merged_seq_file_location='C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\fitseq_sample_data\\1_anc.file_1.M.fq'
    # # merged_seq_file_location='C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\fitseq_sample_data\\1_anc_discarded_trimmed_21-04-15-1134.fq'
    # fragment_umi_dict_location = 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\fitseq_sample_data\\frag_umi_dict.shelve'
    # before = time.time()
    # fragment_umi_counts,trimmed,designed = get_umi_counts(sample_name, merged_seq_file_location,discarded_trim_csv,discarded_trim_fasta,
    #                                                        0,design_umi_dictionary,discarded_design_csv,
    #                                                        discarded_design_fasta,0)
    # print time.time()-before
    # print trimmed
    # print designed
    # # print fragment_umi_counts
    #
    # #test design frequencies as well
    # design_frequencies = count_design_umis(design_umi_dictionary)


    # test align and trim primers
    # seq ='TTCTTGTTGCAGCTCTTCGCCTTTACGCATATGCGCATCGTTGGTCAGTTCCATCGGTTCAAACATCTAGTATTTCTCCTCTTTAATCGGTGGCTAGCATTATACCTAGGACTGAGCTAGCTGTCAGGGCGCGCCATGACTAAGCTTTTCATTGTC'
    # print align_and_trim_primers(seq)
    # print align_and_trim_primers_nosw(seq)

    # test align_and_trim_primers_nosw
    # fasta_file = 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\fitseq_sample_data\\1_anc_discarded_trimmed_21-04-15-1134.fq'
    # for record in SeqIO.parse(open(fasta_file,'rb'),'fastq'):
    #     align_and_trim_primers_nosw(record.seq)

    #real run
    sample_name =  argv[0]
    sample_location = argv[1]
    home_dir = argv[2]
    ref_file = argv[3]
    res_file = argv[4]
    discarded_trim_file = argv[5]
    discarded_design_file= argv[6]
    summary_file = argv[7]
    umi_output_location = argv[8]
    record_umi_sets = argv[9]

    design_frequency = sample_design_frequency(sample_name,sample_location,ref_file,discarded_trim_file,discarded_design_file,
                                             summary_file,umi_output_location,record_umi_sets)
    result_csv = csv.writer(open(res_file,'wb'))
    result_csv.writerow(['',sample_name])
    for design,frequency in design_frequency.iteritems():
        result_csv.writerow([design,frequency])
    print 'run completed'


if __name__ == "__main__":
    #call main giving arguments as so 1-home directory 2-reference file name 3-sample directory 4-result directory
    # 5- result file 6-discarded trim file 7-discarded design file 7-summary_file

     main(sys.argv[1:])


