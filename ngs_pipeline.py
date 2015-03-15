__author__ = 'avivro'
import os
import subprocess
import re
from Bio import SeqIO
import gzip
from Bio.Seq import Seq
import csv
import getopt


#primers act as the adapters in the seqprep tool, need to see how this acts with umi's
# DEFAULT_FORWARD_PRIMER='NNNNNNNNCAGCTCTTCGCCTTTACGCATATG'
# DEFAULT_REVERSE_PRIMER ='ATGAAAAGCTTAGTCATGGCG'
DEFAULT_FORWARD_PRIMER='GGCGCGCCATGACTAAGCTTTTCATTGTCATGC'
DEFAULT_REVERSE_PRIMER = 'CATATGCGTAAAGGCGAAGAGCTGCTGTGTAGATCT'

#home directory where everthing happens
DEFAULT_HOME_DIR = '/home/labs/pilpel/avivro/'
#the variants refrences in the format of <id> <sequence> in tab delimited format
DEFAULT_REFERENCE_FILE = DEFAULT_HOME_DIR + 'FitSeq/data/reference_variant_full_sequences.tab'

#the location of the directory with the fastq files of each sample
# DEFAULT_BIN_DIR = DEFAULT_HOME_DIR  + 'FitSeq/data/ngs_sample_data/Project_avivro'
DEFAULT_BIN_DIR = DEFAULT_HOME_DIR  + 'FitSeq/data/goodman_raw_data/Project_goodman'

#the directory of the results of the pipeline
# DEFAULT_RESULT_DIR = DEFAULT_HOME_DIR  + 'FitSeq/data/ngs_sample_data/ngs_pipeline_result'
DEFAULT_RESULT_DIR = DEFAULT_HOME_DIR  + 'FitSeq/data/goodman_raw_data/ngs_pipeline_result/Project_goodman'

#the file with the result frequency matrix
# DEFAULT_RESULT_FILE = DEFAULT_RESULT_DIR + '/Project_avivro_result.csv'
DEFAULT_RESULT_FILE = DEFAULT_RESULT_DIR + '/goodman_raw_data_result.csv'

#the files that were discarded at the variant count stage
# DEFAULT_DISCARDED_FILE = DEFAULT_RESULT_DIR + '/Project_avivro_discarded.csv'
DEFAULT_DISCARDED_FILE = DEFAULT_RESULT_DIR + '/goodman_raw_data_discarded.csv'



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



def run_seqprep(file_set,location,file_number,bin,
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

        #call the seqprep tool using the subprocesses library which integratges with linux os
        subprocess.check_call([
        'SeqPrep',
        '-f', file_set['F'],
        '-r', file_set['R'],
        '-1', f_prefix[:-6]+'.fq.gz',
        '-2', r_prefix[:-6]+'.fq.gz',
        '-3', f_prefix[:-6]+'.disc.fq.gz',
        '-4', r_prefix[:-6]+'.disc.fq.gz',
        '-s', merged_result_file,
        '-A', adapter_1, '-B', adapter_2,
        '-X', '1', '-g', '-L', '5'])

        #to print seqprep run params
        # print ' '.join([
        # 'SeqPrep',
        # '-f', file_set['F'],
        # '-r', file_set['R'],
        # '-1', f_prefix,
        # '-2', r_prefix,
        # '-3', f_prefix[:-6]+'.disc.fq.gz',
        # '-4', r_prefix[:-6]+'.disc.fq.gz',
        # '-s', merged_result_file,
        # '-A', adapter_1, '-B', adapter_2,
        # '-X', '1', '-g', '-L', '5'])
        return merged_result_file


def count_variants(result_file_location,reference_dictionary, discarded_csv):
    #a method that receives the trimmed merged sequences of one sequencing result file and counts the variant frequency
    #returns a dictionary with the variant names and a count of them in this file
    # the method saves  all the sequences that don't map to a discarded sequences csv
    #recieves the trimmed merged file path, the reference sequence dictionary and a csv writer object to write
    # the discarded sequences
    # TODO maybe create a method to summarize the discarded and mapped sequences
    # TODO create an output varaint frequencies to file mode

    #open trimmed sequence file and loop over it
    with open(result_file_location, 'rb') as results_file:
        results_csv = csv.reader(results_file)
        #for every trimmed file
        for result in results_csv:
            #parse bin, sequnce name in sequencing output file and sequence from the csv
            bin = result[0]
            name = result[1]
            result = Seq(result[2].strip())
            #if the reference dictionary has a variant with this sequence, count it
            if reference_dictionary.has_key(str(result)):
                reference_dictionary[str(result)][1] = reference_dictionary[str(result)][1] + 1
            #if not check reverse complement
            elif reference_dictionary.has_key(str(result.reverse_complement())):
                reference_dictionary[str(result.reverse_complement())][1] = reference_dictionary[str(result.reverse_complement())][1] + 1
            #else write the bin, name and sequence to discarded csv
            else:
                # print 'doesn\'t exist in reference'
                discarded_csv.writerow([bin,name,result])
                # print '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
        #take the sequnce to id:frequency dictionary and transform it to a id:frequency dictionary
        variant_count = {variant:frequency for variant, frequency in reference_dictionary.values()}
    return variant_count


def bin_variant_frequency(bin_name,location,files,reference_dictionary,discarded_csv,output_to_file = False):
    #method that goes over the contents of a directory finds F R file pairs and sends them to be merged, trimemed and
    #counted. this method will return a dictionary of variant ids and numbers of counts for a bin
    #receives the anme of the bin, the directory of the sequncing fils. the list of files, the reference dictionary
    #the discarded sequenc csv writer object, and the write to file flag which isn't used
    # TODO output variant frequency to file option

    #create a list for the files that will be merged which is the size of all the R F pairs but starting from 1
    mergable_files =[{} for _ in xrange(len(files)/2 + 1)]

    #go over all the files
    for file in files:
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

    # if the file set doesn't have any files than don't loop
    for index,file_set in  enumerate(mergable_files):
        if len(file_set)<1:
            continue
        #take the mergable files list of sublists and run seqprep on it
        merged_result_file =  run_seqprep(file_set,location,index,bin_name)
        #trim merged file results from seqprep
        trimmed_merged_result_file = trim_all_merged_sequences(bin_name,merged_result_file)
        #take the trimmed frequencines and count the variants in each file, put the resulting dictionary in the index of the file number
        variant_frequencies[index] = count_variants(trimmed_merged_result_file,reference_dictionary,discarded_csv)
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
    print bin_name
    print sum(bin_var_freq_dict.values())
    return bin_var_freq_dict



def go_over_bins(bin_dir  = DEFAULT_BIN_DIR, result_dir = DEFAULT_RESULT_DIR, result_file = DEFAULT_RESULT_FILE,discarded_file = DEFAULT_DISCARDED_FILE):
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
    #opne the discarded csv writer that will be available for all bins to write to
    discarded_csv = csv.writer(open(discarded_file,'wb'))
    reference_dictionary = load_reference_dict()

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
            bin_freq_dict[bins[walk_count-1]] = bin_variant_frequency(bins[walk_count-1],root,files,reference_dictionary,discarded_csv)
    # a variable to sum all reads mapped
    all_sum = 0
    #for each dictionary of each bin sum all values and put into the sum vairable and prin it
    for freq_dict in bin_freq_dict.values():
        all_sum = all_sum + sum(freq_dict.values())
    print 'all'
    print all_sum

    # create the result file csv and open it as a csv writer
    result_csv = csv.writer(open(result_file, 'wb'))
    result_csv.writerow(['variant']+ bins)
    # for each variant (which should exist in all bins) go over the each freq dictionary find the variants frequency
    # and add it to a list which will then be writen as a row in the csv
    for variant in  bin_freq_dict[bins[1]].keys():
        variant_freq_list = [variant]
        for freq_dict in bin_freq_dict.values():
            all_sum = all_sum + sum(freq_dict.values())
            variant_freq_list.append(freq_dict[variant])
    #     #create a matrix with variants as rows and bins as columns
        result_csv.writerow( variant_freq_list)





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
    go_over_bins()



if __name__ == "__main__":
    main(sys.argv[1:])