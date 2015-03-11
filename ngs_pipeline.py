__author__ = 'avivro'
import os
import subprocess
import re
from Bio import SeqIO
import gzip
from Bio.Seq import Seq
import csv

# DEFAULT_FORWARD_PRIMER='NNNNNNNNCAGCTCTTCGCCTTTACGCATATG'
# DEFAULT_REVERSE_PRIMER ='ATGAAAAGCTTAGTCATGGCG'
DEFAULT_FORWARD_PRIMER='GGCGCGCCATGACTAAGCTTTTCATTGTCATGC'
DEFAULT_REVERSE_PRIMER = 'CATATGCGTAAAGGCGAAGAGCTGCTGTGTAGATCT'

DEFAULT_HOME_DIR = '/home/labs/pilpel/avivro/'
DEFAULT_REFERENCE_FILE = DEFAULT_HOME_DIR + 'FitSeq/data/reference_variant_full_sequences.tab'

# DEFAULT_BIN_DIR = DEFAULT_HOME_DIR  + 'FitSeq/data/ngs_sample_data/Project_avivro'
DEFAULT_BIN_DIR = DEFAULT_HOME_DIR  + 'FitSeq/data/goodman_raw_data/Project_goodman'

# DEFAULT_RESULT_DIR = DEFAULT_HOME_DIR  + 'FitSeq/data/ngs_sample_data/ngs_pipeline_result'
DEFAULT_RESULT_DIR = DEFAULT_HOME_DIR  + 'FitSeq/data/goodman_raw_data/ngs_pipeline_result/Project_goodman'

# DEFAULT_RESULT_FILE = DEFAULT_RESULT_DIR + '/Project_avivro_result.csv'
DEFAULT_RESULT_FILE = DEFAULT_RESULT_DIR + '/goodman_raw_data_result.csv'


# DEFAULT_DISCARDED_FILE = DEFAULT_RESULT_DIR + '/Project_avivro_discarded.csv'
DEFAULT_DISCARDED_FILE = DEFAULT_RESULT_DIR + '/goodman_raw_data_discarded.csv'



def load_reference_dict(reference_location = DEFAULT_REFERENCE_FILE):
    handle = open(reference_location, "rb")
    reference_dict = {}
    for record in SeqIO.parse(handle, "tab") :
        reference_dict[str(record.seq.upper().reverse_complement())] = [record.id,0]
    handle.close()
    return reference_dict

def trim_to_restricition_site(seq):

    re1 = re.compile('(?<=CATATG)(.+)', flags=re.IGNORECASE)
    re2 = re.compile('(.+?)(?=GGCGCGCC)', flags=re.IGNORECASE)
    re3 = re.compile('(?<=CATATG)(.+?)(?=GGCGCGCC)', flags=re.IGNORECASE)

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
    merged_seq_file = gzip.open(merged_seq_file_location,'rb')
    trimmed_merged_seq_file_location = merged_seq_file_location[:-7] + 'T.fq'
    trimmed_merged_seq_csv = csv.writer(open(trimmed_merged_seq_file_location, 'wb'))
    for record in SeqIO.parse(merged_seq_file, "fastq"):
        trimmed_seq = trim_to_restricition_site(str(record.seq))
        trimmed_merged_seq_csv.writerow([bin, record.name, trimmed_seq[1]])
    return trimmed_merged_seq_file_location



def run_seqprep(file_set,location,file_number,bin,
                adapter_1=DEFAULT_FORWARD_PRIMER, adapter_2=DEFAULT_REVERSE_PRIMER):
        output_dir = '/'.join([location,'out/'])
        if not os.path.exists(output_dir):
          os.makedirs(output_dir)

        f_prefix = ''.join([output_dir,'output.',os.path.split(file_set['F'])[1]])
        r_prefix = ''.join([output_dir,'output.',os.path.split(file_set['R'])[1]])
        merged_result_file = output_dir+bin + '.file_' + str(file_number)+'.M.fq.gz'

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

        #this method will also have an output to file mode,
        # and will save all the sequences that don't grep (maybe create a method to summarize)
        #the method will return a dictionary of variant names and counts of them
        #lalalala


def count_variants(result_file_location,reference_dictionary, discarded_csv):
    #this method will also have an output to file mode,
        # and will save all the sequences that don't grep (maybe create a method to summarize)
        #the method will return a dictionary of variant names and counts of them
    # command = ' | '.join([
    #     'grep -ixFf {reference_library}',
    #     'sort',
    #     'uniq -c',
    #     'sort -nr > {output_file}']).format(
    #         reference_library=trimmed_ref_seq_file,
    #         output_file=sequence_output)
    with open(result_file_location, 'rb') as results_file:
        results_csv = csv.reader(results_file)
        for result in results_csv:
            bin = result[0]
            name = result[1]
            result = Seq(result[2].strip())
            if reference_dictionary.has_key(str(result)):
                reference_dictionary[str(result)][1] = reference_dictionary[str(result)][1] + 1
            elif reference_dictionary.has_key(str(result.reverse_complement())):
                reference_dictionary[str(result.reverse_complement())][1] = reference_dictionary[str(result.reverse_complement())][1] + 1

            else:
                # print 'doesn\'t exist in reference'
                discarded_csv.writerow([bin,name,result])
                # print '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

        variant_count = {variant:frequency for variant, frequency in reference_dictionary.values()}
    return variant_count


def bin_variant_frequency(bin_name,location,files,reference_dictionary,discarded_csv,output_to_file = False):
    #this method will return a dictionary of variant ids and numbers of counts, will also have an output to file method
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

    variant_frequencies = [{} for _ in xrange(len(mergable_files)+1)]

    for index,file_set in  enumerate(mergable_files):
        if len(file_set)<1:
            continue
        merged_result_file =  run_seqprep(file_set,location,index,bin_name)
        trimmed_merged_result_file = trim_all_merged_sequences(bin_name,merged_result_file)

        variant_frequencies[index] = count_variants(trimmed_merged_result_file,reference_dictionary,discarded_csv)
    bin_var_freq_dict ={}
    for var_freq_dict in variant_frequencies:
        if len(var_freq_dict.items()) > 0:
            for variant,frequency in var_freq_dict.items():
                if bin_var_freq_dict.has_key(variant):
                    bin_var_freq_dict[variant] = bin_var_freq_dict[variant] + frequency
                else:
                    bin_var_freq_dict[variant]  = frequency
    print bin_name
    print sum(bin_var_freq_dict.values())
    return bin_var_freq_dict
        #then go over the dictionary list, and in each dictionary go over the keys and values
        #create another dictionary, for every key (variant id) see if it exists in the dictionary
        #if it does the new value will be the old one plus the value from the current dictionary
        #if it doesn't exist than create a new entry with the key of the current variant and the currernt value
        #return the dictionary with the aggregated results



def go_over_bins(bin_dir  = DEFAULT_BIN_DIR, result_dir = DEFAULT_RESULT_DIR, result_file = DEFAULT_RESULT_FILE,discarded_file = DEFAULT_DISCARDED_FILE):
    #dictionary for end result of table with bin to variant frequencies
    bin_freq_dict = {}
    if not os.path.exists(result_dir):
        os.makedirs(result_dir)
    discarded_csv = csv.writer(open(discarded_file,'wb'))
    reference_dictionary = load_reference_dict()
    #counting which point we are at in the walk
    walk_count = -1
    for root, dirs, files in os.walk(bin_dir, followlinks=True):
        if os.path.split(root)[1] == 'out':
            continue

        walk_count = walk_count + 1
        #if there is anything in the dirs list it means we are at the top level and should be parsing the bin names
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

    all_sum = 0
    for freq_dict in bin_freq_dict.values():
        all_sum = all_sum + sum(freq_dict.values())
    print 'all'
    print all_sum
    if not os.path.exists(result_dir):
        os.makedirs(result_dir)
    result_csv = csv.writer(open(result_file, 'wb'))
    result_csv.writerow(['variant']+ bins)
    for variant in  bin_freq_dict[bins[1]].keys():
        variant_freq_list = [variant]
        for freq_dict in bin_freq_dict.values():
            all_sum = all_sum + sum(freq_dict.values())
            variant_freq_list.append(freq_dict[variant])
    #     #create a matrix with variants as rows and bins as columns
        result_csv.writerow( variant_freq_list)

    #     #(can think of creating third dimension for seperating time point and repeat)
    #     #go over the binXvariant id X frequency dictionary and input into the matrix
    #     #analyze



if __name__ == "__main__":
    go_over_bins()