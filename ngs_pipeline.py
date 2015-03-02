__author__ = 'avivro'
import os
import subprocess
import re
from Bio import SeqIO



DEFAULT_FORWARD_PRIMER='NNNNNNNNCAGCTCTTCGCCTTTACGCATATG'
DEFAULT_REVERSE_PRIMER ='ATGAAAAGCTTAGTCATGGCG'
DEFAULT_REFERENCE_FILE = 'C:\Users\dell7\Documents\Tzachi\workspace\data\\reference_variant_full_sequences.fasta'

BIN_DIR = 'C:\Users\dell7\Documents\Tzachi\workspace\data\\ngs_sample_data\Project_avivro'


def load_reference_dict(reference_location = DEFAULT_REFERENCE_FILE):
    handle = open("opuntia.aln", "rU")
    reference_dict = {}
    for record in SeqIO.parse(handle, "tab") :
        reference_dict[record.seq] = [record.id,0]
    handle.close()
    return reference_dict

def trim_to_restricition_site(seq):

    re1 = re.compile('(?<=CATATG)(.+)', flags=re.IGNORECASE)
    re2 = re.compile('(.+?)(?=GGCGCGCC)', flags=re.IGNORECASE)
    re3 = re.compile('(?<=CATATG)(.+?)(?=GGCGCGCC)', flags=re.IGNORECASE)

    trim = [None, None, None]

    # Check if both restriction sites are there
    try:
        hit1 = re3.search(seq)
        hit1.groups()
        trim[0] = 1
        trim[1] = str(hit1.group(1))
        trim[2] = 1
        return trim

    except:
        # Check if only CATATG is present
        try:
            hit2 = re1.search(seq)
            hit2.groups
            ()
            trim[0] = 1
            trim[1] = str(hit2.group(1))
            trim[2] = None
            return trim

        except:
            # Check if only GGCGCGCC is present
            try:
                hit3 = re2.search(seq)
                hit3.groups()
                trim[0] = None
                trim[1] = str(hit3.group(1))
                trim[2] = 1
                return trim

            # No restriction sites present
            except:
                trim[0] = None
                trim[1] = str(seq)
                trim[2] = None
                return trim


def trim_all_merged_sequences(merged_seq_file_location):
    merged_seq_file = open(merged_seq_file_location,'rb')
    trimmed_merged_seq_file_location = merged_seq_file_location[:-7] + 'T' + merged_seq_file_location[:-5]
    trimmed_merged_seq_file = open(trimmed_merged_seq_file_location, 'wb')
    for record in SeqIO.parse(merged_seq_file, "fastq"):
        trimmed_seq = trim_to_restricition_site(record.seq)
        trimmed_merged_seq_file.writelines(trimmed_seq)
    return trimmed_merged_seq_file_location



def run_seqprep(file_set,location,file_number,bin,
                adapter_1=DEFAULT_FORWARD_PRIMER, adapter_2=DEFAULT_REVERSE_PRIMER):
        output_dir = '\\'.join([location,'out\\'])


        f_prefix = ''.join([output_dir,'output.',os.path.split(file_set['F'])[1]])
        r_prefix = ''.join([output_dir,'output.',os.path.split(file_set['R'])[1]])
        merged_result_file = output_dir+bin + '.file_' + str(file_number)+'.M.fq.gz'
        #
        # subprocess.check_call([
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


        print ' '.join([
        'SeqPrep',
        '-f', file_set['F'],
        '-r', file_set['R'],
        '-1', f_prefix,
        '-2', r_prefix,
        '-3', f_prefix[:-6]+'.disc.fq.gz',
        '-4', r_prefix[:-6]+'.disc.fq.gz',
        '-s', merged_result_file,
        '-A', adapter_1, '-B', adapter_2,
        '-X', '1', '-g', '-L', '5'])
        return merged_result_file

        #this method will also have an output to file mode,
        # and will save all the sequences that don't grep (maybe create a method to summarize)
        #the method will return a dictionary of variant names and counts of them
        #lalalala


def count_variants(result_file,reference_dictionary):
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

    with open(result_file, 'rb') as results:
        for result in results:
            if reference_dictionary.has_key(result):
                reference_dictionary[reference_dictionary][1] = reference_dictionary[reference_dictionary][1] + 1
            else:
                print 'doesn\'t exist in reference'
                print result
                print '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

        variant_count = {variant:frequency for variant, frequency in reference_dictionary.values()}
        return variant_count


def bin_variant_frequency(bin_name,location,files,output_to_file = False):
    #this method will return a dictionary of variant ids and numbers of counts, will also have an output to file method
    #create a list for the files that will be merged which is the size of all the R F pairs but starting from 1
    mergable_files =[{}]*(len(files)/2 + 1)
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
        mergable_files[int(file_num)][file_direction] = ''.join([location,file])
    print mergable_files
    variant_frequencies = [{}]*(len(mergable_files)+1)

    for index,file_set in  enumerate(mergable_files):
        index = index + 1
        merged_result_file =  run_seqprep(file_set,location,index,bin_name)
        trimmed_merged_result_file = trim_all_merged_sequences(merged_result_file)
        reference_dictionary = load_reference_dict()
        variant_frequencies[index] = count_variants(trimmed_merged_result_file,reference_dictionary)
    bin_var_freq_dict ={}
    for var_freq_dict in variant_frequencies:
        if len(var_freq_dict.items()) > 0:
            for variant,frequency in var_freq_dict.items():
                if bin_var_freq_dict.has_key(variant):
                    bin_var_freq_dict[variant] = bin_var_freq_dict[variant] + frequency
                else:
                    bin_var_freq_dict[variant]  = frequency
    return bin_var_freq_dict
        #then go over the dictionary list, and in each dictionary go over the keys and values
        #create another dictionary, for every key (variant id) see if it exists in the dictionary
        #if it does the new value will be the old one plus the value from the current dictionary
        #if it doesn't exist than create a new entry with the key of the current variant and the currernt value
        #return the dictionary with the aggregated results



def go_over_bins(bin_dir  = BIN_DIR):
    #dictionary for end result of table with bin to variant frequencies
    bin_freq_dict = {}
    #counting which point we are at in the walk
    walk_count = -1
    for root, dirs, files in os.walk(bin_dir, followlinks=True):
        walk_count = walk_count + 1

        #if there is anything in the dirs list it means we are at the top level and should be parsing the bin names
        if len(dirs) > 0:
            bins = []
            for dir in dirs:
                #take the directory name and parse out the bin name, put it in the bin name list which we will use later
                bin_name = '_'.join(dir.split('_')[1:])
                bins.append(bin_name)
                # bin_freq_dict[bin_name] = variant_frequency(bin_name, )
            # print bins
        # if there are files in the file list this means we are inside a bin and should be running the variant count
        if len(files) > 0:
            # we will use the bin name as the key for this bin in the bin-variant_freq dictionary
            bin_freq_dict[bins[walk_count-1]] = bin_variant_frequency(bins[walk_count-1],root,files )
    print bin_freq_dict
        #create a matrix with variants as rows and bins as columns
        #(can think of creating third dimension for seperating time point and repeat)
        #go over the binXvariant id X frequency dictionary and input into the matrix
        #analyze



if __name__ == "__main__":
    go_over_bins()