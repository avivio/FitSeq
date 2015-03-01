__author__ = 'avivro'
import os


BIN_DIR = 'C:\Users\dell7\Documents\Tzachi\workspace\data\\ngs_sample_data\Project_avivro'


def run_seqprep(file_set):
        #this method will also have an output to file mode,
        # and will save all the sequences that don't grep (maybe create a method to summarize)
        #the method will return a dictionary of variant names and counts of them
        #lalalala
    pass

def grep_to_reference(result_file,reference_file):
    #this method will also have an output to file mode,
        # and will save all the sequences that don't grep (maybe create a method to summarize)
        #the method will return a dictionary of variant names and counts of them
    pass

def variant_frequency(bin_name,location,files,output_to_file = False):
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
    variant_frequencies = [{}]*len(mergable_files)
    for file_set, index in  enumerate(mergable_files):
        result_file = run_seqprep(file_set)
        variant_frequencies[index] = grep_to_reference(result_file,reference_file)

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
            bin_freq_dict[bins[walk_count-1]] = variant_frequency(bins[walk_count-1],root,files )
    print bin_freq_dict
        #create a matrix with variants as rows and bins as columns
        #(can think of creating third dimension for seperating time point and repeat)
        #go over the binXvariant id X frequency dictionary and input into the matrix
        #analyze



if __name__ == "__main__":
    go_over_bins()