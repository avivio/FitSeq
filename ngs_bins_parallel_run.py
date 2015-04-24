__author__ = 'dell7'

import subprocess
import time
import csv
import os
import sys

#the location of the sample directory and the home directory for runs on the sample test set
# DEFAULT_HOME_DIR =  "/home/labs/pilpel/avivro/workspace/data/fitseq_sample_data/multiple_bins_example/"
# DEFAULT_SAMPLE_DIR = DEFAULT_HOME_DIR + 'bins/'

#the location of the sample directory and the home directory for runs on the fitseq raw data
DEFAULT_HOME_DIR =  "/home/labs/pilpel/avivro/workspace/data/fitseq_raw_data/150330_D00257_0179_AC6FFDANXX/"
DEFAULT_SAMPLE_DIR = DEFAULT_HOME_DIR + 'Unaligned_fastq/'

#the location of the reference table
DEFAULT_REFERENCE_FILE = "/home/labs/pilpel/avivro/workspace/data/reference_variant_full_sequences.tab"


def go_over_samples(record_umi_sets,sample_dir  = DEFAULT_BIN_DIR, home_dir = DEFAULT_HOME_DIR, ref_file = DEFAULT_REFERENCE_FILE):
    #runs over the files in the bin directory and sends each to be processed as a separete job on wexac
    #receives raw data directory, the home directory of the entire pipeline and the location of the reference file
    #returns nothing but creates two files in the final result directory containing lists of the result and summary file locations

    #date time string for tracking files
    date_time =  (time.strftime("%d-%m-%y-%H%M"))

    #create the directory to store all the results fot this run
    all_results_dir = home_dir + 'results_' + date_time + '/'
    if not os.path.exists(all_results_dir):
          os.makedirs(all_results_dir)

    #directory storing all the output files of each run, so they can be counted at the end and see how many returned
    counter_directory = all_results_dir + 'counter' + date_time + '/'
    if not os.path.exists(counter_directory):
          os.makedirs(counter_directory)

    #lists of the files that need to be procesessed after the run
    result_files_list = []
    summary_files_list = []
    counter_file_list = []


    #counting which point we are at in the walk over all the directories in the sample directory
    walk_count = -1
    #walking over the sample directory
    for root, dirs, files in os.walk(sample_dir, followlinks=True):
        #if the directory we're on is one of the working output directories than skip it
        if os.path.split(root)[1] == 'out':
            continue

        walk_count = walk_count + 1
        #if the walk count is zero then we are at the top level and should be parsing the sample names
        if walk_count == 0:
            samples = []
            for dir in dirs:
                #take the directory name and parse out the sample name, put it in the sample name list which we will use later
                sample_name = '_'.join(dir.split('_')[1:])
                samples.append(sample_name)


        # if there are files in the file list this means we are inside a sample and should be running the design count
        if len(files) > 0:
            #the sample name is the name before this walk went into the direcotry of the sample
            sample_name = samples[walk_count-1]
            #create the result directory fot this sample
            res_dir = all_results_dir + sample_name + '_results_' + date_time + '/'
            if not os.path.exists(res_dir):
                os.makedirs(res_dir)
            #create the result file where we'll record the design frequencies for this sample
            res_file = res_dir + sample_name + '_frequency_' + date_time +'.csv'

            # add the file location to the result file list
            result_files_list.append(res_file)

            #create discarded read file names for the reads discarded at trim stage and at match to design stage
            discarded_trim_file = res_dir + sample_name + '_discarded_trimmed_' + date_time
            discarded_match_file = res_dir + sample_name + '_discarded_match_' + date_time

            #create the summary file where we write the run stats, then add the file to the list of locations
            summary_file = res_dir + sample_name + '_summary_' + date_time +'.txt'
            summary_files_list.append(summary_file)

            #the counter file is in fact the log file where the run output is directed too
            counter_file = counter_directory + sample_name + '_counter_' + date_time +'.txt'
            counter_file_list.append(counter_file)

            #the umi output directory is where the fastas of each designs umis are written to if the user requires it
            umi_output_file = res_dir + 'umi_output/'

            #create the command for the wexac run by joining all the params together
            # for more memory -R "rusage[mem=4000]" change 4000 to how much memory you want
            # -N is to send the email at the end of the run -o it to write to a new output file which is in fact the counter logfile
            # the params for the python script are the sample anme, the root as the sample directory, home directory,
            # the reference file, the result file location, the discarded file location, the summary file location,
            #the umi output direcotry, and a boolean if needed to record the umi output
            command = ' '.join(['bsub  -N -o' ,counter_file,"-q new-all.q /apps/RH6U4/blcr/0.8.5/bin/cr_run python ./ngs_pipeline.py",
                sample_name, root, home_dir, ref_file,res_file ,  discarded_trim_file , discarded_match_file ,
                summary_file,umi_output_file,record_umi_sets])

            #print and send the command to the shell
            print command
            subprocess.call(command, shell = True)

    #create the directory where the final results will be collected
    final_result_directory = all_results_dir + 'final_result_' + date_time + '/'
    print final_result_directory
    if not os.path.exists(final_result_directory):
      os.makedirs(final_result_directory)

    #write the list of both the resutl file locations and the summary file locations to the directory by joining the list on the return character
    all_result_file_location = final_result_directory + 'all_results.txt'
    all_result_file = open(all_result_file_location,'wb')
    all_result_file.write("\n".join(result_files_list))
    all_summary_file_location = final_result_directory + 'all_summaries.txt'
    all_summary_file = open(all_summary_file_location,'wb')
    all_summary_file.write("\n".join(summary_files_list))


if __name__ == "__main__":
    #recieves a boolean value if you want to record the strings of the umis of every design to fasta files
    record_umi_sets = sys.argv[1]
    #runs the go over samples method
    go_over_samples(record_umi_sets)
