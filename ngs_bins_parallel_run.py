__author__ = 'dell7'

import subprocess
import time
import csv
import os

# DEFAULT_HOME_DIR =  "/home/labs/pilpel/avivro/workspace/data/fitseq_sample_data/multiple_bins_example/"
# DEFAULT_BIN_DIR = DEFAULT_HOME_DIR + 'bins/'

DEFAULT_HOME_DIR =  "/home/labs/pilpel/avivro/workspace/data/fitseq_raw_data/150330_D00257_0179_AC6FFDANXX/"
DEFAULT_BIN_DIR = DEFAULT_HOME_DIR + 'Unaligned_fastq/'

DEFAULT_REFERENCE_FILE = "/home/labs/pilpel/avivro/workspace/data/reference_variant_full_sequences.tab"


def wait_for_results(counter_directory,counter_file_list):
    # count =  len([name for name in os.listdir(counter_directory) if os.path.isfile(name)])
    count =  len(os.listdir(counter_directory))
    while count < len(counter_file_list):
        # count =  len([name for name in os.listdir(counter_directory) if os.path.isfile(name)])
        count =  len(os.listdir(counter_directory))
    success_counter = 0
    time.sleep(5)
    for filename in os.listdir(counter_directory):
        open_file = open(counter_directory + filename,'rb')
        line = open_file.read()
        if line.find('run completed')>-1:
            success_counter = success_counter + 1
    if success_counter == len(counter_file_list):
        return True





def go_over_bins(bin_dir  = DEFAULT_BIN_DIR, home_dir = DEFAULT_HOME_DIR, ref_file = DEFAULT_REFERENCE_FILE):
    #go over the bin directory, parse bin name from directory names and get the file contents of each
    #create a dictionary of the variant count of each bin and the aggregate them to a matrix of varaint to count per bin
    #receives raw data directory, result directory, result file name, discarded file name
    # can think of creating third dimension for seperating time point and repeat maybe better to do straigh in R
    # or do that by writing a list in one of the each of the csv cells

    #dictionary for end result of table with bin to variant frequencies
    date_time =  (time.strftime("%d-%m-%y-%H%M"))
    all_results_dir = home_dir + 'results_' + date_time + '/'
    if not os.path.exists(all_results_dir):
          os.makedirs(all_results_dir)


    counter_directory = all_results_dir + 'counter' + date_time + '/'
    if not os.path.exists(counter_directory):
          os.makedirs(counter_directory)


    bin_freq_dict = {}
    result_files_list = []
    summary_files_list = []
    counter_file_list = []


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
            bin_name = bins[walk_count-1]
            res_dir = all_results_dir + bin_name + '_results_' + date_time + '/'
            if not os.path.exists(res_dir):
                os.makedirs(res_dir)
            res_file = res_dir + bin_name + '_frequency_' + date_time +'.csv'
            result_files_list.append(res_file)
            discarded_trim_file = res_dir + bin_name + '_discarded_trimmed_' + date_time
            discarded_variant_file = res_dir + bin_name + '_discarded_variant_' + date_time
            summary_file = res_dir + bin_name + '_summary_' + date_time +'.txt'
            summary_files_list.append(summary_file)
            counter_file = counter_directory + bin_name + '_counter_' + date_time +'.txt'
            counter_file_list.append(counter_file)
            command = ' '.join(['bsub -R "rusage[mem=4000]" -N -o' ,counter_file,"-q new-all.q /apps/RH6U4/blcr/0.8.5/bin/cr_run python ./ngs_pipeline.py",
                bin_name, root, home_dir, ref_file,res_file ,  discarded_trim_file , discarded_variant_file ,summary_file])
            print command
            # subprocess.call(command, shell = True)
    # print result_files_list
    # print summary_files_list
    # # # success = wait_for_results(counter_directory,counter_file_list)
    #
    # if success:
    #     collect_all_results(date_time,all_results_dir,result_files_list,summary_files_list)
    # else:
    #     print 'run failed'

def collect_all_results(date_time,all_results_dir,result_files_list,summary_files_list):
    final_result_directory = all_results_dir + 'final_result' + date_time + '/'
    if not os.path.exists(final_result_directory):
      os.makedirs(final_result_directory)
    final_result_csv_location = final_result_directory + 'final_result_' + date_time + '.csv'
    final_result_csv = csv.writer(open(final_result_csv_location,'wb'))
    first_result_csv = csv.reader(open(result_files_list[0],'rb'))
    result_csvs = []
    for result_file in result_files_list[1:]:
        result_csvs.append(csv.reader(open(result_file,'rb')))
    for line in first_result_csv:
        row = []
        row= line
        for result_csv in result_csvs:
            result_line = result_csv.next()
            row.append(result_line[1])
            final_result_csv.writerow(row)
    final_summary_csv_location = final_result_directory + 'final_summary_' + date_time + '.csv'
    final_summary_csv = csv.writer(open(final_summary_csv_location,'wb'))
    first_summary_csv = csv.reader(open(summary_files_list[0],'rb'))
    summary_csvs = []
    for summary_file in summary_files_list[1:]:
        summary_csvs.append(csv.reader(open(summary_file,'rb')))
    for line in first_summary_csv:
        row = []
        row = line
        for summary_csv in summary_csvs:
            sumamry_line = summary_csv.next()
            row.append(sumamry_line[1])
            final_summary_csv.writerow(row)

if __name__ == "__main__":
    go_over_bins()
    # counter_directory = "/home/labs/pilpel/avivro/workspace/data/fitseq_sample_data/multiple_bins_example/results_19-04-15-1632/counter19-04-15-1632/"
    # counter_file_list =    ["/home/labs/pilpel/avivro/workspace/data/fitseq_sample_data/multiple_bins_example/results_19-04-15-1632/counter19-04-15-1632/D_8_counter_19-04-15-1632.txt",
    #  "/home/labs/pilpel/avivro/workspace/data/fitseq_sample_data/multiple_bins_example/results_19-04-15-1632/counter19-04-15-1632/1_anc_counter_19-04-15-1632.txt"]
    # print wait_for_results(counter_directory,counter_file_list)
    # date_time = '19-04-15-1259'
    # all_results_dir = "/home/labs/pilpel/avivro/workspace/data/fitseq_sample_data/multiple_bins_example/results_19-04-15-1259/"
    # result_file_list = ['/home/labs/pilpel/avivro/workspace/data/fitseq_sample_data/multiple_bins_example/results_19-04-15-1259/D_8_results_19-04-15-1259/D_8_frequency_19-04-15-1259.csv',
    #                     '/home/labs/pilpel/avivro/workspace/data/fitseq_sample_data/multiple_bins_example/results_19-04-15-1259/1_anc_results_19-04-15-1259/1_anc_frequency_19-04-15-1259.csv']
    # summary_file_list =['/home/labs/pilpel/avivro/workspace/data/fitseq_sample_data/multiple_bins_example/results_19-04-15-1259/D_8_results_19-04-15-1259/D_8_summary_19-04-15-1259.txt',
    #                     '/home/labs/pilpel/avivro/workspace/data/fitseq_sample_data/multiple_bins_example/results_19-04-15-1259/1_anc_results_19-04-15-1259/1_anc_summary_19-04-15-1259.txt']

    # collect_all_results(date_time,all_results_dir,result_file_list,summary_file_list)