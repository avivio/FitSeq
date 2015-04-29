__author__ = 'dell7'
import os
import csv
import sys

#names of the files to read in the final result directory
DEFAULT_RESULT_FILE_NAME = '/all_results.txt'
DEFAULT_SUMMARY_FILE_NAME = '/all_summaries.txt'

#names of the files to create and write to in the final result directory
DEFAULT_FINAL_RESULT_FILE_NAME = '/final_result.csv'
DEFAULT_FINAL_SUMMARY_FILE_NAME = '/final_summary.csv'


DEFAULT_MATCH_COUNT_FILE_NAME = '/all_match_counts.txt'
DEFAULT_FINAL_MATCH_COUNT_FILE_NAME = '/final_match_count.csv'


def collect_all_results(final_result_directory,all_result_file_name = DEFAULT_RESULT_FILE_NAME,
                        all_summary_file_name = DEFAULT_SUMMARY_FILE_NAME,
                        final_result_file_name = DEFAULT_FINAL_RESULT_FILE_NAME,
                        final_summary_file_name = DEFAULT_FINAL_SUMMARY_FILE_NAME,
                        all_match_count_file_name = DEFAULT_MATCH_COUNT_FILE_NAME,final_match_count_file_name = DEFAULT_FINAL_MATCH_COUNT_FILE_NAME):
    #method that locates the result and summary files of a fitseq run, collects all the results and writes a unified result and summary file
    #accepts a result directory which should include a list of result files and a list of summary files.
    #does not return anything but creates a matrix of sample to design frequency, and a unified table of run stats

    #open all results file and read all lines into a list for looping
    all_result_file_location = final_result_directory + all_result_file_name
    all_result_file = open(all_result_file_location,'rb')
    result_files_list = all_result_file.readlines()

    #open all summaries file and read all lines into a list for looping
    all_summary_file_location = final_result_directory + all_summary_file_name
    all_summary_file = open(all_summary_file_location,'rb')
    summary_files_list = all_summary_file.readlines()

    #open all match count file and read all lines into a list for looping
    all_match_count_file_location = final_result_directory + all_match_count_file_name
    all_match_count_file = open(all_match_count_file_location,'rb')
    match_count_files_list = all_match_count_file.readlines()


    #open final result matrix as a csv location in the directory specified by the user
    final_result_csv_location = final_result_directory + final_result_file_name
    final_result_csv = csv.writer(open(final_result_csv_location,'wb'))

    #open the first csv of the result csvs as a seperate file so it will be used to write the design id column
    first_result_csv = csv.reader(open(result_files_list[0].strip(),'rb'))

    # create a list of all the result csvs already open by going over the rest of the result file locations and opening
    # them as csv objects
    result_csvs = []
    for result_file in result_files_list[1:]:
        result_csvs.append(csv.reader(open(result_file.strip(),'rb')))

    #go over the lines of the first result csv, have it be the base of the final matrix row
    for line in first_result_csv:
        row= line
        #add the value of the design frequency from each result file to the final matrix current row
        for result_csv in result_csvs:
            result_line = result_csv.next()
            row.append(result_line[1])
        #write the row in the final csv and continue loop
        final_result_csv.writerow(row)

    #do exactly the same thing for the summaries consider creating method for this
    #create final summary file
    final_summary_csv_location = final_result_directory + final_summary_file_name
    final_summary_csv = csv.writer(open(final_summary_csv_location,'wb'))

    #open first csv
    first_summary_csv = csv.reader(open(summary_files_list[0].strip(),'rb'))
    #create csv list
    summary_csvs = []
    for summary_file in summary_files_list[1:]:
        summary_csvs.append(csv.reader(open(summary_file.strip(),'rb')))
    #use first csv row as base of final row and add the values from each csv
    for line in first_summary_csv:
        row = line
        for summary_csv in summary_csvs:
            summary_line = summary_csv.next()
            row.append(summary_line[1])
        #write the row in the final csv and continue loop
        final_summary_csv.writerow(row)


    #do exactly the same thing for the match count files consider creating method for this
    #create final summary file
    final_match_count_csv_location = final_result_directory + final_match_count_file_name
    final_match_count_csv = csv.writer(open(final_match_count_csv_location,'wb'))

    #open first csv
    first_match_count_csv = csv.reader(open(match_count_files_list[0].strip(),'rb'))
    #create csv list
    match_count_csvs = []
    for match_count_file in match_count_files_list[1:]:
        match_count_csvs.append(csv.reader(open(match_count_file.strip(),'rb')))
    #use first csv row as base of final row and add the values from each csv
    for line in first_match_count_csv:
        row = line
        for match_count_csv in match_count_csvs:
            match_count_line = match_count_csv.next()
            row.append(match_count_line[1])
        #write the row in the final csv and continue loop
        final_match_count_csv.writerow(row)

if __name__ == "__main__":
    #receives the final result directory as command line argument and pasees to method
     collect_all_results(sys.argv[1])


