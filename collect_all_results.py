__author__ = 'dell7'
import os
import csv



def collect_all_results(final_result_directory):
    all_result_file_location = final_result_directory + 'all_results.txt'
    all_result_file = open(all_result_file_location,'rb')
    result_files_list = all_result_file.readlines()
    all_summary_file_location = final_result_directory + 'all_summaries.txt'
    all_summary_file = open(all_summary_file_location,'wb')
    summary_files_list = all_summary_file.readlines()


    final_result_csv_location = final_result_directory + 'final_result.csv'
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
    final_summary_csv_location = final_result_directory + 'final_summary.csv'
    final_summary_csv = csv.writer(open(final_summary_csv_location,'wb'))
    first_summary_csv = csv.reader(open(summary_files_list[0],'rb'))
    summary_csvs = []
    for summary_file in summary_files_list[1:]:
        summary_csvs.append(csv.reader(open(summary_file,'rb')))
    for line in first_summary_csv:
        row = []
        row = line
        for summary_csv in summary_csvs:
            summary_line = summary_csv.next()
            row.append(summary_line[1])
            final_summary_csv.writerow(row)

