#!/usr/bin/python
import sys
import traceback #Dv
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
__author__="mariamrizkallah"
__date__ ="$Jan 18, 2012 4:29:17 PM$"
try:
    print "Welcome to abi2fasta.py"
    count = 0
    path = sys.argv[1]
    output_dir = sys.argv[2]
    output_file = sys.argv[3]
    handle_fasta = open(output_dir + output_file+".fasta", "a")
    handle_fastq = open(output_dir + output_file+".fastq", "a")
    handle_qual = open(output_dir + output_file+".qual", "a")
    #Ref: http://bogdan.org.ua/2007/08/12/python-iterate-and-read-all-files-in-a-directory-folder.html
    ls = os.listdir(path)
    for file in ls:
        ext = os.path.splitext(file)[1]
        rt = os.path.splitext(file)[0]
        if ext == ".ab1":
            print path+file
            for seq_record in SeqIO.parse(path+file, "abi"):
                seq_record.id = path.strip("/")+'_'+seq_record.id+'_'+seq_record.name
                #print seq_record.id, seq_record.name, seq_record.seq
                #record = SeqRecord(seq=seq_record.seq, id=path.strip("/")+'_'+seq_record.id+'_'+seq_record.name, description="")
                handle_fasta.write(seq_record.format("fasta"))
                handle_fastq.write(seq_record.format("fastq"))
                handle_qual.write(seq_record.format("qual"))
            #Generator returns a dict {'MySeq': SeqRecord(seq=Seq('AATTGGCC', IUPACUnambiguousDNA()), id='MySeq', name='abi_name', description=", dbxrefs=[])}
            #record_dict = SeqIO.to_dict(SeqIO.parse(path+file, "abi"))
            #seq_record = str(record_dict.values()).strip('[]')
            #for index, record in enumerate(SeqIO.parse(path+file, "abi")):
                #print "ID = %s, length %i, with %i features" % (record.id, len(record.seq), len(record.features))
            #append to one file
            #SeqIO.convert(path+file, "abi", handle_fasta, "fasta")
            #SeqIO.convert(path+file, "abi", handle_fastq, "fastq")
            #SeqIO.convert(path+file, "abi", handle_qual, "qual")
            #separate files
            #SeqIO.convert(path+file, "abi", output_dir+rt+".fastq", "fastq")
            #SeqIO.convert(path+file, "abi", output_dir+rt+".fasta", "fasta")
            count += 1
    print "Converted %i records" % count
except traceback, t:
    print "[Usage]: ./abi2fasta abi_dir output_dir output_file"
    print print_tb(t)