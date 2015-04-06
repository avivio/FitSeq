__author__ = 'dell7'
import re

class SamRecord:
    """A record of a sam format output from bowtie2"""

    def __init__(self, record):
        self.record = record
        self.split_record = record.split('\t')
        self.ref = self.get_ref()
        if self.is_sam_record():
            self.cigar = self.get_cigar()
            self.mismatch =  self.get_mismatch()
            self.length = self.get_length()
            self.edit_distance = self.get_edit_distance()
            self.name = self.get_name()

    def get_ref(self):
        return self.split_record[2]

    def get_cigar(self):
        return inner_cigar.match(self.split_record[5]).group(2)

    def get_length(self):
        number = ''
        sum = 0
        for i in self.cigar:
            if i.isdigit():
                number = number + i
            else:
                sum = sum + int(number)
                number = ''
        return sum

    def get_mapq(self):
        return self.split_record[4]

    def get_mismatch(self):
        for field in self.split_record:
            if field.startswith('MD:Z:'):
                return field.replace('MD:Z:','')

    def get_aligment_score(self):
        for field in self.split_record:
            if field.startswith('AS:i:'):
                return int(field.replace('AS:i:',''))

    def get_next_aligment_score(self):
        for field in self.split_record:
            if field.startswith('XS:i:'):
                return int(field.replace('XS:i:',''))

    def get_edit_distance(self):
        for field in self.split_record:
            if field.startswith('NM:i:'):
                return int(field.replace('NM:i:',''))
    def is_sam_record(self):
        return self.ref.find('*') < 0

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            if str(self.ref + self.cigar + self.mismatch).__hash__()== (other.ref + other.cigar + other.mismatch).__hash__():
                if self.cigar == other.cigar:
                    if self.mismatch == other.mismatch:
                        return True
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return str(self.ref + self.cigar + self.mismatch).__hash__()
    def get_name(self):
        return '|'.join([self.ref ,self.cigar ,self.mismatch])


inner_cigar = re.compile('([0-9]*S)?(([0-9]+[^0-9S])+)([0-9]*S)?')



sam_file_location = 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\discarded_merged_reads.sa'
print sam_file_location
sam_file = open(sam_file_location,'rb')
all_counter = 0
length_counter = 0
dist_counter = 0
mutation_dict = {}
for line in sam_file:
    if not line.startswith('@'):
        current_sam_record = SamRecord(line)
        if current_sam_record.is_sam_record():
            # split_line = line.split('\t')
            # seq_dict = {}
            # ref = split_line[2]
            # if ref == '*':
            #     #no result from bowtie
            #     continue
            # print 'ref'
            # print ref
            # seq_dict['cigar'] = inner_cigar.match(split_line[5]).group(2)
            # print 'cigar'
            # print seq_dict['cigar']
            # number = ''
            # sum = 0
            # for i in seq_dict['cigar']:
            #     if i.isdigit():
            #         number = number + i
            #     else:
            #         sum = sum + int(number)
            #         number = ''
            # print sum
            # for field in split_line:
            #     if field.startswith('MD:Z:'):
            #         seq_dict['mismatch'] = field.replace('MD:Z:','')
            #     if field.startswith('AS:i:'):
            #         seq_dict['aligment_score'] = field.replace('AS:i:','')
            #     if field.startswith('XS:i:'):
            #         seq_dict['next_aligment_score'] = field.replace('XS:i:','')
            #     if field.startswith('NM:i:'):
            #         seq_dict['edit_distance'] = field.replace('NM:i:','')
            # print 'mismatch'
            # print seq_dict['mismatch']
            #
            #
            #
            #
            # print 'aligment score'
            # print seq_dict['aligment_score']
            # print 'next aligment score'
            # print seq_dict['next_aligment_score']
            # print 'edit distance'
            # print seq_dict['edit_distance']
            #
            # seq_dict['mapq'] =  split_line[4]
            # print 'mapq'
            # print seq_dict['mapq']
            if current_sam_record.length == mutation_dict[current_sam_record.ref]['length']:
                length_counter = length_counter + 1
                if current_sam_record.name in mutation_dict[current_sam_record.ref]['mutations']:
                    mutation_dict[current_sam_record.ref]['mutations'][current_sam_record.name]['counter'] = \
                        mutation_dict[current_sam_record.ref]['mutations'][current_sam_record.name]['counter'] + 1
                    mutation_dict[current_sam_record.ref]['mutations'][current_sam_record.name]['sam_records'].append(current_sam_record)
                    print current_sam_record.name
                    print  mutation_dict[current_sam_record.ref]['mutations'][current_sam_record.name]
                else:
                    mutation_dict[current_sam_record.ref]['mutations'].setdefault\
                        (current_sam_record.name,{'counter':1,'sam_records':[current_sam_record]})
                if current_sam_record.edit_distance<=3:
                    dist_counter = dist_counter + 1
            all_counter =all_counter + 1
    else:
        split_ref = line.split('	')
        if split_ref[1].startswith('SN:'):
            variant = split_ref[1].replace('SN:','')
            length = split_ref[2].replace('LN:','')
            mutation_dict.setdefault(variant, {'length': int(length), 'mutations': {}})

for mutant,dictionary in mutation_dict.items():
    print mutant
    print len(dictionary['mutations'])
    print dictionary
print length_counter
print all_counter
print dist_counter


