# -*- coding: utf-8 -*-
"""
Created on Wed Dec  3 16:20:15 2014

@author: ernestmordret
"""

from SuffixTree import SubstringDict
from Bio import SeqIO
from os import chdir

chdir('/home/labs/pilpel/avivro/workspace/data/') # don't forget to change the path variable
SD=SubstringDict()
record_list=[""]

"""
    This is where you load your fasta database and create the suffix tree.
    If you want to look at the reverse complement, you should add these sequences
    explicitely to the db. record_list will store the designs in the right order,
    which will help you to fish the matched sequences using their index.
"""

for record in SeqIO.parse(open('/home/labs/pilpel/avivro/workspace/data/reference_variant_full_sequences.fa','rU'),'fasta'):
    record_list.append(record)
    SD._addToTree(unicode(record.seq))

S1=SD._trees[0]
r1=S1.root()
S2=SD._trees[1]
r2=S2.root()

S=SD._trees[0]
# """Here, you should check the SD.trees variable. If there are more than 1 element there, call me, otherwise proceed"""
r=S.root()

def getAllIndices(node):
    """
    Helper function, not interesting.
    Given a node, return all string ids of the node and its
    children.  Written nonrecursively to avoid potential stack
    problems.  If performance suffers, we can code this in C."""
    stack = [node]
    results = []

    while stack:
        node = stack.pop()
        results.extend(getLeafIds(node))
        child, num_children = node.children(), node.num_children()

        for i in range(num_children):
            stack.append(child)
            child = child.next()

    return results


def getLeafIds(node):
    """
    Helper function, not interesting.
    """
    results = []

    for i in range(node.num_leaves()):
        results.append(node.leaf(i+1)[2])

    return results


def match(seq_quality,threshold=55,node=r,current_score=0,matched_seq='',edge_string=''):
    """ seq_quality : seq zipped with quality

    Match is for matching without allowing for the detection of new mutants.
    You specify a threshold of quality (threshold), and the function will
    return db sequences that match the input read, such that the sum of the
    qualities of the mismatched bases is <= threshold.
    """

    if not seq_quality:
        return [(matched_seq,current_score,node)]

    inp_base,score = seq_quality[0]
    descendants=[]

    if score+current_score <= threshold or inp_base=='N':
        possible_mismatches=set('ATCG').symmetric_difference(inp_base)
        if not edge_string:
            for i in possible_mismatches:
                new_node=node.find_child(i)
                if new_node is not None:
                    descendants.extend(match(seq_quality[1:],threshold,new_node,score+current_score,matched_seq+i,new_node.edgestr()))
        else:
            for i in possible_mismatches:
                if i==edge_string[0]:
                    descendants.extend(match(seq_quality[1:],threshold,node,score+current_score,matched_seq+i,edge_string[1:]))

    if not edge_string:
        new_node=node.find_child(inp_base)
        if new_node is not None:
            descendants.extend(match(seq_quality[:],threshold,new_node,current_score,matched_seq,new_node.edgestr()))
    else:
        if inp_base == edge_string[0]:
            descendants.extend(match(seq_quality[1:],threshold,node,current_score,matched_seq+inp_base,edge_string[1:]))
    return descendants

def match_error_number_with_tree_num(seq,node_number,max_error=2,nodes=(r1,r2)):
    node = nodes[node_number-1]
    return match_error_number(seq,max_error,node)

def match_error_number(seq,max_error=2,node=r,error_number=0,matched_seq='',edge_string=''):
    """ seq_quality : seq zipped with quality

    I modified the script so that it accepts a discrete max_error parameter.
    I haven't tested it yet, so a bug is possible.
    for an input sequence, it returns : the matched sequence, the number of errors associated with the match,
    and the "node". The node can be used to then retrieve the db sequence the read was matched to.
    """

    if not seq:
        return [(matched_seq,error_number,node)]

    inp_base = seq[0]
    descendants=[]

    if error_number < max_error or inp_base=='N':
        possible_mismatches=set('ATCG').symmetric_difference(inp_base)
        if not edge_string:
            for i in possible_mismatches:
                new_node=node.find_child(i)
                if new_node is not None:
                    descendants.extend(match_error_number(seq[1:],max_error,new_node,error_number+1,matched_seq+i,new_node.edgestr()))
        else:
            for i in possible_mismatches:
                if i==edge_string[0]:
                    descendants.extend(match_error_number(seq[1:],max_error,node,error_number+1,matched_seq+i,edge_string[1:]))

    if not edge_string:
        new_node=node.find_child(inp_base)
        if new_node is not None:
            descendants.extend(match_error_number(seq[:],max_error,new_node,error_number,matched_seq,new_node.edgestr()))
    else:
        if inp_base == edge_string[0]:
            descendants.extend(match_error_number(seq[1:],max_error,node,error_number,matched_seq+inp_base,edge_string[1:]))
    return descendants


def match_with_error(seq_quality,threshold=55,node=r,current_score=0,matched_seq='',edge_string='',mut=False):
    """ seq_quality : seq zipped with quality

    Same as match, but this time also allowing for the detection of mutated designs
    """

    if not seq_quality:
        return [(matched_seq,current_score,node)]

    inp_base,score = seq_quality[0]
    descendants=[]

    if score+current_score <= threshold or inp_base=='N' or mut==False:
        possible_mismatches=set('ATCG').symmetric_difference(inp_base)
        if not edge_string:
            for i in possible_mismatches:
                new_node=node.find_child(i)
                if new_node is not None:
                    if not mut:
                        descendants.extend(match_with_error(seq_quality[1:],threshold,new_node,current_score,matched_seq+i,new_node.edgestr(),True))
                    descendants.extend(match_with_error(seq_quality[1:],threshold,new_node,score+current_score,matched_seq+i,new_node.edgestr(),mut))
        else:
            for i in possible_mismatches:
                if i==edge_string[0]:
                    if not mut:
                        descendants.extend(match_with_error(seq_quality[1:],threshold,node,current_score,matched_seq+i,edge_string[1:],True))
                    descendants.extend(match_with_error(seq_quality[1:],threshold,node,score+current_score,matched_seq+i,edge_string[1:],mut))

    if not edge_string:
        new_node=node.find_child(inp_base)
        if new_node is not None:
            descendants.extend(match_with_error(seq_quality[:],threshold,new_node,current_score,matched_seq,new_node.edgestr(),mut))
    else:
        if inp_base == edge_string[0]:
            descendants.extend(match_with_error(seq_quality[1:],threshold,node,current_score,matched_seq+inp_base,edge_string[1:],mut))
    return descendants



def my_filter(records):
    for rec in records:
        if 20<len(rec)<33: #here, apply your own length based filterning if needed
            seq_quality=zip(rec.seq,rec.letter_annotations["phred_quality"]) #use only rec.seq as input if you don't care about the quality
            m=match(seq_quality,threshold=10) # pick the threshold wisely if you use the quality, otherwise use m=match_error_number(rec.seq,max_error=2)
            if rRNAs.isdisjoint(set(sum([getAllIndices(i[2]) for i in m],[]))):
                yield rec


# fastq_parser = SeqIO.parse('/Users/ernestmordret/Desktop/st_alignment/xah', "fastq") #insert path to the fastq file
# if you simply want to use the match_error_number functions in a loop, just call : m=match_error_number(rec.seq,max_error=2)
# m will be a list of tuples where each tuple contains (matched_seq,error_number,node).
# to get the db matches related to your read, you should call "getAllIndices(m[2])"


# SeqIO.write(my_filter(fastq_parser), open('test.fastq','w',0), "fastq")

