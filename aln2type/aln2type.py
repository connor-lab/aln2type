#!/usr/bin/env python3

import argparse
import csv
import os
import re
import sys

from .readfq import readfq
from .scheme import read_scheme_yaml
from .vcf_annotator import Annotator
from operator import itemgetter
from itertools import groupby, product
from collections import OrderedDict
import json
import gzip

## define ambiguous bases, if the base in question is ambigious, then return the set of possible bases or else return the base itself (non-abiguous)
def iupac_to_base(base):
    
    lookup_iupac = { 'R' : ['A', 'G'],
                     'Y' : ['C', 'T'],
                     'S' : ['G', 'C'],
                     'W' : ['A', 'T'],
                     'K' : ['G', 'T'],
                     'M' : ['A', 'C'],
                     'B' : ['C', 'G', 'T'],
                     'D' : ['A', 'G', 'T'],
                     'H' : ['A', 'C', 'T'],
                     'V' : ['A', 'C', 'G'] }

    if lookup_iupac.get(base):
        return lookup_iupac.get(base)

    else:
        return base

##
def deconvolute_IUPAC(var):

    var_iupac = []
    # for each base in the variant bases for the sample
    for base in var['variant-base']:
        # append the IUPAC bases to the var_iupac list
        var_iupac.append(iupac_to_base(base))
    
    var_product = []
    # for cases where there are multiple bases involved in a variant (i.e. indels), transform list of bases into single string of bases
    for element in product(*var_iupac):
        var_product.append(''.join(element))
    # add as new entry to variants
    var['iupac-variant-bases'] = var_product

    return var

## opens the sample's MSA using readfq(), finds the reference from the MSA using the cli arg ref_name, converts sequence to uppercase
def get_ref_seq(msa, refname):

    with open(msa) as f:
        for name, seq, qual in readfq(f):
            if name == refname:
                return seq.upper()
                break

# create list of tuples containing the start position of deletions and their length
def group_dels(dels):
    grouped_dels = []
    # x[0] is the eumeration, x[1] is the deletion. Groups cases where enumeration - deletion position are the same, i.e. 
    # (1 - 6512) == (2 - 6513) == (3 - 6514) 
    for k,g in groupby(enumerate(dels),lambda x:x[0]-x[1]):
        # group is item 1 of each grouping, i.e. the base position
        group = (map(itemgetter(1),g))
        # get the integers of the base positions for each group, turn each group to a list
        group = list(map(int,group))
        # append each group[0] - 1 and the length of the group + 1 to the grouped dels list
        # i.e. the position of the last base before the deletion, plus the size of the deletion
        grouped_dels.append((group[0] - 1, len(group) + 1))
    
    return grouped_dels

def group_ins(ins):
    # res is a list of tuples where
    # tuple[0] is each element in ins (each position where there's an insertion) - 1
    # and tuple[1] is the number of times that element appears, +1 for all elements in ins
    res =  [(el - 1, ins.count(el) + 1) for el in ins]

    # return the set i.e. remove duplicates
    return set(res)

# as above for deletions, but for SNPs compared to the reference genome
def group_snps(snps):  
    grouped_snps = []
    
    for k,g in groupby(enumerate(snps),lambda x:x[0]-x[1]):
        group = (map(itemgetter(1),g))
        group = list(map(int,group))
        grouped_snps.append((group[0], len(group)))

    return grouped_snps

# fix overlapping variant calls
def fix_overlapping_vars(all_vars, rseq, qseq):
    # positions is the positions of all variants in the sample
    positions = [i['position'] for i in all_vars]
    # get the set of all duplicated positions
    duplicate_positions = set([x for x in positions if positions.count(x) > 1])

    indel_positions = []
    # for all variants in the list of variants sorted by type
    for variant in sorted(all_vars, key=lambda k: k['type']):
        # if the variant position is in the set of duplicate positions,
        if variant['position'] in duplicate_positions:
            # if the variant type is a snp, remove it
            if variant['type'] == 'snp':
                all_vars.remove(variant)
                continue
            # widen the variant so that it spans +2 and starts at -2 bases - why?
            variant['var-length'] += 2
            variant['position'] -= 2
            variant['ins-corrected-position'] -= 2
            # get the new start and end positions of this widened variant
            var_position = variant['ins-corrected-position']
            var_end = variant['ins-corrected-position'] + variant['var-length']
            
            # if the variant is a deletion
            if variant['type'] == "del":
                # the variant base is now the base of the sample's sequence at the corrected position
                variant['variant-base'] = rseq[var_position]
                # the variant bases are now the bases of the reference at the corrected positions
                variant['reference-base'] = rseq[var_position:var_end]
                # add the range that the new "deletion" spans to the indel_positions list
                indel_positions.extend(range(var_position, var_end))
            
            # else if the variant type is an insertion
            elif variant['type'] == "ins":
                # as above for the deletion, change ref base and bases inserted
                variant['reference-base'] = rseq[var_position]
                variant['variant-base'] = qseq[var_position:var_end]
                # add the range of this to the indel_positions list
                indel_positions.extend(range(var_position, var_end))
        # if the variant position isn't in the list of duplicated calls and it's a snp
        else:
            if variant['type'] == 'snp':
                # if the snp position is in the newly created list of indels, remove the call
                if variant['position'] in indel_positions:
                    all_vars.remove(variant)

    return all_vars

#takes query seq and reference seq
def update_variants(qseq, rseq):

    # create lists for all
    dels = []
    ins = []
    snps = []
    nocalls = []

    # iterate over the bases (qbase) with enumerate to generate qpos. For each base,
    for qpos, qbase in enumerate(qseq):
        # ins_aware_pos is the base position minus the length of ins (initially 0), i.e. the insertion-aware position
        ins_aware_pos = qpos - len(ins)
        # if the base is not equal to the same base in the reference
        if qbase != rseq[qpos]:
            # if the base is a deletion
            if qbase == "-":
                # deletion
                # add the position (minus current length of ins) to the dels list
                dels.append(ins_aware_pos)
            # else, if the reference base at the same position is a deletion
            elif rseq[qpos] == "-":
                # insertion
                # add this position (minus the current length of ins) to the list of insertion positions
                ins.append(ins_aware_pos)
            # else if it does not match IUPAC code, add this position (minus the current length of ins list) to the list of no calls
            elif re.match('[^ATGCRYSWKMBDHV]', qbase):
                nocalls.append(ins_aware_pos)
            # if none of these, it's a SNP. Add the position (minus the current length of ins list) to the SNPs list
            else:
                # snp
                snps.append(ins_aware_pos)
                
        # if both the sample and ref are - at this position, 
        elif qbase == "-" and rseq[qpos] == "-":
            # insertion but not in this sequence
            # append position (minus current length of ins) to list of insertion positions
            ins.append(ins_aware_pos)

    all_vars = []
    # if there are deletions
    if dels:
        #group_dels creates a list of tuples of the start position and the lengths of deletions
        for start,length in group_dels(dels):
            #correct for insertion if the starting position of a deletion is after an insertion
            #ins_correction is the length of the number of times i in ins is less than the starting position of a deletion
            # i.e. the number of inserted bases before a deletion
            ins_correction = len([i for i in ins if i < start])
            # add ins_correction to starting position
            c_position = start + ins_correction
            # end position with insertion correction
            c_end = start + ins_correction + length
            # reference sequence from corrected start to corrected end positions, i.e. the bases of the deletion
            var_ref_seq = rseq[c_position:c_end]
            # alt seq is the starting base of the deletion
            var_alt_seq = rseq[c_position]
            # make set of data for this deletion
            var = { "type" : "del", 
                    "reference-base" : var_ref_seq, 
                    "position" : start, 
                    "variant-base" : var_alt_seq, 
                    "ins-corrected-position" : c_position, 
                    "var-length" : length }
            # append to list of sets of info on all variants for this sample
            all_vars.append(var)
    
    # if there are no calls, i.e. ambiguous bases, again adjust for insertions and then add to cumulative set of variants for this sample
    if nocalls:
        for start in nocalls:
            ins_correction = len([i for i in ins if i < start])

            c_position = start + ins_correction
            
            var_ref_seq = rseq[c_position]
            var_alt_seq = qseq[c_position]

            var = { "type" : "no-call", 
                    "reference-base" : var_ref_seq, 
                    "position" : start, 
                    "variant-base" : var_alt_seq, 
                    "ins-corrected-position" : c_position, 
                    "var-length" : 1 }

            all_vars.append(var)
    # if there are snps,
    if snps:
        # group sequential snps
        grouped_snps = group_snps(snps)
        # as with deletions above, adjust position for insertions
        for start,length in grouped_snps:
            ins_correction = len([i for i in ins if i < start])

            c_position = start + ins_correction
            c_end = start + ins_correction + length
            # get bases from start to end positions of the (groups of) SNP(s)
            var_ref_seq = rseq[c_position:c_end]
            var_alt_seq = qseq[c_position:c_end]
            # if the SNP length is more than 1, mark as multi nucleotide poly
            if length > 1:
                var = {'type' : 'mnp'}
            # else it's a snp
            else:
                var = {'type' : 'snp'}

            var.update({ "reference-base" : var_ref_seq,
                    "position" : start,
                    "variant-base" : var_alt_seq,
                    "ins-corrected-position" : c_position,
                    "var-length" : length })
            # append to cumulative set of all variant positions
            all_vars.append(var)
            
            # if there is an mnp,split it into individual snp calls but note each is an mnp-snp
            if var['type'] == 'mnp':
            
                for idx,mnp_variant_base in enumerate(var['variant-base'][1:], 1):
                    # get reference base at the variant position
                    mnp_snp_ref_seq = var['reference-base'][idx]
                    # get alt sequence in sample
                    mnp_snp_alt_seq = mnp_variant_base
                    # get mnp start position + how far along the single snp it is (i.e. true position)
                    mnp_position = var['position'] + idx
                    # correct for insertions + position within mnp
                    mnp_c_position = var['ins-corrected-position'] + idx
                    # create entry for each snp in the mnp
                    mnpvar = { "type" : "mnp-snp", 
                            "reference-base" : mnp_snp_ref_seq, 
                            "position" : mnp_position, 
                            "variant-base" : mnp_snp_alt_seq, 
                            "ins-corrected-position" : mnp_c_position, 
                            "var-length" : 1 }
                    # append to set
                    all_vars.append(mnpvar)

    if ins:
        # if there are insertions, get the start position and length from the insertions grouping
        for start,length in group_ins(ins):
            # as with deletions above, ins_correction is the length of the number of times i in ins is less than the starting position of a deletion
            # i.e. the number of inserted bases that occur before the current insertion
            ins_correction = len([i for i in ins if i < start])
            # insertion position correction
            c_position = start + ins_correction
            # corrected end position
            c_end = start + ins_correction + length
            # reference base at the corrected position
            var_ref_seq = rseq[c_position]
            # alternative positions in sample, i.e. bases of the insertion
            var_alt_seq = qseq[start:start+length]
            # if the reference bases at the corrected positions are the same as the insertion sequence, 
            if rseq[c_position:c_end] == var_alt_seq:
                # insertion but not in this sequence
                continue
            # make entry
            var = { "type" : "ins",
                    "reference-base" : var_ref_seq,
                    "position" : start,
                    "variant-base" : var_alt_seq,
                    "ins-corrected-position" : c_position,
                    "var-length" : length }
            # append to set 
            all_vars.append(var)
    # if there is at least one variant in this sample
    if len(all_vars) >= 1:
        # fix overlapping variant calls, see function for notes
        fixed_vars = fix_overlapping_vars(all_vars, rseq, qseq)
        # for each variant in the newly corrected variant calls
        for var in fixed_vars:
            # remove the ins-corrected-position out of the entry for this variant
            ob_ins_c_pos = var.pop('ins-corrected-position')
            # replace with value + 1
            var['ob-ins-corrected-position'] = ob_ins_c_pos + 1 
            # remove the position entry
            ob_pos = var.pop('position')
            # replace with the entry +1
            var['one-based-reference-position'] = ob_pos + 1

            deconvolute_IUPAC(var)
        # return variants sorted by one-based ref position
        return sorted(fixed_vars, key=lambda x: x['one-based-reference-position'])

    else:
        return None


def annotate_translation(gb, vars):
    
    annotator = Annotator(gb_file=gb, vcf_records=vars)
    annotator.annotate_records()

    annotated_vars = annotator.get_vcf_records()

    return annotated_vars

# types variants, taking the sample name, the sequence, the ymls and whether to include no_call_deletions
def type_variants(name, f_variants, variant_types, no_call_deletion=False):

    # below is structured as, for the sample:
    # { position: {variant dict}}
    formatted_vars_s = dict(zip([ i['one-based-reference-position'] for i in f_variants ], f_variants))

    sample_type = { 'sample_id' : name, 'typing' : [] }

    # create nested lists of the positions of deletions, starting at one-based ref position +1, e.g. [[6513, 6514, 6515], [11283, 11284, 11285, 11286, 11287, 11288, 11289, 11290, 11291]]
    deleted_positions = [ 
        list(
            range(
                i['one-based-reference-position'] + 1, 
                i['one-based-reference-position'] + i['var-length']
            )
        ) 
        for i in f_variants if i['type'] == 'del'
        ]
    
    # for each of the variant types (i.e. the ymls), here called vtype
    for vidx,vtype in enumerate(variant_types):
        # summary is everything that isn't the variants and calling-definition sections XXXXX CW this is where the required line would be filtered out
        summary = { i:vtype[i] for i in vtype if i not in [ 'variants', 'calling-definition' ]}
        # append to the typing section of above, the summary info and an empty dict of the sample typing result
        sample_type['typing'].append( {'sample-typing-summary' : summary, 'sample-typing-result' : {} } )

        # empty dictionary
        variant_lists = {}
        # populate with keys and populate variants key with the variants section of the yml
        variant_lists['master'] = { 'variants' : vtype['variants'], 'calling-definition' : {} }
        # for each item in calling-definitions, i.e. each confirmed, probable, low-qc for each yml
        # calling id is confirmed, probable etc and calling definition is the mutations required etc
        for calling_id, calling_definition in vtype['calling-definition'].items():
            
            if 'variants' in calling_definition:
                # has own variants that this call specifically requires
                # example is carnival-shimmy
                # in these cases, exclude variant lines and pull out the different requirements as dicts e.g. {'mutations-required': 5, 'indels-required': 0, 'allowed-wildtype': 4}
                calling_def_dict = { i:calling_definition[i] for i in calling_definition if i != 'variants' }
                # populate each yml's variant_lists with the variants info, and the calling definitions as nested dicts 
                # calling id is confirmed, probable etc 
                variant_lists[calling_id] = {'variants' : calling_definition['variants'], 
                                             'calling-definition' : { calling_id : calling_def_dict }}

            else:
                variant_lists['master']['calling-definition'].update( { calling_id : calling_definition } )

        ## see notes for dict structure at this point, repeated for each yml    
        for name,props in variant_lists.items():
            # for each "name" i.e. master, alt-probable and so on, and the defs for each,
            #make a set of call keys
            calls_keys = [
                'mutation_ref_calls', 
                'indel_ref_calls', 
                'mutation_calls', 
                'mutation_mixed_calls', 
                'indel_calls', 
                'no_calls',
                'no_call_deletion'
                ]
            # populate with 0s
            calls = dict.fromkeys(calls_keys, 0)
            # for each set of variants for each name (aka master, alt-probable)
            for idx,var in enumerate(props['variants']):
                # get the lines in 
                sample_var = formatted_vars_s.get(var['one-based-reference-position'])
                # if there are variant bases in the sample
                if sample_var:
                    
                    # MNP is found in sample and subordinate SNP is part of variant definition - lowercase = sample, UPPER = yml definition
                    if sample_var['type'] == 'mnp' and var['type'] == 'SNP':
                        # sample_variant_base is the set of bases
                        sample_variant_base = [ i[0] for i in sample_var['variant-base'] ]


                    # MNP in definition contains a reference base and sample variant base is subordinate to it
                    elif ( sample_var['type'] == 'mnp' or sample_var['type'] == 'snp' ) and var['type'] == 'MNP':

                        # Definition MNP is longer than the one from sample - grab the reference bases from the definition 
                        if len(var['variant-base']) > len(sample_var['variant-base']):

                            extra_mnp_length = len(sample_var['variant-base']) - len(var['variant-base'])

                            extra_mnp_bases = var['variant-base'][-extra_mnp_length]

                            sample_variant_base = [ i + extra_mnp_bases for i in sample_var['iupac-variant-bases'] ]

                        else: 
                            sample_variant_base = sample_var['iupac-variant-bases']


                    # All other non-special cases
                    else:
                        sample_variant_base = sample_var['iupac-variant-bases']
                    
                    variant_lists[name]['variants'][idx]['sample-call'] = ','.join(sample_variant_base)
               
                    if sample_var['type'] == 'no-call':
                        variant_lists[name]['variants'][idx]['status'] = 'no-call'
                        calls['no_calls'] += 1

                    elif var['variant-base'] in sample_variant_base:

                        if len(sample_variant_base) > 1:
                            
                            variant_lists[name]['variants'][idx]['status'] = 'detect-mixed'

                            calls['mutation_mixed_calls'] += 1

                        elif len(sample_variant_base) == 1:
                            
                            variant_lists[name]['variants'][idx]['status'] = 'detect'

                            if var['type'] == 'insertion' or var['type'] == 'deletion':
                                
                                calls['indel_calls'] += 1

                            else:
                                calls['mutation_calls'] += 1
                        
                                           
                else:
                    if var['reference-base'] == var['variant-base']:
                        variant_lists[name]['variants'][idx]['status'] = 'detect'

                        calls['mutation_calls'] += 1

                    elif no_call_deletion and var['one-based-reference-position'] in [i for s in deleted_positions for i in s]:
                        variant_lists[name]['variants'][idx]['status'] = 'no-call-deletion'
                    
                        calls['no_call_deletion'] += 1

                    else:
                        variant_lists[name]['variants'][idx]['status'] = 'no-detect'

                        if var['type'] == 'insertion' or var['type'] == 'deletion':
                            calls['indel_ref_calls'] += 1
                        
                        else:
                            calls['mutation_ref_calls'] += 1
                     
            sample_type['typing'][vidx]['sample-typing-result'][name] = { 'calls' : calls, 
                                                  'calling-definition' : props['calling-definition'],
                                                  'variants' : variant_lists[name]['variants'] }
          
    return sample_type


def score_typing(sample_type, output_unclassified=False):

    sample_summary = []
    for vtype in sample_type['typing']:

        sample_type_summary = { 'sample_id' : sample_type['sample_id'] }

        defs = []

        for name, properties in vtype['sample-typing-result'].items():
            for cdef,vals in properties['calling-definition'].items():
                defs.append({ 'name' : name, 'level': cdef, 'calling-definition' : vals })

        s_defs = sorted(defs, key=lambda x:x['calling-definition']['mutations-required'])
        
        potential_calling_defs = []
        for definition in s_defs:
             if (vtype['sample-typing-result'][definition['name']]['calls']['mutation_ref_calls'] <= definition['calling-definition']['allowed-wildtype'] and
                    vtype['sample-typing-result'][definition['name']]['calls']['mutation_calls'] >= definition['calling-definition']['mutations-required'] and 
                        vtype['sample-typing-result'][definition['name']]['calls']['indel_calls'] >= definition['calling-definition']['indels-required'] ):
                            potential_calling_defs.append(definition['level'])

        if potential_calling_defs:
            vtype['sample-typing-result']['variant-status'] = potential_calling_defs[-1]

            sample_type_summary.update(
                {
                    'phe-label' : vtype['sample-typing-summary']['phe-label'],
                    'unique-id' : vtype['sample-typing-summary']['unique-id'],
                    'prelim_status' : vtype['sample-typing-result']['variant-status'],
                    'mutation-ref-calls': vtype['sample-typing-result'][definition['name']]['calls']['mutation_ref_calls'],
                    'mutation-mixed-calls': vtype['sample-typing-result'][definition['name']]['calls']['mutation_mixed_calls'],
                    'mutation-calls': vtype['sample-typing-result'][definition['name']]['calls']['mutation_calls'],
                    'indel-ref-calls': vtype['sample-typing-result'][definition['name']]['calls']['indel_ref_calls'],
                    'indel-calls': vtype['sample-typing-result'][definition['name']]['calls']['indel_calls'],
                    'no-calls': vtype['sample-typing-result'][definition['name']]['calls']['no_calls'],
                    'no-calls-deletion': vtype['sample-typing-result'][definition['name']]['calls']['no_call_deletion']
                }
            )

            sample_summary.append(sample_type_summary)

        else:
            vtype['sample-typing-result']['variant-status'] = None

    if output_unclassified:
        
        called_types = [t['sample-typing-result']['variant-status'] for t in sample_type['typing']]
      
        if called_types.count(None) == len(called_types):
            sample_summary.append(
                {   'sample_id' : sample_type['sample_id'],
                    'phe-label' : 'unclassified',
                    'unique-id' : 'unclassified',
                    'prelim_status' : 'unclassified',
                    'mutation-ref-calls': None,
                    'mutation-mixed-calls': None,
                    'mutation-calls': None,
                    'indel-ref-calls': None,
                    'indel-calls': None,
                    'no-calls': None,
                    'no-calls-deletion': None
                }
            )

    return sample_summary, sample_type

def normalise_fn(name):
    return re.sub(r'[^A-Za-z0-9.-]', '_', name)

def make_dir(fp):
    if not os.path.exists(fp):
        os.makedirs(fp)

def write_json(name, typing_data, json_outdir, no_gzip_json=False):

    json_outpath = os.path.abspath(json_outdir)

    make_dir(json_outpath) 

    if no_gzip_json:
        jsonfilename = os.path.join(json_outpath, normalise_fn(name) + '.json')
        with open(jsonfilename, 'w', encoding='UTF-8') as zipfile:
            json.dump(typing_data, zipfile, indent = 4)
    
    else:
        jsonfilename = os.path.join(json_outpath, normalise_fn(name) + '.json.gz')
        with gzip.open(jsonfilename, 'wt', encoding='UTF-8') as zipfile:
            json.dump(typing_data, zipfile, indent = 4)

def write_variant_types(typing_summary, output_csv):
    # It is theoretically possible that a set of genomes will not contain any VoCs
    # This checks to see if the typing summary contains information before writing the csv.
    if typing_summary:
        csvfilename = os.path.abspath(output_csv)

        fieldnames = list(typing_summary[0].keys())

        with open(csvfilename, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            for row in typing_summary:
                writer.writerow(row)
    else:
        print ("No Variants Identified for these sequences from these definitions.")

def write_sample_variant_csv(name, variants, sample_csv_outdir, csv_N=False):

    csv_outpath = os.path.abspath(sample_csv_outdir)

    make_dir(csv_outpath)

    csvfilename = os.path.join(csv_outpath, normalise_fn(name) + '.csv')

    all_fieldnames = list(max(variants, key=len).keys())

    rejected_fields = ['ob-ins-corrected-position']

    fieldnames = [ i for i in all_fieldnames if i not in rejected_fields ]

    fieldnames.insert(0, 'sampleID')

    with open(csvfilename, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in variants:
            #if row['type'] == 'mnp-snp':
            #    continue

            if row['type'] == 'no-call' and csv_N == False:
                continue
            
            csv_row = { k: row[k] for k in row if k not in rejected_fields }
            csv_row['iupac-variant-bases'] = '|'.join(csv_row['iupac-variant-bases'])
            csv_row['sampleID'] = name

            writer.writerow(csv_row)

# removes any "-"s or Ns from the end of the sequence
def remove_terminal_gapns(seq):
    return re.sub(r'(N|-)*$', '', seq)

# simply removes any "-"s from the end of the sequence
def remove_terminal_gaps(seq):
    return re.sub(r'-*$', '', seq)


def check_required(variant_types):
    required_types = {'unclassified':[None]}

    for variant in variant_types:
        req_v = variant.get('requires')
        if type(req_v) != list:
            required_types[variant['unique-id']] = [req_v]
        else:
            required_types[variant['unique-id']] = req_v

    def recursive_check(unique_id, required_dict):

        required_ids = required_dict[unique_id]
       
        if required_ids != [None]:
            sub_requires = [y for x in required_ids for y in required_dict[x]]
        else:
            sub_requires = []
       
        if sub_requires != [] and sub_requires != [None]:
            for req in sub_requires:
                if req != None:
                    recurse = recursive_check(req, required_dict)
                    if recurse != [] and recurse != [None]:
                        sub_requires.extend(recurse)
 
        return sub_requires       
        
    for key in required_types.keys():
        rec_sub_req = recursive_check(key, required_types)        
        if rec_sub_req != [None]:
            required_types[key].extend(rec_sub_req)        
    
    return required_types

def format_requires(typing_summary, required_types):
    added_requires = typing_summary.copy()
    for i in range(len(added_requires)):
        if required_types.get(added_requires[i]['unique-id']) == [None]:
            added_requires[i]['met_required'] = '' 
            added_requires[i]['status'] = added_requires[i]['prelim_status']

        else:
            # get the other variant calls for a sample 
            # if not the same call, add to dict with status 
            sample_other_calls = {}
            for j in range(len(added_requires)):
                if added_requires[j]['sample_id'] == added_requires[i]['sample_id']:
                    if added_requires[j]['unique-id'] == added_requires[i]['unique-id']:
                        pass 
                    else:
                        sample_other_calls[added_requires[j]['unique-id']] = added_requires[j]['prelim_status']          
            final_status_bool = []

            for v in required_types[added_requires[i]['unique-id']]:
                if v in sample_other_calls:
                    if "confirmed" in sample_other_calls[v] or "probable" in sample_other_calls[v]:
                        final_status_bool.append(True)

                    else:
                        final_status_bool.append(False)

                else:
                    final_status_bool.append(False)

            if False not in final_status_bool:
                added_requires[i]['met_required'] = 'TRUE'
                added_requires[i]['status'] = added_requires[i]['prelim_status']
            else:
                added_requires[i]['met_required'] = 'FALSE'
                added_requires[i]['status'] = 'low-qc'
      
    return added_requires

def go(args):

    # reference sequence is get_ref_seq taking msa and reference name, reading the MSA and pulling out the ref sequence (in upper case)
    print("Get Refseq")
    refseq = get_ref_seq(args.msa, args.ref_name)
    print(f"Get Refseq: {args.ref_name}")

    # if it can't find the reference in the MSA, state not found and exit
    if not refseq:
        print(args.ref_name + " not found\n") 
        sys.exit(1)

    # make empty list
    print("Get variants")
    variant_types = []
    # for variant yml in the arg typing_yml, read the scheme and append to the list
    for scheme in args.typing_yaml:
        # read_scheme_yaml just reads each of the yml files in safe_load mode.safe_load recognizes only standard YAML tags https://pyyaml.org/wiki/PyYAMLDocumentation
        variant_types.append(read_scheme_yaml(scheme))
    print("Check for required flags")
    variant_requires = check_required(variant_types)

    # make empty list
    typing_summary = []
    print("Start typing")
    # open the sample's MSA
    with open(args.msa) as f:
        # for each sample in the MSA, if the name = the reference name arg supplied, continue, otherwise
        for name, seq, qual in readfq(f):
            if name == args.ref_name:
                continue
            else:
                # if the arg no_trum_terminal_N is used, then 
                if args.no_trim_terminal_N:
                    # run remove_terminal_gaps on the uppercase sequence
                    # removes any "-"s from the end of the sequence
                    # update_variants - calls variants, removes duplicate contrasting variant calls, fixes iupac
                    variants = update_variants(remove_terminal_gaps(seq.upper()), refseq)
                else:
                    # else if the arg is not used, run remove_terminal_gapns on the uppercase sequence
                    # removes any "-"s or Ns from the end of the sequence
                    # as above with update_variants
                    variants = update_variants(remove_terminal_gapns(seq.upper()), refseq)
                # print the sample's name
                #print(name)
                # if there are variants identified
                if variants:
                    # run type_variants with sample name, the set of variants called, the ymls and whether no_call_deletion is true or false
                    typed_variants = type_variants(name, variants, variant_types, args.no_call_deletion)
                    
                    sample_summary, scored_variants = score_typing(typed_variants, args.output_unclassified)

                    typing_summary.extend(sample_summary)

                    write_json(name, scored_variants, args.json_outdir, args.no_gzip_json)

                    if args.gb:
                        annotated_variants = annotate_translation(args.gb, variants)
                        
                        write_sample_variant_csv(name, annotated_variants, args.sample_csv_outdir, args.csv_N)

                    else:
                        write_sample_variant_csv(name, variants, args.sample_csv_outdir, args.csv_N)
    print("writing out")    
    typing_summary = format_requires(typing_summary, variant_requires)
    write_variant_types(typing_summary, args.summary_csv_outfile)
    
                    
                    

