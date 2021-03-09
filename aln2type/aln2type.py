#!/usr/bin/env python3

import argparse
import csv
import os
import re
import sys
from .readfq import readfq
from .scheme import read_scheme_yaml
from operator import itemgetter
from itertools import groupby, product
from collections import OrderedDict
import json
import gzip


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

def deconvolute_IUPAC(var):

    var_iupac = []
    for base in var['variant-base']:
        var_iupac.append(iupac_to_base(base))
    
    var_product = []
    for element in product(*var_iupac):
        var_product.append(''.join(element))

    var['variant-base'] = var_product

    return var

def get_ref_seq(msa, refname):

    with open(msa) as f:
        for name, seq, qual in readfq(f):
            if name == refname:
                return seq.upper()
                break

def group_dels(dels):
    grouped_dels = []

    for k,g in groupby(enumerate(dels),lambda x:x[0]-x[1]):
        group = (map(itemgetter(1),g))
        group = list(map(int,group))
        grouped_dels.append((group[0] - 1, len(group) + 1))

    return grouped_dels

def group_ins(ins):

    res =  [(el - 1, ins.count(el) + 1) for el in ins]

    return set(res)

def group_snps(snps):  
    grouped_snps = []

    for k,g in groupby(enumerate(snps),lambda x:x[0]-x[1]):
        group = (map(itemgetter(1),g))
        group = list(map(int,group))
        grouped_snps.append((group[0], len(group)))

    return grouped_snps

def fix_overlapping_vars(all_vars, rseq, qseq):
    
    positions = [i['position'] for i in all_vars]

    duplicate_positions = set([x for x in positions if positions.count(x) > 1])

    indel_positions = []

    for variant in sorted(all_vars, key=lambda k: k['type']):
        
        if variant['position'] in duplicate_positions:

            if variant['type'] == 'snp':
                all_vars.remove(variant)
                continue

            variant['var-length'] += 2
            variant['position'] -= 2
            variant['ins-corrected-position'] -= 2
            var_position = variant['ins-corrected-position']
            var_end = variant['ins-corrected-position'] + variant['var-length']
            
            if variant['type'] == "del":
                
                variant['variant-base'] = rseq[var_position]
                variant['reference-base'] = rseq[var_position:var_end]

                indel_positions.extend(range(var_position, var_end))
            
            
            elif variant['type'] == "ins":
             
                variant['reference-base'] = rseq[var_position]
                variant['variant-base'] = qseq[var_position:var_end]

                indel_positions.extend(range(var_position, var_end))

        else:
            if variant['type'] == 'snp':
                if variant['position'] in indel_positions:
                    all_vars.remove(variant)

    return all_vars

def update_variants(qseq, rseq):

    dels = []
    ins = []
    snps = []
    nocalls = []

    for qpos, qbase in enumerate(qseq):
        ins_aware_pos = qpos - len(ins)

        if qbase != rseq[qpos]:
            if qbase == "-":
                # deletion
                dels.append(ins_aware_pos)

            elif rseq[qpos] == "-":
                # insertion
                ins.append(ins_aware_pos)
            
            elif re.match('[^ATGCRYSWKMBDHV]', qbase):
                nocalls.append(ins_aware_pos)

            else:
                # snp
                snps.append(ins_aware_pos)
                

        elif qbase == "-" and rseq[qpos] == "-":
            # insertion but not in this sequence
            ins.append(ins_aware_pos)

    all_vars = []

    if dels:

        for start,length in group_dels(dels):
            ins_correction = len([i for i in ins if i < start])
            
            c_position = start + ins_correction
            c_end = start + ins_correction + length

            var_ref_seq = rseq[c_position:c_end]
            var_alt_seq = rseq[c_position]

            var = { "type" : "del", 
                    "reference-base" : var_ref_seq, 
                    "position" : start, 
                    "variant-base" : var_alt_seq, 
                    "ins-corrected-position" : c_position, 
                    "var-length" : length }

            all_vars.append(var)

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

    if snps:
        grouped_snps = group_snps(snps)

        for start,length in grouped_snps:
            ins_correction = len([i for i in ins if i < start])

            c_position = start + ins_correction
            c_end = start + ins_correction + length

            var_ref_seq = rseq[c_position:c_end]
            var_alt_seq = qseq[c_position:c_end]

            if length > 1:
                var = {'type' : 'mnp'}
            
            else:
                var = {'type' : 'snp'}

            var.update({ "reference-base" : var_ref_seq,
                    "position" : start,
                    "variant-base" : var_alt_seq,
                    "ins-corrected-position" : c_position,
                    "var-length" : length })

            all_vars.append(var)

    if ins:
        for start,length in group_ins(ins):
            ins_correction = len([i for i in ins if i < start])
            
            c_position = start + ins_correction
            c_end = start + ins_correction + length

            var_ref_seq = rseq[c_position]
            var_alt_seq = qseq[start:start+length]

            if rseq[c_position:c_end] == var_alt_seq:
                # insertion but not in this sequence
                continue
            
            var = { "type" : "ins",
                    "reference-base" : var_ref_seq,
                    "position" : start,
                    "variant-base" : var_alt_seq,
                    "ins-corrected-position" : c_position,
                    "var-length" : length }
            
            all_vars.append(var)

    if len(all_vars) >= 1:
        fixed_vars = fix_overlapping_vars(all_vars, rseq, qseq)

        for var in fixed_vars:
           
            deconvolute_IUPAC(var)

            # Make 1-based
            ob_ins_c_pos = var.pop('ins-corrected-position')
            var['ob-ins-corrected-position'] = ob_ins_c_pos + 1 

            ob_pos = var.pop('position')
            var['one-based-reference-position'] = ob_pos + 1 

        # sort by positions
        fixed_vars_s = sorted(fixed_vars, key=lambda x: x['one-based-reference-position'])

        return fixed_vars_s

    else:
        return None

def type_variants(name, f_variants, variant_types):

    # { position: {variant dict}}
    formatted_vars_s = dict(zip([ i['one-based-reference-position'] for i in f_variants ], f_variants))

    sample_type = { 'sample_id' : name, 'typing' : [] }

    for vidx,vtype in enumerate(variant_types):

        summary = { i:vtype[i] for i in vtype if i not in [ 'variants', 'calling-definition' ]}

        sample_type['typing'].append( {'sample-typing-summary' : summary, 'sample-typing-result' : {} } )

        variant_lists = {}

        variant_lists['master'] = { 'variants' : vtype['variants'], 'calling-definition' : {} }

        for calling_id, calling_definition in vtype['calling-definition'].items():

            if 'variants' in calling_definition:
                # has own variants
                calling_def_dict = { i:calling_definition[i] for i in calling_definition if i != 'variants' }

                variant_lists[calling_id] = {'variants' : calling_definition['variants'], 
                                             'calling-definition' : { calling_id : calling_def_dict }}

            else:
                variant_lists['master']['calling-definition'].update( { calling_id : calling_definition } )
        
        for name,props in variant_lists.items():

            calls_keys = [
                'mutation_ref_calls', 
                'indel_ref_calls', 
                'mutation_calls', 
                'mutation_mixed_calls', 
                'indel_calls', 
                'no_calls'
                ]

            calls = dict.fromkeys(calls_keys, 0)
     
            for idx,var in enumerate(props['variants']):
                
                sample_var = formatted_vars_s.get(var['one-based-reference-position'])

                if sample_var:

                    variant_lists[name]['variants'][idx]['sample-call'] = ','.join(sample_var['variant-base'])
               
                    if sample_var['type'] == 'no-call':
                        variant_lists[name]['variants'][idx]['status'] = 'no-call'
                        calls['no_calls'] += 1

                    elif var['variant-base'] in sample_var['variant-base']:

                        if len(sample_var['variant-base']) > 1:
                            
                            variant_lists[name]['variants'][idx]['status'] = 'detect-mixed'

                            calls['mutation_mixed_calls'] += 1

                        elif len(sample_var['variant-base']) == 1:
                            
                            variant_lists[name]['variants'][idx]['status'] = 'detect'

                            if var['type'] == 'insertion' or var['type'] == 'deletion':
                                
                                calls['indel_calls'] += 1

                            else:
                                calls['mutation_calls'] += 1
                        
                                           
                else:
                    if var['reference-base'] == var['variant-base']:
                        variant_lists[name]['variants'][idx]['status'] = 'detect'

                        calls['mutation_calls'] += 1

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


def score_typing(sample_type):

    typing_summary = { 'sample_id' : sample_type['sample_id'] }  
    for vtype in sample_type['typing']:

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
        
        else:
            vtype['sample-typing-result']['variant-status'] = None

        typing_summary_header = vtype['sample-typing-summary']['phe-label'] + '|' + vtype['sample-typing-summary']['unique-id']

        typing_summary.update({ typing_summary_header : vtype['sample-typing-result']['variant-status'] } )

    return typing_summary, sample_type

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
    csvfilename = os.path.abspath(output_csv)

    fieldnames = list(typing_summary[0].keys())

    with open(csvfilename, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in typing_summary:
            writer.writerow(row)

def write_sample_variant_csv(name, variants, sample_csv_outdir, csv_N=False):

    csv_outpath = os.path.abspath(sample_csv_outdir)

    make_dir(csv_outpath)

    csvfilename = os.path.join(csv_outpath, normalise_fn(name) + '.csv')

    fieldnames = [
        'one-based-reference-position',
        'reference-base',
        'variant-base',
        'type'
    ]

    with open(csvfilename, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in variants:
            
            csv_row = { k: row[k] for k in fieldnames }
            csv_row['variant-base'] = '|'.join(csv_row['variant-base'])
            
            if csv_row['type'] == 'no-call':
                if csv_N == True:
                    writer.writerow(csv_row)
                else:
                    continue
            else:
                writer.writerow(csv_row)

def remove_terminal_gapns(seq):
    return re.sub(r'(N|-)*$', '', seq)

def go(args):

    refseq = get_ref_seq(args.msa, args.ref_name)

    if not refseq:
        print(args.ref_name + " not found\n") 
        sys.exit(1)

    variant_types = []
    for scheme in args.typing_yaml:
        variant_types.append(read_scheme_yaml(scheme))

    typing_summary = []
    with open(args.msa) as f:
        
        for name, seq, qual in readfq(f):
            if name == args.ref_name:
                continue
            else:
                variants = update_variants(remove_terminal_gapns(seq.upper()), refseq)

                print(name)
                if variants:
                    typed_variants = type_variants(name, variants, variant_types)
                    
                    sample_typing_summary, scored_variants = score_typing(typed_variants)

                    typing_summary.append(sample_typing_summary)

                    write_json(name, scored_variants, args.json_outdir, args.no_gzip_json)

                    write_sample_variant_csv(name, variants, args.sample_csv_outdir, args.csv_N)

    write_variant_types(typing_summary, args.output_csv)
    
                    
                    

