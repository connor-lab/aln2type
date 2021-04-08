import csv
from Bio import SeqIO
from Bio.Seq import Seq

"""
   Much of this from https://github.com/rpetit3/vcf-annotator
"""

class Annotator(object):
    """Annotate a given VCF file according to the reference GenBank."""

    def __init__(self, gb_file=False, vcf_records=False):
        """Initialize variables."""
        self.__annotated_features = ["mat_peptide", "CDS", "tRNA", "rRNA", "ncRNA", "misc_feature"]
        self.__gb = GenBank(gb_file)
        self.__gb.accession = self.get_accession_from_gbk(gb_file)
        self.__vcf_records = self.fix_codon_mnp(vcf_records)

    def get_accession_from_gbk(self, gb_file):

        accession = []

        with open(gb_file, 'r') as gb_fh:
            for record in SeqIO.parse(gb_fh, 'genbank'):
                accession.append(record.name)

        if len(set(accession)) > 1:
            raise NotImplementedError("Found more than one record in your GenBank file - multi-GenBank files aren't supported by aln2type")

        return accession[0]


    def fix_codon_mnp(self, vcf_records):
        fixed_records = []
        for record in vcf_records:
            
            self.__gb.index = record['one-based-reference-position']
            
            if record['type'] == 'mnp' and self.__gb.feature_exists:

                # Determine codon information
                codon = self.__gb.codon_by_position(record['one-based-reference-position'])

                if codon[1] != 0 or record['var-length'] > 3:
                    
                    codon_aligned_record = { 
                        'type': 'mnp',
                        'reference-base': codon[0][:codon[1]] + record['reference-base'],
                        'variant-base': codon[0][:codon[1]] + record['variant-base'],
                        'var-length': codon[1] + record['var-length'],
                        'ob-ins-corrected-position': record['ob-ins-corrected-position'] - codon[1], 
                        'one-based-reference-position': record['one-based-reference-position'] - codon[1],
                        'iupac-variant-bases': [ codon[0][:codon[1]] + i for i in record['iupac-variant-bases'] ]
                        }

                    for codon in range(0, codon_aligned_record['var-length'], 3):
                        new_ref = codon_aligned_record['reference-base'][codon:codon+3]
                        new_var = codon_aligned_record['variant-base'][codon:codon+3]

                        new_var_length = len(new_var)
                        
                        if new_var_length == 1:
                            new_type = 'snp'
                        else:
                            new_type = 'mnp'
                        
                        new_one_based_reference_position = codon_aligned_record['one-based-reference-position'] + codon

                        new_ob_ins_corrected_position = codon_aligned_record['ob-ins-corrected-position'] + codon

                        new_iupac_variant_bases = [ i[codon:codon+3] for i in codon_aligned_record['iupac-variant-bases'] ]

                        new_record = { 
                        'type': new_type,
                        'reference-base': new_ref,
                        'variant-base': new_var,
                        'var-length': new_var_length,
                        'ob-ins-corrected-position': new_ob_ins_corrected_position, 
                        'one-based-reference-position': new_one_based_reference_position,
                        'iupac-variant-bases': new_iupac_variant_bases
                        }

                        fixed_records.append(new_record)

                        
                else:
                    fixed_records.append(record)

            else:
                fixed_records.append(record)

         

        return fixed_records

    def annotate_records(self):
        """Annotate each record acording to the input GenBank."""

        updated_records = []
        for record in self.__vcf_records[:]:
            if record['type'] == 'no-call' or record['type'] == 'mnp-snp':
                continue

            if record['type'] == 'del' or record['type'] == 'ins':
                self.__gb.index = int(record['one-based-reference-position'] + 1)

            else:
                self.__gb.index = int(record['one-based-reference-position'])

            print(record)
            # Set defaults
            record['RefCodon'] = None
            record['AltCodon'] = None
            record['RefAminoAcid'] = None
            record['AltAminoAcid'] = None
            record['CodonPosition'] = None
            record['SNPCodonPosition'] = None
            record['AminoAcidChange'] = None
            record['IsSynonymous'] = 9
            record['Comments'] = None
            record['IsGenic'] = '0'
            record['IsPseudo'] = '0'
            record['LocusTag'] = None
            record['Gene'] = None
            record['Inference'] = None
            record['Product'] = None
            record['ProteinID'] = None
            record['FeatureType'] = 'inter_genic'

            # Get annotation info
            if self.__gb.feature_exists:
                record['FeatureType'] = self.__gb.feature.type
                if self.__gb.feature.type in self.__annotated_features:
                    feature = self.__gb.feature
                    if feature.type == "CDS" or feature.type == 'mat_peptide':
                        record['IsGenic'] = '1'

                    record.update( 
                        {
                        'LocusTag': ';'.join(self.__gb.feature.qualifiers['locus_tag']),
                        'Gene': ';'.join(self.__gb.feature.qualifiers['gene']), 
                        'Product': ';'.join(self.__gb.feature.qualifiers['product']),
                        'ProteinID': ';'.join(self.__gb.feature.qualifiers['protein_id'])
                        }
                    )

                    if feature.type == "tRNA":
                        qualifiers['Note'] = 'anticodon'
                    
                    if 'pseudo' in feature.qualifiers:
                        record['IsPseudo'] = '1'

            # Support codon replacements
            if record['type'] != 'snp' and len(str(record['variant-base'])) == 3 and len(record['reference-base']) == 3:
                is_codon_mnp = True
            else:
                is_codon_mnp = False

            # Determine variant type
            if record['type'] == 'del':
                record['VariantType'] = 'Deletion'
            elif record['type'] == 'ins':
                record['VariantType'] = 'Insertion'
            else:

                if len(record['iupac-variant-bases']) > 1:
                    record['VariantType'] = 'Ambiguous_SNP'
                else:
                    record['VariantType'] = 'SNP'

                if int(record['IsGenic']):
                    alt_base = str(record['variant-base'])

                    # Determine codon information
                    codon = self.__gb.codon_by_position(int(record['one-based-reference-position']))
                    record['RefCodon'] = ''.join(list(codon[0]))
                    record['SNPCodonPosition'] = codon[1]
                    record['CodonPosition'] = codon[2]

                    if is_codon_mnp:
                        if record['SNPCodonPosition'] != 0:
                            raise ValueError("Do not support non-codon-aligned multinucleotide polymorphisms")

                    # Adjust for ambiguous base and negative strand.
                    if feature.strand == -1:
                        alt_base = str(
                            Seq(alt_base).complement()
                        )

                        record['Comments'] = 'Negative:{0}->{1}'.format(
                            Seq(record['reference-base']).complement(),
                            alt_base
                        )

                    # Determine alternates
                    record['AltCodon'] = list(record['RefCodon']) 

                    if not is_codon_mnp:
                        for idx, base in enumerate(alt_base):
                            record['AltCodon'][idx + codon[1]] = alt_base[idx]

                    else:
                        record['AltCodon'] = alt_base

                    record['AltCodon'] = ''.join(record['AltCodon'])
                    record['RefAminoAcid'] = Seq(
                        record['RefCodon']
                    ).translate()
                    record['AltAminoAcid'] = Seq(
                        record['AltCodon']
                    ).translate()
                    record['AminoAcidChange'] = '{0}{1}{2}'.format(
                        str(record['RefAminoAcid']),
                        record['CodonPosition'],
                        str(record['AltAminoAcid'])
                    )

                    if record['VariantType'] != 'Ambiguous_SNP':
                        ref = str(record['RefAminoAcid'])
                        alt = str(record['AltAminoAcid'])
                        if ref == alt:
                            record['IsSynonymous'] = 1
                        else:
                            record['IsSynonymous'] = 0

    def get_vcf_records(self):
        return self.__vcf_records


class GenBank(object):
    """A class for parsing GenBank files."""

    def __init__(self, gb=False):
        """Inititalize variables."""
        self.genbank_file = gb
        self.records = {}
        self.record_index = {}
        self.__gb = None
        self._index = None
        self._accession = None
        self.__position_index = None
        self.feature = None
        self.features = ["mat_peptide", "CDS", "rRNA", "tRNA", "ncRNA", "repeat_region",
                         "misc_feature"]
        self.gene_codons = {}
        self.parse_genbank()

    @property
    def accession(self):
        """Accession for records."""
        return self._index

    @accession.setter
    def accession(self, value):
        self._accession = value
        self.__gb = self.records[value]
        self.__position_index = self.record_index[value]

    @property
    def index(self):
        """Postion index for features."""
        return self._index

    @index.setter
    def index(self, value):
        self._index = self.__position_index[value - 1]
        self.__set_feature()

    def parse_genbank(self):
        with open(self.genbank_file, 'r') as gb_fh:
            for record in SeqIO.parse(gb_fh, 'genbank'):
                self.records[record.name] = record
                self.gene_codons[record.name] = {}
                self.record_index[record.name] = [None] * len(record.seq)
                for i in range(len(record.features)):
                    if record.features[i].type in self.features:
                        start = int(record.features[i].location.start)
                        end = int(record.features[i].location.end)
                        self.record_index[record.name][start:end] = [i] * (end - start)


    def __set_feature(self):
        if self._index is None:
            self.feature_exists = False
            self.feature = None
        else:
            self.feature_exists = True
            self.feature = self.records[self._accession].features[self._index]

    def codon_by_position(self, pos):
        """Retreive the codon given a postion of a CDS feature."""
        if self._index not in self.gene_codons[self._accession]:
            self.split_into_codons()
        gene_position = self.position_in_gene(pos)
        codon_position = gene_position // 3
        return [self.gene_codons[self._accession][self._index][codon_position],
                gene_position % 3,
                codon_position + 1]

    def split_into_codons(self):
        """Split the complete CDS feature in to a list of codons."""
        start = self.feature.location.start
        end = self.feature.location.end
        seq = ''.join(list(self.feature.extract(self.__gb.seq))) 

        if self.feature.strand == -1:
            seq = Seq(seq).reverse_complement()

        self.gene_codons[self._accession][self._index] = [
            seq[i:i + 3] for i in range(0, len(seq), 3)
        ]

    def position_in_gene(self, pos):
        """Return a codon postion in a gene."""
        positions = list(self.feature.location)
        if self.feature.strand == -1:
             raise ValueError('not yet tested')
        return positions.index(pos-1)

    def base_by_pos(self, pos):
        """Print the base by position."""
        print(self.__gb.seq[pos - 1])

    def is_transition(self, ref_base, alt_base):
        """
        Identify SNP as being a transition or not.

        1: Transition, 0:Transversion
        """
        substitution = ref_base + alt_base
        transition = ['AG', 'GA', 'CT', 'TC']

        if substitution in transition:
            return 1

        return 0
        