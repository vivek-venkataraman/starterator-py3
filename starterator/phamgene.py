# Copyright (c) 2013, 2014 All Right Reserved, Hatfull Lab, University of Pittsburgh
#
# THIS CODE AND INFORMATION ARE PROVIDED "AS IS" WITHOUT WARRANTY OF ANY
# KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A
# PARTICULAR PURPOSE.  USE AT YOUR OWN RISK.
#
# Marissa Pacey
# April 4, 2014
# Class and functions for pham genes

from .phage import new_phage
import copy
from .database import DB
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastpCommandline as Blastp
from Bio.Blast.Applications import BlastallCommandline
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
import re
from .database import get_db
from . import utils
from .utils import StarteratorError
import subprocess
import math
import os

def get_protein_sequences():
    proteins = []
    results = get_db().query('SELECT GeneID, translation from gene')
    for row in results:
        gene_id = row[0].replace("-", "_")
        protein = SeqRecord(Seq(row[1].replace('-', ''), IUPAC.protein), id=gene_id+"_", name=row[0], description=gene_id)
        proteins.append(protein)
    return proteins

def update_protein_db():
    proteins = get_protein_sequences()
    try:
        fasta_file = os.path.join(utils.PROTEIN_DB, "Proteins.fasta")
        count = SeqIO.write(proteins, fasta_file, 'fasta')
    except:
        print("creating proteins folder in correct place")
        utils.create_folders()
        fasta_file = os.path.join(utils.PROTEIN_DB, "Proteins.fasta")
        count = SeqIO.write(proteins, fasta_file, 'fasta')
    if True:
        blast_db_command = [utils.BLAST_DIR + 'makeblastdb',
                    '-in',"\""+ fasta_file+ "\"",
                    "-dbtype","prot", "-title", "Proteins",
                     "-out", "%s"% fasta_file]
        print(blast_db_command)
    # else:
    #     blast_db_command = [BLAST_DIR + 'formatdb',
    #                 '-i', "\""+ fasta_file+ "\"",
    #                 '-o', 'T',
    #                 "-t", "Proteins"]
    #     print blast_db_command
    subprocess.check_call(blast_db_command)

def check_protein_db(count):
    results = get_db().query('SELECT count(*) from gene')
    new_count = results[0][0]
    print(new_count)
    if int(new_count) != int(count):
        update_protein_db()
        config = utils.get_config()
        config["count"] = new_count
        utils.write_to_config_file(config)

def get_pham_no(phage_name, gene_number):
    """
        Gets the pham number of a gene, given the phage name and the gene number
    """
    print(phage_name, gene_number)
    db = DB()
    query = "SELECT pham.name \n\
            FROM gene JOIN pham ON gene.GeneID = pham.GeneID \n\
            JOIN phage ON gene.PhageID = phage.PhageID \n\
            WHERE (phage.Name LIKE %s or phage.PhageID = %s) AND gene.Name RLIKE %s \n\
            "% (phage_name+ "%", phage_name, '^[:alpha:]*(_)*%s$' % str(gene_number))
    print(query)
    try:
        results = db.query("SELECT pham.name \n\
            FROM gene JOIN pham ON gene.GeneID = pham.GeneID \n\
            JOIN phage ON gene.PhageID = phage.PhageID \n\
            WHERE (phage.Name LIKE %s or phage.PhageID = %s) AND gene.Name RLIKE %s", 
            (phage_name+"%", phage_name, '^([[:alnum:]]*_)*([[:alpha:]])*%s$' % str(gene_number)))
        print(results)
        row = results[0]
        pham_no = row[0]
        return str(pham_no)
    except:
        raise StarteratorError("Gene %s of Phage %s not found in database!" % (gene_number, phage_name))


def find_upstream_stop_site(start, stop, orientation, phage_sequence):
    """
        Given the coordinates of a gene, the sequence of the phage it is in, and 
        the orientation of the gene, returns a seqeunce that contains the gene 
        and upstream sequence before a stop site.
    """
    ahead_of_start = 0
    stop_site_found = False
    stop_codons = ['AGT', 'AAT', 'GAT']
    while not stop_site_found:
        ahead_of_start += 99
        if orientation == 'R':
            sequence = Seq(phage_sequence[stop-1:(start+ahead_of_start)],
                    IUPAC.unambiguous_dna)
            sequence = sequence.reverse_complement()
            if stop < 400:
                return sequence, ahead_of_start
        else:
            if start < ahead_of_start:
                ahead_of_start = start - start % 3
                sequence = Seq(phage_sequence[(start-ahead_of_start):stop],
                     IUPAC.unambiguous_dna)
                return sequence, ahead_of_start
            if stop < start:
                end_sequence = phage_sequence[(start-ahead_of_start):]
                start_sequence = phage_sequence[:stop]
                sequence = Seq(end_sequence+start_sequence, IUPAC.unambiguous_dna)
            else:
                sequence = Seq(phage_sequence[(start-ahead_of_start):stop],
                 IUPAC.unambiguous_dna)
        sequence_ahead_of_start = sequence[:ahead_of_start]
        sequence_ahead_of_start = sequence_ahead_of_start[::-1]
        
        for index in range(0, len(sequence_ahead_of_start), 3):
            codon = str(sequence_ahead_of_start[index:index+3])
            if codon in stop_codons:
                new_ahead_of_start = index
                new_sequence = sequence[(ahead_of_start - index):]
                return new_sequence, new_ahead_of_start


class Gene(object):
    def __init__(self, phage, name, start, stop, orientation, db_id=None):
        self.phage = phage
        self.name = name
        self.start = start
        self.stop = stop
        self.orientation = orientation
        self.db_id = db_id
   
    def gene_no():
        get_gene_number(self.name)


pham_genes = {}
def new_PhamGene(db_id, start, stop, orientation, phage_id, phage_sequence=None):
    if db_id == None:
        return UnPhamGene(db_id, start, stop, orientation, phage_id, phage_sequence)
    if pham_genes.get(db_id, True):
        pham_genes[db_id] = PhamGene(db_id, start, stop,
                                 orientation, phage_id)
    return pham_genes[db_id]

def get_gene_number(gene_name):
    """ Given a gene_name, returns the number of the gene
    """
    # NAMING IN THIS DATBASE DRIVES ME CRAZY!!!
    # GeneID in database: form of <PhageID>_(<PhageName>([_-]Draft*))*_(gene)*(gp)*<GeneNo>
    #   where PhageID can be a number or the name of the phage (with _Draft perhaps)
    
    # ...and I don't think this function is needed. Precisely because of this!
    # match = re.search(r'^(\w+)([_]*\w*)_([])')
    match = re.search(r'^((\w+)([_-]*\w*)_)*([a-zA-Z]*)([0-9]+)+$', gene_name)
    gene_number = match.groups()[-1]
    try:
        return int(gene_number)
    except:
        # it is one of the 3 horrible genes that do not have an actual number
        # !!! WHAT DO I DO HERE?????  doesn't look like it is terrible for it to be a string?
        # so it will return hypothetical or null - don't ask me why
        return gene_number


class PhamGene(Gene):
    def __init__(self, db_id, start, stop, orientation, phage_id, pham_no=None):
        self.db_id = db_id
        self.phage_id = phage_id
        self.start = start
        self.stop = stop
        self.orientation = orientation
        self.pham_no = pham_no
        # self.translation
        # self.ahead_of_start
        self.sequence = self.make_gene()
        self.candidate_starts = self.add_candidate_starts()
        self.alignment = None
        self.alignment_start = None
        self.alignment_candidate_starts = None
        self.suggested_start = {}


    def make_gene(self):
        """
           makes the gene 
        """
        phage = new_phage(phage_id=self.phage_id)
        self.phage_name = phage.get_name()
        gene_no = self.db_id.split("_")[-1]
        gene_no = gene_no.split(" ")[0]
        self.gene_id = self.phage_name + "_" + gene_no
        self.gene_id = self.gene_id.replace('-', "_")
        phage_sequence = phage.get_sequence()
        if self.orientation == 'R':
            temp_start = self.stop
            self.stop = self.start
            self.start = temp_start
        sequence, self.ahead_of_start = find_upstream_stop_site(
                                self.start, self.stop, self.orientation, phage_sequence)
        gene = SeqRecord(sequence, id=self.gene_id , name=self.gene_id,
                 description="|%i-%i| %s" %(self.start, self.stop, self.orientation))
        return gene



    def add_candidate_starts(self):
        """
            Finds all the possible start site of the gene and returns a list of indexes of start sites
        """
        gene_sequence = self.sequence.seq
        starts = []
        start_codons = ['ATG', 'GTG', 'TTG']
        for index in range(0, len(gene_sequence), 3):
            codon = str(gene_sequence[index:index+3])
            if codon in start_codons:
                starts.append(index)
        return sorted(starts)

    def add_alignment_start_site(self):
        """
            Gives the coordinate the called start site in the alignment sequence
        """
        count = -1
        i = 0
        for index, letter in enumerate(self.alignment.seq):
                if letter in ['A', 'G', 'T', 'C']:
                    count += 1
                    i = index
                    if count >= self.ahead_of_start:
                        break
        if count > self.ahead_of_start:
            i -= 1
        self.alignment_start_site = i
        return i

    def add_alignment_candidate_starts(self):
        """
            Creates a list of candidate starts of the alignment based on the candidate starts
            of the gene
        """
        count = -1 # starts at -1 because the count starts at 0
        aligned_starts = []
        for index, char in enumerate(self.alignment.seq):
            if char != '-':
                count += 1
            if count in self.candidate_starts and char != '-':
                aligned_starts.append(index)
        self.alignment_candidate_starts = aligned_starts
        return aligned_starts

    def alignment_index_to_coord(self, index):
        """
                Given an index of the alignment
                finds the coordinates of the index on the phage sequence.
        """
        new_start_index = 0
        for i in range(0, index):
            if self.alignment.seq[i] != '-':
                new_start_index += 1
        if self.orientation == 'R':
            new_start_coords = (self.start + 
                self.ahead_of_start - new_start_index)
        else:
            new_start_coords = (self.start - 
                self.ahead_of_start + new_start_index + 1)
        return new_start_coords

    def add_gaps_as_features(self):
            gap = False
            gap_count = 0
            seq_count = 0
            for index, char in enumerate(self.alignment.seq):
                if char == '-' and gap == False:
                    gap = True
                    gap_count = 0
                    if seq_count > 0:
                        seq_feature = SeqFeature(FeatureLocation(index-seq_count, index-1),
                            type='seq', strand=None)
                        self.alignment.features.append(seq_feature)
                elif char == '-' and gap == True:
                    gap_count += 1
                elif char != '-' and gap == True:
                    gap = False
                    seq_count = 0
                    gap_feature = SeqFeature(FeatureLocation(index-gap_count-1, index-1),
                            type='gap', strand=None)
                    self.alignment.features.append(gap_feature)
                else: #gap = False and char != '-'
                    seq_count +=1
            if gap == True:
                gap_feature = SeqFeature(FeatureLocation(index-gap_count-1, index), 
                            type='gap', strand=None) 
                self.alignment.features.append(gap_feature)
            else:
                seq_feature = SeqFeature(FeatureLocation(index-seq_count, index), 
                            type='seq', strand=None)
                self.alignment.features.append(seq_feature)

    def is_equal(self, other):
        """
            Checks if another PhamGene is equal to this one
            PhamGenes are equal if they have the same called start, if the amount ahead of start
            (amount of sequence before the previous stop site) is the same, if the candidate
            starts of the genes are the same, and if the alignment gaps or not are the same
            (This is essentially, they would look the same on the graph output)
        """
        if self.start != other.start:
            return False
        if self.ahead_of_start != other.ahead_of_start:
            return False
    
        if set(self.alignment_candidate_starts) != set(other.alignment_candidate_starts):
            return False
        if len(self.sequence.features) != len(other.sequence.features):
            return False
        self.sequence.features.sort()
        other.sequence.features.sort()
        for feature1, feature2 in zip(self.sequence.features, other.sequence.features):
            if feature1.location.start != feature2.location.start:
                return False
            if feature1.location.end != feature2.location.end:
                return False
            if feature1.type != feature2.type:
                return False
        return True


class UnPhamGene(PhamGene):
    def __init__(self, number, start, stop, orientation, phage_name, phage_sequence):
        self.number = number
        self.phage_name = phage_name
        self.gene_id = "%s_%s" % (phage_name, number)
        self.start = start-1
        self.stop = stop
        self.orientation = orientation
        self.sequence = self.make_gene(phage_sequence)
        self.candidate_starts = self.add_candidate_starts()
        self.alignment = None
        self.alignment_start = None
        self.alignment_candidate_starts = None
        self.suggested_start = {}


    def make_gene(self, phage_sequence):
        if self.orientation == 'R':
            temp_start = self.stop
            self.stop = self.start
            self.start = temp_start
        sequence, self.ahead_of_start = find_upstream_stop_site(
                                self.start, self.stop, self.orientation, phage_sequence)
        gene = SeqRecord(sequence, id=self.gene_id, name=self.gene_id)
        return gene

    def blast(self):
        # not sure where to put this... this makes more sense, 
        # but I wanted to keep the Genes out of file making...

        try:
            result_handle =  open("%s/%s.xml" % (utils.INTERMEDIATE_DIR, self.gene_id))
            result_handle.close()
        except:
            protein = SeqRecord(self.sequence.seq.translate(), id=self.gene_id)
            print(protein, self.sequence)
            e_value = math.pow(10, -30)
            SeqIO.write(protein, '%s/%s.fasta' % (utils.INTERMEDIATE_DIR, self.gene_id), 'fasta')
            blast_command = Blastp(
                            query='%s%s.fasta' % (utils.INTERMEDIATE_DIR, self.gene_id),
                            db="\"%s/\"" % (os.path.abspath(utils.PROTEIN_DB)), evalue=e_value, outfmt=5,
                            out="%s.xml" % (os.path.join(utils.INTERMEDIATE_DIR, self.gene_id)))
            # print self.gene_id, "\"%sProteins\"" % (utils.PROTEIN_DB)
            blast_args = ["%sblastp"  % utils.BLAST_DIR, 
                "-out", '%s/%s.xml' % (utils.INTERMEDIATE_DIR, self.gene_id),
                "-outfmt", "5",
                "-query", '%s/%s.fasta' % (utils.INTERMEDIATE_DIR, self.gene_id),
                "-db", "\"%s/Proteins.fasta\"" % (utils.PROTEIN_DB),
                "-evalue", str(e_value)
                ]
            print(" ".join(blast_args))
            try:
                subprocess.check_call(blast_args)
            except:
                raise StarteratorError("Blast could not run!")
        # print blast_command
        # stdout, stderr = blast_command()
        return self.parse_blast()

    def parse_blast(self):
        result_handle = open("%s/%s.xml" % (utils.INTERMEDIATE_DIR, self.gene_id))

        try:
            blast_record = NCBIXML.read(result_handle)
        except:
            result_handle.close()
            result_handle = open('%s/%s.xml' % (self.output_dir, self.name))
            blast_records = NCBIXML.parse(result_handle)
            blast_record = next(blast_records)

        if len(blast_record.descriptions) > 0:
            first_result = blast_record.descriptions[0].title.split(',')[0].split(' ')[-1]
            print(first_result)
            phage_name = first_result.split("_")[0]
            gene_number = first_result.split("_")[-1]
            print(phage_name, gene_number)
            pham_no = get_pham_no(phage_name, gene_number)
            self.pham_no = pham_no
            return pham_no
        else:
            self.pham_no = None
            return None
