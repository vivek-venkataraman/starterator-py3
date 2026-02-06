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
from .database import DB, get_db
from Bio.Blast import NCBIXML
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
import re
from itertools import groupby
from . import utils
from .utils import StarteratorError, clean_up_files
import subprocess
import math
import os


def get_protein_sequences():
    proteins = []
    results = get_db().query('SELECT GeneID, Translation from gene')
    for row in results:
        gene_id = row[0].replace("-", "_")
        translation = utils.decode_if_bytes(row[1])
        protein = SeqRecord(Seq(translation.replace('-', '')),
                            id=gene_id+"_", name=row[0], description=gene_id)
        proteins.append(protein)
    return proteins


def update_protein_db():
    clean_up_files(utils.INTERMEDIATE_DIR)
    proteins = get_protein_sequences()
    try:
        fasta_file = os.path.join(utils.PROTEIN_DB, "Proteins.fasta")
        SeqIO.write(proteins, fasta_file, 'fasta')
    except:
        print("creating proteins folder in correct place")
        utils.create_folders()
        fasta_file = os.path.join(utils.PROTEIN_DB, "Proteins.fasta")
        SeqIO.write(proteins, fasta_file, 'fasta')

    blast_db_command =  [
           'makeblastdb',
           '-in', fasta_file,
           '-dbtype', 'prot',
           '-title', 'Proteins',
           '-out', fasta_file
    ]
    # print blast_db_command
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
    # print new_count
    if int(new_count) != int(count):
        update_protein_db()
        config = utils.get_config()
        config["count"] = new_count
        utils.write_to_config_file(config)


def get_pham_no(phage_name, gene_number):
    """
        Gets the pham number of a gene, given the phage name and the gene number
    """
    # print phage_name, gene_number
    db = DB()
    query = "SELECT pham.Name \n\
            FROM gene JOIN pham ON gene.GeneID = pham.GeneID \n\
            JOIN phage ON gene.PhageID = phage.PhageID \n\
            WHERE (phage.Name LIKE %s or phage.PhageID = %s) AND gene.Name RLIKE %s \n\
            "% (phage_name + "%", phage_name, '^[:alpha:]*(_)*%s$' % str(gene_number))
    # print query
    try:
        results = db.query("SELECT pham.Name \n\
            FROM gene JOIN pham ON gene.GeneID = pham.GeneID \n\
            JOIN phage ON gene.PhageID = phage.PhageID \n\
            WHERE (phage.Name LIKE %s or phage.PhageID = %s) AND gene.geneid RLIKE %s",
            (phage_name + "%", phage_name, '^([[:alnum:]]*_)*([[:alpha:]])*%s$' % str(gene_number)))
        # print "DB query 1"
        if len(results) < 1:
            # print "DB query 1 failed, try search 2"
            results = db.query("SELECT pham.Name \n\
                FROM gene JOIN pham ON gene.GeneID = pham.GeneID \n\
                JOIN phage ON gene.PhageID = phage.PhageID \n\
                WHERE (phage.Name LIKE %s or phage.PhageID = %s) AND gene.geneID RLIKE %s",
                (phage_name + "%", phage_name, '^([[:alnum:]]*_)*([[:alpha:]])*%s$' % str(gene_number)))
        if len(results) < 1:
            #try to determine root of gene names since they are
            # print "DB query 2 failed, try search 3"
            results = db.query("SELECT pham.Name \n\
                FROM gene JOIN pham ON gene.GeneID = pham.GeneID \n\
                JOIN phage ON gene.PhageID = phage.PhageID \n\
                WHERE gene.geneid LIKE %s AND gene.geneID RLIKE %s",
                   (phage_name + "%", '^([[:alnum:]]*_)*([[:alpha:]])*%s$' % str(gene_number)))

        # print results
        row = results[0]
        pham_no = row[0]
        return str(pham_no)
    except:
        raise StarteratorError("Gene %s of Phage %s not found in database!" % (gene_number, phage_name))


def find_upstream_stop_site(start, stop, orientation, phage_sequence):
    """
        Given the coordinates of a gene, the sequence of the phage it is in, and 
        the orientation of the gene, returns a sequence that contains the gene
        and upstream sequence before a stop site.
    """
    ahead_of_start = 0
    stop_site_found = False
    stop_codons = ['AGT', 'AAT', 'GAT']
    while not stop_site_found:
        ahead_of_start += 99
        if orientation == 'R':
            if start + ahead_of_start > len(phage_sequence):     # i.e. hit end of phage while looking for stop
                ahead_of_start = len(phage_sequence) - start   # start is zero based counting
                ahead_of_start = ahead_of_start - ahead_of_start % 3
                sequence = Seq(phage_sequence[stop:(start+ahead_of_start)])
                sequence = sequence.reverse_complement()
                return sequence, ahead_of_start

            sequence = Seq(phage_sequence[stop:(start+ahead_of_start)])
            sequence = sequence.reverse_complement()
            if stop < 400:
                return sequence, ahead_of_start
        else:
            if start < ahead_of_start:
                ahead_of_start = start - start % 3
                sequence = Seq(phage_sequence[(start-ahead_of_start):stop])
                return sequence, ahead_of_start
            if stop < start:
                end_sequence = phage_sequence[(start-ahead_of_start):]
                start_sequence = phage_sequence[:stop]
                sequence = Seq(end_sequence+start_sequence)
            else:
                sequence = Seq(phage_sequence[(start-ahead_of_start):stop])
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
   
    def gene_no(self):
        get_gene_number(self.name)


pham_genes = {}

PRINTED_BAD_STARTS = set()


def new_PhamGene(db_id, start, stop, orientation, phage_id, name, phage_sequence=None):
    if db_id is None:
        return UnPhamGene(db_id, start, stop, orientation, phage_id, phage_sequence)
    if pham_genes.get(db_id, True):
        pham_genes[db_id] = PhamGene(db_id, start, stop, orientation, phage_id, name)
    return pham_genes[db_id]


def get_gene_number(gene_name):
    """ Given a gene_name, returns the number of the gene
    """
    # NAMING IN THIS DATABASE DRIVES ME CRAZY!!!
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
    def __init__(self, db_id, start, stop, orientation, phage_id, name, pham_no=None):
        self.db_id = db_id
        self.gene_id = db_id
        self.phage_id = phage_id
        self.start = start
        self.stop = stop
        self.cluster = None
        self.subcluster = None
        self.cluster_hash = None
        self.locustag = None
        self.gene_no = name
        self.full_name = self.phage_id + "_" + self.gene_no

        self.status = None # 'draft' = auto-annotated, 'final' = final/approved, 'gbk' imported non Pitt phage

        if orientation == 'R':
            self.start_codon_location = stop
            self.stop_codon_location = start + 1
        else:
            self.start_codon_location = start + 1
            self.stop_codon_location = stop

        self.orientation = orientation
        self.pham_no = pham_no
        self.pham_size = None
        # self.translation
        self.ahead_of_start = None
        self.sequence = self.make_gene()
        self.candidate_starts = self.add_candidate_starts()

        # --- QC: adjacent-start clusters and "bad starts" (all-but-last in each cluster) ---
        self.adjacent_candidate_start_groups = self._find_adjacent_start_groups(self.candidate_starts)

        # Flatten clusters into a single "bad starts" list:
        # for each adjacent run [a,b,c], treat [a,b] as bad (drop last), then merge across runs.
        bad = []
        for grp in self.adjacent_candidate_start_groups:
            if len(grp) >= 2:
                bad.extend(grp[:-1])
        self.bad_adjacent_candidate_starts = sorted(set(bad))
        self.has_bad_adjacent_candidate_starts = bool(self.bad_adjacent_candidate_starts)

        if self.has_bad_adjacent_candidate_starts:
            print(
                f"[QC] bad adjacent starts (offsets) gene={getattr(self, 'gene_no', getattr(self, 'number', '?'))}: {self.bad_adjacent_candidate_starts}")

        self.alignment = None
        self.alignment_start_site = None
        self.alignment_candidate_starts = None
        self.alignment_candidate_start_nums = None
        self.alignment_candidate_start_counts = None
        self.alignment_annot_start_nums = None
        self.alignment_annot_start_counts = None
        self.alignment_annot_start_fraction = None
        self.alignment_annot_counts_by_start = {}
        self.alignment_start_num_called = None
        self.alignment_start_conservation = None
        self.calls_most_annotated = None
        self.has_most_annotated = None
        self.suggested_start = {}

    def make_gene(self):
        """
           makes the gene which is a SeqRecord from Biopython. In this case the "gene" should
           include all the sequence upstream of the annotated start all the way to the first
           in frame stop codon.
        """
        phage = new_phage(phage_id=self.phage_id)
        self.phage_name = phage.get_name()
        self.name = self.phage_name + "_" + self.gene_no
        self.cluster = phage.cluster
        if self.cluster is None:
            self.cluster = "singleton"
        if phage.subcluster:
            self.subcluster = phage.subcluster
        else:
            self.subcluster = self.cluster
        self.cluster_hash = sum([pow(ord(elem), i+1) for i, elem in enumerate(self.subcluster)])
        status = phage.get_status()
        if status == 'final':        # values of 'draft' or 'gbk' considered draft quality by starterator
            self.draftStatus = False
        else:
            self.draftStatus = True

        phage_sequence = phage.get_sequence()
        if self.orientation == 'R':
            temp_start = self.stop
            self.stop = self.start
            self.start = temp_start
        sequence, self.ahead_of_start = find_upstream_stop_site(
                                self.start, self.stop, self.orientation, phage_sequence)
        self.ahead_of_start_coord = self.start - self.ahead_of_start
        gene = SeqRecord(sequence, id=self.gene_id, name=self.name,
                         description="|%i-%i| %s" % (self.start, self.stop, self.orientation))
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

    def _find_adjacent_start_groups(self, candidate_starts):
        """Returns groups of start sites that are adjacent in the same ORF (exactly 3 apart)

        groups neighboring start sites, ignoring ones that aren't adjacent

        bad_starts should NOT be called. Starts where there is at least 1 start codon immediately following it.
        """

        starts_sorted = sorted(candidate_starts)
        bad_groups = []
        current = [starts_sorted[0]]

        for s in starts_sorted[1:]:
            if s - current[-1] == 3:
                current.append(s)
            else:
                if len(current) >= 2:
                    bad_groups.append(current)
                current = [s]

        if len(current) >= 2:
            bad_groups.append(current)
        return bad_groups



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
        count = -1  # starts at -1 because the count starts at 0
        aligned_starts = []
        for index, char in enumerate(self.alignment.seq):
            if char != '-':
                count += 1
            if count in self.candidate_starts and char != '-':
                aligned_starts.append(index)
        self.alignment_candidate_starts = aligned_starts
        return aligned_starts

    def add_alignment_start_stats(self, pham):
        annotated = [gene.full_name for gene in pham.stats['most_common']['annot_list']]
        self.alignment_candidate_start_nums = []
        self.alignment_candidate_start_counts = []
        self.alignment_annot_start_nums = []
        self.alignment_annot_start_counts = []
        self.alignment_start_conservation = []

        if self.pham_no is None:
            self.pham_no = pham.pham_no

        self.pham_size = len(pham.genes)

        num_gene_in_pham = len(pham.genes)

        for startnum, genelist in pham.stats['most_common']['possible'].items():
            if self.full_name in genelist:
                self.alignment_candidate_start_nums.append(startnum)

        for num in self.alignment_candidate_start_nums:
            conservation_count = len(pham.stats['most_common']['possible'][num])
            self.alignment_candidate_start_counts.append(conservation_count)

            conserved_fraction = float(conservation_count) / float(num_gene_in_pham)
            self.alignment_start_conservation.append(conserved_fraction)

            annot_count = 0
            for gene in pham.stats['most_common']['called_starts'][num]:
                if gene in annotated:
                    annot_count += 1

            if annot_count > 0:
                self.alignment_annot_start_nums.append(num)
                self.alignment_annot_start_counts.append(annot_count)

        for startnum, genelist in pham.stats['most_common']['called_starts'].items():
            if self.full_name in genelist:
                self.alignment_start_num_called = startnum

        if len(self.alignment_annot_start_counts) > 0:
            most_annot_count = max(self.alignment_annot_start_counts)
            most_annot_index = self.alignment_annot_start_counts.index(most_annot_count)
            most_annotated_start_num = self.alignment_annot_start_nums[most_annot_index]
        else:
            most_annotated_start_num = None

        if most_annotated_start_num is None:
            self.calls_most_annotated = None
            self.has_most_annotated = None
        else:
            if self.alignment_start_num_called == most_annotated_start_num:
                self.calls_most_annotated = True
            else:
                self.calls_most_annotated = False

            if most_annotated_start_num in self.alignment_candidate_start_nums:
                self.has_most_annotated = True
            else:
                self.has_most_annotated = False

        total_annots = sum(self.alignment_annot_start_counts)
        self.alignment_annot_start_fraction = [float(count)/float(total_annots) for count in self.alignment_annot_start_counts]

        self.alignment_annot_counts_by_start = dict(zip(self.alignment_annot_start_nums, self.alignment_annot_start_counts))

        return

    def alignment_index_to_coord(self, index):
        """
                Given an index of the alignment
                finds the coordinates of the index on the phage sequence.
                The coordinate is 1 based count, not zero based
        """
        new_start_index = 0
        for i in range(0, index):
            if self.alignment.seq[i] != '-':
                new_start_index += 1
        if self.orientation == 'R':
            new_start_coords = (self.start + self.ahead_of_start - new_start_index)
        else:
            new_start_coords = (self.start - self.ahead_of_start + new_start_index + 1)
        return new_start_coords

    def add_gaps_as_features(self):
        # start by counting blocks of either bases or gap
        # use groupby() to give list of sizes of blocks of gap or sequence characters in block_length
        # and labels of type of block in block_type

        block_type = [k for k,g in groupby(self.alignment.seq, lambda x: 'seq' if x in ['A', 'C', 'G', 'T'] else 'gap')]
        block_length = [len(list(g)) for k, g in groupby(self.alignment.seq, lambda x: x in ['A', 'C', 'G', 'T'])]
        breakpoints = []
        found_start = False
        for i, length in enumerate(block_length):
            if i == 0:
                if length > self.alignment_start_site:
                    found_start = True
                    block_type.insert(0, 'seq')
                    breakpoints.append(self.alignment_start_site)
                breakpoints.append(length)
            else:
                if breakpoints[-1] + length > self.alignment_start_site and not found_start:
                    found_start = True
                    block_type.insert(i, 'seq')
                    breakpoints.append(self.alignment_start_site)
                    breakpoints.append(breakpoints[-2] + length)
                else:
                    breakpoints.append(breakpoints[-1] + length)

        for type_of_block, end_point in zip(block_type, breakpoints):
            end_point_index = breakpoints.index(end_point)
            if end_point_index == 0:
                start_point = 0
            else:
                start_point = breakpoints[end_point_index - 1]

            seq_feature = SeqFeature(FeatureLocation(start_point, end_point), type=type_of_block)
            self.alignment.features.append(seq_feature)

    def has_valid_start(self):
        return self.ahead_of_start in self.candidate_starts

    def is_equal(self, other):
        """
            Checks if another PhamGene is equal to this one
            PhamGenes are equal if they have the same called start, if the amount ahead of start
            (amount of sequence before the previous stop site) is the same, if the candidate
            starts of the genes are the same, and if the alignment gaps or not are the same
            (This is essentially, they would look the same on the graph output)

        """
        if self.alignment_start_site != other.alignment_start_site:
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
            print("phamgene.is_equal comparing features")
            if feature1.location.start != feature2.location.start:
                return False
            if feature1.location.end != feature2.location.end:
                return False
            if feature1.type != feature2.type:
                return False
        return True

    def get_locustag(self):
        db_return = get_db().get(
                "SELECT phage.annotationauthor, phage.status, gene.locustag from gene JOIN phage on gene.phageid=phage.phageid where gene.geneid = %s",
                self.db_id)

        self.annot_author = db_return[0]
        self.status = db_return[1]
        self.locustag = db_return[2]
        return

#test
    def __repr__(self):
        return 'Phamgene for %s' % self.gene_id


class UnPhamGene(PhamGene):
    def __init__(self, number, start, stop, orientation, phage_name, phage_sequence):
        self.number = number
        self.phage_name = phage_name
        self.gene_id = "%s_%s" % (phage_name, number)
        self.full_name = self.gene_id
        self.name = number
        self.start = start-1
        self.stop = stop
        self.orientation = orientation
        self.pham_size = None
        self.pham_no = None
        self.cluster = "Unassigned"
        self.subcluster = "Unassigned"
        self.cluster_hits = None
        self.subcluster_hits = None
        self.cluster_hash = sum([pow(ord(elem), i + 1) for i, elem in enumerate(self.subcluster)])


        if orientation == 'R':
            self.start_codon_location = stop
            self.stop_codon_location = start
        else:
            self.start_codon_location = start
            self.stop_codon_location = stop

        self.sequence = self.make_gene(phage_sequence)
        self.candidate_starts = self.add_candidate_starts()

        # --- QC: adjacent-start clusters and "bad starts" (all-but-last in each cluster) ---
        self.adjacent_candidate_start_groups = self._find_adjacent_start_groups(self.candidate_starts)

        bad = []
        for grp in self.adjacent_candidate_start_groups:
            if len(grp) >= 2:
                bad.extend(grp[:-1])
        self.bad_adjacent_candidate_starts = sorted(set(bad))
        self.has_bad_adjacent_candidate_starts = bool(self.bad_adjacent_candidate_starts)





# test change

        if self.has_bad_adjacent_candidate_starts:
            print(
                f"[QC] bad adjacent starts (offsets) gene={getattr(self, 'gene_no', getattr(self, 'number', '?'))}: {self.bad_adjacent_candidate_starts}")

        self.alignment = None
        self.alignment_start = None
        self.alignment_candidate_starts = None
        self.alignment_candidate_start_nums = None
        self.alignment_annot_start_nums = None
        self.alignment_annot_start_counts = None
        self.alignment_start_num_called = None
        self.calls_most_annotated = None
        self.has_most_annotated = None
        self.suggested_start = {}
        self.draftStatus = True

    def make_gene(self, phage_sequence):
        if self.orientation == 'R':
            temp_start = self.stop
            self.stop = self.start
            self.start = temp_start
        sequence, self.ahead_of_start = find_upstream_stop_site(
                                self.start, self.stop, self.orientation, phage_sequence)
        gene = SeqRecord(sequence, id=self.gene_id, name=self.gene_id)
        return gene

    def phambymatch(self):
        #try to must make a perfect match to the translation field
        protein = str(self.sequence[self.ahead_of_start:].seq.translate())
        #repair translations if start codon was TTG or GTG and remove stop codon
        protein = "M" + protein[1:-1]
        db = DB()
        result = db.query("SELECT GeneID FROM gene WHERE gene.Translation = %s", protein)
        if len(result) < 1:
            return None
        else:
            result2 = db.query("SELECT phamid FROM gene WHERE geneid = %s", result[0])
            print("pham %s by exact match to gene %s"%(result2[0],result[0]))
            number, = result2[0]
            self.pham_no = number

            return self.pham_no

    def blast(self):
        # not sure where to put this... this makes more sense, 
        # but I wanted to keep the Genes out of file making...
        # print "Running BLASTp"
        try:
            result_handle = open("%s/%s.xml" % (utils.INTERMEDIATE_DIR, self.gene_id))
            result_handle.close()
        except:
            protein = SeqRecord(self.sequence[self.candidate_starts[0]:].seq.translate(), id=self.gene_id)
            # print protein, self.sequence
            # short proteins need lower e_value
            query_len = (self.stop - self.start) / 3
            if query_len < 50:
                e_value = math.pow(10, -5)
            else:
                e_value = math.pow(10, -20)

            SeqIO.write(protein, '%s/%s.fasta' % (utils.INTERMEDIATE_DIR, self.gene_id), 'fasta')
            # Using subprocess approach instead of deprecated Bio.Application
            blast_args = ["%sblastp" % utils.BLAST_DIR,
                          "-out", '%s/%s.xml' % (utils.INTERMEDIATE_DIR, self.gene_id),
                          "-outfmt", "5",
                          "-query", '%s/%s.fasta' % (utils.INTERMEDIATE_DIR, self.gene_id),
                          "-db", "\"%s/Proteins.fasta\"" % (utils.PROTEIN_DB),
                          "-evalue", str(e_value)
                          ]
            # print " ".join(blast_args)
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
            # print first_result
            if "_" not in first_result:
                first_result = blast_record.descriptions[1].title.split(',')[0].split(' ')[-1]
            # Try to get pham directly from gene name
            db = DB()
            results = db.query("SELECT phamid from gene where geneID = %s", first_result)
            if len(results) == 1:
                number, = results[0]
                self.pham_no = number
                return number
            else:
                phage_name = first_result.split("_")[0]
                #exception for error in locus tags specific to this phage:
                if phage_name == 'FRIAPREACHER':
                    phage_name = 'FRIARPREACHER'
                if phage_name.lower() == "draft":
                    phage_name = first_result.split("_")[-3]
                gene_number = first_result.split("_")[-1]
                # print phage_name, gene_number
                pham_no = get_pham_no(phage_name, gene_number)
                self.pham_no = pham_no
                return pham_no
        else:
            self.pham_no = None
            return None

    def add_cluster_hits(self):
        self.cluster_hits = []
        db = DB()
        result = db.query("SELECT distinct(phage.cluster) FROM phage JOIN gene on gene.phageid = phage.phageid\
                           WHERE gene.phamid = %s", self.pham_no)

        for item in result:
            hit = list(item)[0]
            if hit is not None:
                self.cluster_hits.append(hit)

        self.subcluster_hits = []
        result2 = db.query("SELECT distinct(phage.subcluster) FROM phage JOIN gene on gene.phageid = phage.phageid\
                           WHERE gene.phamid = %s", self.pham_no)
        for item in result2:
            hit2 = list(item)[0]
            if hit2 is not None:
                self.subcluster_hits.append(hit2 )

