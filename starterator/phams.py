import MySQLdb
from .database import DB, get_db
from .phamgene import new_PhamGene
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO
from Bio import SeqIO
from collections import Counter
from . import utils
import subprocess
import os
from .utils import StarteratorError



def get_pham_number(phage_name, gene_number):
    try:
        db = DB()
        results = db.query("SELECT pham.name \n\
            FROM gene JOIN pham ON gene.GeneID = pham.Gene \n\
            JOIN phage ON gene.PhageID = phage.PhageID \n\
            WHERE phage.Name LIKE %s AND gene.Name LIKE %s \n\
            ESCAPE '!'", (phage_name+"%", '%'+str(gene_number)))
        row = results[0]
        pham_no = row[0]
        return str(pham_no)
    except:
        raise StarteratorError("Gene %s of Phage %s not found in database!" % (gene_number, phage_name))

def get_pham_colors():
    db = DB()
    results = db.query("SELECT `name`, `color` from `pham_color`");
    pham_colors = {}
    for row in results:
        pham_colors[str(row[0])] = row[1]
    return pham_colors

class Pham(object):
    def __init__(self, pham_no, genes=None):
        self.pham_no = pham_no
        self.genes = self.get_genes()
        self.color = self.get_color()
        self.stats = {}
        self.file = ""
        if genes:
            for gene in genes:
                self.add(gene)
            whole = "All" if len(genes) > 1 else "One"
            self.file = "%s%s" % (genes[0].phage_name, whole)


    def get_genes(self):
        """
            Get the genes of the Phamily
        """
        results = get_db().query("SELECT `gene`.`GeneID`, `gene`.`phageID`, " +
               " `Length`, `Start`, `Stop`, `Orientation`" +
               " FROM `gene`"
                " JOIN `pham` ON `gene`.`GeneID` = `pham`.`GeneID`" + 
               " WHERE `pham`.`name` =%s; ", self.pham_no)
        genes = {}
        for gene_info in results:
            gene_id = gene_info[0]
            phage_id = gene_info[1]
            start = gene_info[3]
            stop = gene_info[4]
            orientation = gene_info[5]
            gene = new_PhamGene(gene_id, start, stop, orientation, phage_id)
            genes[gene.gene_id] = gene
        if len(genes) < 1:
            raise StarteratorError("Pham Number %s not found!" % self.pham_no)
        return genes
    
    def get_phage_genes():
        pass

    def add(self, gene):
        """
            Add an unphameratored gene to the pham
        """
        self.genes[gene.gene_id] = gene

    def get_color(self):
        """
            Get the color of the phamily from the database
        """
        try:
            result = get_db().get("SELECT `name`, `color`\n\
                FROM `pham_color` WHERE `name` = %s;", self.pham_no)
            return result[1]
        except:
            raise StarteratorError("Pham number %s not found in database!" % self.pham_no)

    def add_alignment(self, alignment):
        """
            Using the alignment, add the alignment to the each gene in the pham
        """ 
        for record in alignment:
            gene = self.genes[record.id]
            gene.alignment = record
            gene.add_alignment_start_site()
            gene.add_alignment_candidate_starts()
            gene.add_gaps_as_features()

    def call_clustal(self, fasta_file):
        subprocess.check_call(['clustalw', 
        '-infile=%s' % (fasta_file),
        '-quicktree'])

    def make_fasta(self, file_name=None):
        if file_name == None:
            file_name = os.path.join(utils.INTERMEDIATE_DIR, "%sPham%s" % (self.file, self.pham_no))
        genes = [gene.sequence for gene in list(self.genes.values())]
        count = SeqIO.write(genes, "%s.fasta" % file_name, "fasta")

    def align(self):
        """
            Makes a fasta file of the genes in the Pham
            if the alignment already exists, uses that .aln as the alignment
            Otherwise, calls Clustalw from the command line and creates alignment
        """
        # files?
        file_name = os.path.join(utils.INTERMEDIATE_DIR, "%sPham%s" % (self.file, self.pham_no))
        genes = [gene.sequence for gene in list(self.genes.values())]
        count = SeqIO.write(genes, "%s.fasta" % file_name, "fasta")
        if len(self.genes) == 1:
            alignment = [gene.sequence]
        else:
            try:
                alignment = AlignIO.read(file_name+".aln", "clustal")
            except:
                # cline =  ClustalwCommandline("clustalw", infile=("%s.fasta" % file_name))
                # cline()
                self.call_clustal(file_name+".fasta")
                alignment = AlignIO.read(file_name+".aln", "clustal")
        self.add_alignment(alignment)

    def add_total_possible_starts(self):
        """ Returns a list of all the candidate starts from the alignment
        """
        self.total_possible_starts = []
        for gene in list(self.genes.values()):
            for site in gene.alignment_candidate_starts:
                if site not in self.total_possible_starts:
                    self.total_possible_starts.append(site) 
        self.total_possible_starts = sorted(self.total_possible_starts)
        return self.total_possible_starts

    def group_similar_genes(self):
        """
            Groups genes that have the same called start site, the same candidate starts
            and the same alignment (gaps are the same) together
        """
        groups = []
        i = 0
        genes = list(self.genes.values())
        grouped = [False for gene in genes]
        while i < len(self.genes):
            if not grouped[i]:
                # gene is not in a group yet
                gene = genes[i]
                j = i + 1 # genes before index i have been grouped
                group = []
                group.append(gene)  # add gene to this group - first one
                while j < len(self.genes): # see if other genes are similar 
                    if not grouped[j]: # skip genes that have already been grouped 
                        gene_2 = genes[j]
                        if gene.is_equal(gene_2): # if similar, then add to the group
                            grouped[i] = True 
                            grouped[j] = True
                            group.append(gene_2) 
                    j += 1
                groups.append(group)
            i += 1
        return groups

    def find_most_common_start(self, ignore_draft=False):
        """
            From the total candidate strats of each gene in the pham and all the start
            called in each gene, finds the start that is most commonly called.
            Returns a dictionary containing:
                "most_called" : a list of genes currently call the "most common start"
                "most_not_called" : a list of genes that have the "most common start" but do not call it
                "no_most_called" : a list of genes that do no have the "most called start"
                "possible" : a list containing lists of genes with the start of the index of self.total_possible_starts
                "called_start: a list containing lists of genes with the called start of the index of self.total_possible_starts

            Also, for each gene in the pham, a suggested start is given, gene.suggested_start["most_commom"]
            For genes that have the most common start called (or not) a tuple containing the index of the
            most common start and the coordinate of the sequence is given.
            For genes that do not have the most common start, a list of all possible starts, containing the index
            (useful when looking at the graphical output), and the coordinate is given.
        """
        # TODO:
            # add functionality for ignoring DRAFT phages?
        all_start_sites = [gene.alignment_start_site for gene in list(self.genes.values())]
        all_start_sites_set = set([gene.alignment_start_site for gene in list(self.genes.values())])
        start_stats = {}
        # creates two lists each containing a list of gene ids
        # for each candidate start of the pham:
        # start_stats["possible"] contains a list of genes with the candidate starts
            # for the index of each start in the pham
        # start_stats["called_starts"] contains of list of the genes that have the site
        #   of the index called as their start
        start_stats["possible"] = {}
        start_stats["called_starts"] = {}
        # start_stats["most_called"] = {}
        self.add_total_possible_starts()
        for i, site in enumerate(self.total_possible_starts):
            start_stats["possible"][i+1] = []
            # start_stats["most_called"][i+1] = []
            start_stats["called_starts"][i+1] = []
            for gene in list(self.genes.values()):
                if site in gene.alignment_candidate_starts:
                    start_stats["possible"][i+1].append(gene.gene_id)
                if site == gene.alignment_start_site:
                    start_stats["called_starts"][i+1].append(gene.gene_id)

        all_starts_count = Counter(all_start_sites)
        called_starts_count = all_starts_count.most_common()
        most_called_start_index = self.total_possible_starts.index(called_starts_count[0][0])+1
        genes_start_most_called = start_stats["called_starts"][most_called_start_index]
        start_stats["most_called_start"] = most_called_start_index
        # start_stats["most_called"] = start_stats["called_starts"][most_called_start_index]
        start_stats["most_called"] = []
        start_stats["most_not_called"] = []
        start_stats["no_most_called"] = []
        genes_without_most_called = []
        print(genes_start_most_called)
        for gene in list(self.genes.values()):
            if gene.gene_id in start_stats["possible"][most_called_start_index]:
                if gene.gene_id in genes_start_most_called:
                    if gene.orientation == 'F':   #only +1 for forward genes
                        gene.suggested_start["most_called"] =(most_called_start_index, gene.start+1)
                    else:
                        gene.suggested_start["most_called"] =(most_called_start_index, gene.start)
                    start_stats["most_called"].append(gene.gene_id)
                else:
                    start_stats["most_not_called"].append(gene.gene_id)
                    most_called_alignment_index = self.total_possible_starts[most_called_start_index-1]
                    suggested_start = gene.alignment_index_to_coord(most_called_alignment_index) #+1 issue dealt with in function
                    gene.suggested_start["most_called"] = (most_called_start_index, suggested_start)

            else:
                start_stats["no_most_called"].append(gene.gene_id)
                possible_starts_coords = []
                for start in gene.alignment_candidate_starts:
                    index = self.total_possible_starts.index(start) + 1
                    new_start = gene.alignment_index_to_coord(start) +1
                    possible_starts_coords.append((index, new_start))
                gene.suggested_start["most_called"] = possible_starts_coords


        self.stats["most_common"] = start_stats
        return start_stats

