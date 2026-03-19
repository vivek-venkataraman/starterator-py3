# Copyright (c) 2013, 2014 All Right Reserved, Hatfull Lab, University of Pittsburgh
#
# THIS CODE AND INFORMATION ARE PROVIDED "AS IS" WITHOUT WARRANTY OF ANY
# KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A
# PARTICULAR PURPOSE.  USE AT YOUR OWN RISK.
#
# Marissa Pacey
# April 4, 2014
# Making of the Starterator Reports

import subprocess
import pickle
import PyPDF2
import sys

from . import phams
from . import phamgene
from . import phage
from . import utils
from time import sleep
from Bio import SeqIO
import math
import os
from .utils import StarteratorError
import csv
from . import annotate
from collections import Counter

# for the phamerator api pulls:

import json
import os
from .phamerator_api import fetch_phamerator_genes

class Report(object):
    def __init__(self, name=None):
        self.output_dir = utils.INTERMEDIATE_DIR
        self.final_dir = utils.FINAL_DIR
        self.base_name = name

    def make_file(self, specifics, whole_phage=False):
        if whole_phage:
            specifics = ["-a", 'All'] + specifics

        # Get the full path of the current Python interpreter
        python_executable = sys.executable

        args = [python_executable, utils.MAKING_FILES, "-d", utils.INTERMEDIATE_DIR] + specifics
        # sargs = " ".join(args)
        # print sargs
        subprocess.check_call(args)


class PhageReport(Report):
    def __init__(self, name, gui=None, event=None):
        Report.__init__(self, name)
        self.phage_genes = {}
        self.phage_ = []
        self.genes = []
        self._phams = {}
        self.seq_length = None
        self.name = name
        self.gui = gui
        self.event = event

    def final_report(self):
        self.get_phams()

        # ---- Phamerator gap scores (save to JSON for the PDF subprocess)
        try:
            dataset = "Actino_Draft"
            if dataset:
                print("[starterator] pulling Phamerator gaps for", self.name, "dataset", dataset, flush=True)

                genes_json = fetch_phamerator_genes(dataset=dataset, phage=self.name)
                gap_map = gap_map_from_genes(genes_json)

                out_path = os.path.join(self.output_dir, f"{self.name.lower()}_phamerator_gaps.json")
                with open(out_path, "w") as f:
                    json.dump(gap_map, f)

                print("[starterator] wrote", out_path, "entries:", len(gap_map), flush=True)
            else:
                print("[starterator] PHAMERATOR_DATASET not set; skipping Phamerator gaps", flush=True)
        except Exception as e:
            print("[starterator] WARNING: Phamerator gaps failed:", e, flush=True)

        for phm in self._phams.keys():
            if self._phams[phm][0].orientation == 'R' and self._phams[phm][0].start == self.seq_length:
                print('found probable broken gene, deleting ' + self._phams[phm][0].gene_id + ' from list to starterate')
                pass  # change this to delete the entry from self._phams

        self.make_reports()
        self.make_phage_pages()
        final = self.merge_reports()
        return final

    def prepare_report(self):
        self.phage_ = phage.new_phage(name=self.name)
        self.genes = self.phage_.get_genes()

    def get_phams(self):
        self.phage_ = phage.new_phage(name=self.name)
        self._phams = self.phage_.get_phams()
        self.seq_length = self.phage_.length()

    def make_reports(self):
        pham_counter = 0
        total_no = len(self._phams.keys())
        pham_items = list(self._phams.items())
        pham_items.sort(key=lambda item: (item[1][0].start_codon_location, item[0]))
        sorted_keys = [item[0] for item in pham_items]
        for pham_no in sorted_keys:
            if pham_no is not None:
                pham_counter += 1
                gene_name = self._phams[pham_no][0].gene_id
                self.update_gui("Aligning and Graphing %s \nPham %s (%d of %d)"
                                % (gene_name, pham_no, pham_counter, total_no), (float(pham_counter)/total_no))
                gene_report = GeneReport(self.base_name, self.phage_)
                pham_no = gene_report.get_pham(pham_no)
                pham = gene_report.make_report(whole=True)
                # print pham_no
                genes = self._phams[pham_no]
                for gene in genes:
                    # print gene.gene_id
                    phage_gene = pham.genes[gene.gene_id]
                    gene_no = phamgene.get_gene_number(gene.gene_id)
                    self.phage_genes[int(gene_no)] = {'suggested_start': phage_gene.suggested_start["most_called"],
                                                 'gene': pham.genes[gene.gene_id],
                                                 'pham_no': pham_no}

    def make_phage_pages(self):
        pickle_file = os.path.join(self.output_dir, "%s.pickle" % self.name)
        pickle.dump(self.phage_genes, open(pickle_file, "wb"))
        args = ["-p", self.name, "-f", pickle_file, "-l", str(self.seq_length), "-m", "genome"]
        self.make_file(args, True)

    # def make_suggested_starts_page(self):
    #   pickle_file = 
    #   cPickle.dump(self.phage_genes, open(pickle_file, "wb"))

    def merge_reports(self):
        merger = PyPDF2.PdfWriter()
        
        # Add phage genome pages
        with open("%s/%sPhamsGraph.pdf" % (self.output_dir, self.name), 'rb') as phage_genome:
            reader = PyPDF2.PdfReader(phage_genome)
            for page in reader.pages:
                merger.add_page(page)
        
        # Add phage starts pages  
        with open(os.path.join(self.output_dir, "%sSuggestedStarts.pdf" % self.name), 'rb') as phage_starts:
            reader = PyPDF2.PdfReader(phage_starts)
            for page in reader.pages:
                merger.add_page(page)
        phams_added = []
        for gene_no in sorted(self.phage_genes.keys()):
            phage_gene = self.phage_genes[gene_no]
            pham = phage_gene["pham_no"]
            if pham not in phams_added and pham is not None:
                with open(os.path.join(self.output_dir, "%sAllPham%sGraph.pdf" % (self.name, pham)), 'rb') as graph:
                    reader = PyPDF2.PdfReader(graph)
                    for page in reader.pages:
                        merger.add_page(page)
                with open("%s/%sAllPham%sText.pdf" % (self.output_dir, self.name, pham), 'rb') as text:
                    reader = PyPDF2.PdfReader(text)
                    for page in reader.pages:
                        merger.add_page(page)
                phams_added.append(pham)
        merger.write(open(os.path.join(self.final_dir, "%sReport.pdf" % self.name), 'wb'))
        return (os.path.join(self.final_dir,"%sReport.pdf" % self.name), "%sReport.pdf" % self.name)

    def check_stop(self):
        if self.gui:
            if self.event.is_set():
                sys.exit(0)

    def update_gui(self, message, percent):
        if self.gui:
            self.check_stop()
            self.gui.update(message, percent)


class UnPhamPhageReport(PhageReport):
    def __init__(self, name, fasta_file, profile_file, gui=None, event=None):
        PhageReport.__init__(self, name, gui, event)
        self.name = name
        self.fasta = fasta_file
        self.profile = profile_file
        self._phams = None
        self.sequence = None

    def final_report(self):
        self.get_phams()
        self.get_cluster()
        self.make_reports()
        self.make_phage_pages()
        return self.merge_reports()

    def get_sequence(self):
        if not self.sequence:
            try:
                with open(self.fasta, "r") as fasta_file:
                    next(fasta_file)
                    sequence = ""
                    for line in fasta_file:
                        sequence += (line.strip())
                self.sequence = sequence
                self.seq_length = len(sequence)
            except:
                raise StarteratorError("The fasta file (%s) could not be opened!" % self.fasta)
        return self.sequence

    def get_phams(self):
        if not self._phams:
            self._phams = {}
            sequence = self.get_sequence()
            genes = []
            # try:
            if self.profile is None:
                gene_predictions = annotate.auto_annotate(self.fasta)
                for gene in gene_predictions.genes:
                    gene = phamgene.UnPhamGene(gene.id, gene.start, gene.stop, gene.orientation, self.name, sequence)
                    genes.append(gene)
                    pham_no = gene.blast()
                    if pham_no not in self._phams:
                        self._phams[pham_no] = []
                    self._phams[pham_no].append(gene)
            else:
                try:
                    with open(self.profile, "r") as profile:
                        # print self.profile, "has been opened!"
                        first_line = profile.readline()
                        first_word = first_line.split()[0]
                        if first_word == "Profile":
                            csv_reader = csv.reader(profile)
                            line = next(csv_reader)
                            #  print line
                            next(csv_reader)
                            for row in csv_reader:
                                # print row
                                feature_type = row[8].strip()
                                #  print feature_type
                                if feature_type == "ORF":
                                    number = row[1].replace('"', "")
                                    orientation = row[2]
                                    start = int(row[5])
                                    stop = int(row[6])
                                    # print number, start, stop, orientation, self.name
                                    gene = phamgene.UnPhamGene(number, start, stop, orientation, self.name, sequence)
                                    genes.append(gene)

                                    pham_no = gene.phambymatch()
                                    if pham_no is None:
                                        pham_no = gene.blast()

                                    if pham_no not in self._phams:
                                        self._phams[pham_no] = []
                                    self._phams[pham_no].append(gene)

                                    if pham_no is not None:
                                        gene.add_cluster_hits()
                        else:
                            if first_word == "CDS":
                                profile.seek(0)
                                gene_count = 0
                                for line in profile:
                                    if line[0:3] == "CDS":
                                        if line[4:8] == 'join':
                                            continue

                                        else:
                                            gene_count += 1
                                            line2 = line.replace("(", "")
                                            line3 = line2.replace(")", "")
                                            line_items = line3.split()

                                            if line_items[1] == "complement":
                                                gene_orientation = "R"
                                                gene_start = int(line_items[2])
                                                gene_end = int(line_items[4])

                                            else:
                                                gene_orientation = "F"
                                                gene_start = int(line_items[1])
                                                gene_end = int(line_items[3])
                                        gene = phamgene.UnPhamGene(gene_count, gene_start, gene_end, gene_orientation, self.name, sequence)
                                        genes.append(gene)
                                        pham_no = gene.phambymatch()
                                        if pham_no is None:
                                            pham_no = gene.blast()

                                        if pham_no not in self._phams:
                                            self._phams[pham_no] = []
                                        self._phams[pham_no].append(gene)

                                        if pham_no is not None:
                                            gene.add_cluster_hits()

                                    else:
                                        continue
                            else:
                                raise StarteratorError("The profile file (%s) could not be read correctly! Please make sure it is correct." % self.profile)
                except:
                    raise StarteratorError("The profile file (%s) could not be read correctly! Please make sure it is correct." % self.profile)
        return self._phams

    def get_cluster(self):
        cluster_list = []
        subcluster_list = []
        for genelist in self._phams.values():
            for pg in genelist:
                if pg.cluster_hits is not None:
                    for i in pg.cluster_hits:
                        if i is not None:
                            cluster_list.append(i)
                if pg.subcluster_hits is not None:
                    for i in pg.subcluster_hits:
                        if i is not None:
                            subcluster_list.append(i)

        cluster_counts = Counter(cluster_list)
        subcluster_counts = Counter(subcluster_list)
        cluster_guess, guess_count = cluster_counts.most_common(1)[0]
        penult_cluster, penult_count = cluster_counts.most_common(2)[1]
        total_genes = sum([len(gl) for gl in self._phams.values()])
        best_fraction = guess_count/(total_genes * 1.0)
        penult_fraction = penult_count/(total_genes * 1.0)
        if best_fraction > 0.9 and (best_fraction - penult_fraction) > 0.2:
            cluster = cluster_guess
        else:
            cluster = "Unassigned"

        # start of changes
        if cluster == "Unassigned":
            subcluster = "Unassigned"
            cluster = "Unassigned"
        else:
            # If no subcluster hits, cluster = subcluster
            if not subcluster_counts:
                subcluster = cluster
            else:
                subcluster_guess, subguess_count = subcluster_counts.most_common(1)[0]
                # no more error with only one item
                top2 = subcluster_counts.most_common(2)
                penult_subcount = top2[1][1] if len(top2) > 1 else 0

                best_subfraction = subguess_count / (total_genes * 1.0)
                penult_subfraction = penult_subcount / (total_genes * 1.0)
                if best_subfraction > 0.97 or (
                        best_subfraction > 0.9 and (best_subfraction - penult_subfraction) > 0.10):
                    subcluster = subcluster_guess
                    # dont overwrite subcluster
                elif cluster != "Unassigned":
                    subcluster = cluster
                else:
                    subcluster = "Unassigned"
                    cluster = "Unassigned"
        # end of changes

        for genelist in self._phams.values():
            ch = sum([pow(ord(elem), i + 1) for i, elem in enumerate(subcluster)])
            for pg in genelist:
                pg.cluster = cluster
                pg.subcluster = subcluster
                pg.cluster_hash = ch


    def make_reports(self):
        self.phage_genes = {}
        pham_counter = 0
        total_no = len(self._phams)
        for pham_no in self._phams:
            genes = self._phams[pham_no]
            if pham_no is not None:
                self.update_gui("Aligning and Graphing Pham %s \n Pham %d of %d"
                                % (pham_no, pham_counter, total_no), float(pham_counter)/total_no)
                gene_report = GeneReport(self.base_name, self.name, whole_phage=True)
                pham_no = gene_report.get_pham(pham_no, genes)
                pham = gene_report.make_report(True)
                for gene in genes:
                    phage_gene = pham.genes[gene.gene_id]
                    gene_no = phamgene.get_gene_number(gene.gene_id)
                    self.phage_genes[gene_no] = {'suggested_start': phage_gene.suggested_start["most_called"],
                                                 'gene': gene,
                                                 'pham_no': pham_no}
            else:
                for gene in genes:
                    self.phage_genes[int(gene_no)] = {'pham_no': None, "gene": gene, "suggested_start": None}
            pham_counter += 1


class GeneReport(Report):
    def __init__(self, phage_name, number=None, whole_phage=False, fasta_file=None):
        """

        :param phage_name: name of phage
        :param number: gene number
        :param whole_phage: bool
        :param fasta_file: name of file
        """
        Report.__init__(self)
        self.phage_name = phage_name
        self.number = number
        self.all = whole_phage
        self.fasta = fasta_file
        self.pham = None
        self.sequence = None
        self.seq_length = None

    def get_pham(self, pham_no=None, genes=None):
        if pham_no and genes:
            self.pham = phams.Pham(pham_no, genes)
        elif pham_no:
            self.pham = phams.Pham(pham_no)
        else:
            # print self.phage_name, self.number
            pham_no = phamgene.get_pham_no(self.phage_name, self.number)
            self.pham = phams.Pham(pham_no)
        return pham_no

    def get_sequence(self):
        with open(self.fasta) as fasta_file:
            next(fasta_file)
            sequence = ""
            for line in fasta_file:
                sequence += (line.strip())
        self.sequence = sequence
        self.seq_length = len(sequence)
        return self.sequence

    def make_unpham_gene(self, start, stop, orientation):
        self.get_sequence()
        gene = phamgene.UnPhamGene(self.number, start, stop, orientation, self.phage_name, self.sequence)
        self.pham_no = gene.blast()
        self.get_pham(self.pham_no, [gene])

    def make_report(self, whole=False):
        self.pham.align()
        self.pham.find_most_common_start()
        pickle_file = os.path.join(self.output_dir, "%s%s.pickle" % (self.pham.file, self.pham.pham_no)) # TODO: Figure out base name things
        f = open(pickle_file, "wb")
        pickle.dump(self.pham, f)
        f.close()

        args = ["-p", self.phage_name, "-n", str(self.pham.pham_no),  "-f", pickle_file, "-m", "text"]
        self.make_file(args, whole)
        # graph_args = ["-p", self.phage, "-n", pham_no, "-f", pickle_file, "-m", "graph"]
        return self.pham

    def get_specific_gene(self):
        # gene = self.pham.genes[]
        pass

    def merge_report(self):
        """
        """
        merger = PyPDF2.PdfWriter()
        with open(os.path.join(self.output_dir, "%sOnePham%sGraph.pdf" % (self.phage_name, self.pham.pham_no)), "rb") as graph:
            reader = PyPDF2.PdfReader(graph)
            for page in reader.pages:
                merger.add_page(page)
        with open(os.path.join(self.output_dir, '%sOnePham%sText.pdf' % (self.phage_name, self.pham.pham_no)), 'rb') as text:
            reader = PyPDF2.PdfReader(text)
            for page in reader.pages:
                merger.add_page(page)
        file_path = os.path.join(self.final_dir, "%sPham%sReport.pdf" % (self.phage_name, self.pham.pham_no))
        merger.write(open(file_path, 'wb'))
        return file_path, "%sPham%sReport.pdf" % (self.phage_name, self.pham.pham_no)


class UnPhamGeneReport(GeneReport):
    def __init__(self, base_name, phage_name, number, start, stop, orientation, fasta_file):
        GeneReport.__init__(self, base_name, phage_name, number)
        self.start = start
        self.stop = stop
        self.orientation = orientation
        self.fasta_file = fasta_file

    def get_sequence(self):
        if not self.sequence:
            try:
                with open(self.fasta_file) as fasta_file:
                    next(fasta_file)
                    sequence = ""
                    for line in fasta_file:
                        sequence.join(line.strip())
                self.sequence = sequence
                self.seq_length = len(sequence)
            except:
                raise StarteratorError("The fasta file (%s) could not be opened!" % self.fasta_file)
        return self.sequence

    def make_gene(self, start, stop, orientation):
        sequence = self.get_sequence()
        try:
            gene = phamgene.UnPhamGene(self.number, start, stop, orientation, self.name, sequence)
        except:
            raise StarteratorError("The gene could not be made! Check coordinates: Start: %s, Stop: %s, Orientation: %s"
                                   % (start, stop, orientation))
        self.pham_no = gene.blast()
        self.get_pham(self.pham_no, gene)

    def get_pham(self, pham_no, gene):
        self.pham = phams.Pham(pham_no, [gene])
        pham.add(gene)
        return pham_no


class PhamReport(Report):
    def __init__(self, pham_no):
        Report.__init__(self)
        self.pham_no = pham_no

    def final_report(self, save_json=False):
        if save_json:
            self.make_report(save_json=True)
        else:
            self.make_report()
        return self.merge_report()

    def make_report(self, save_json=False):
        self.pham = phams.Pham(self.pham_no)
        self.pham.align()
        self.pham.find_most_common_start()
        pickle_file = os.path.join(self.output_dir, "%s.pickle" % self.pham.pham_no)  # TODO:Figure out base name things
        f = open(pickle_file, "wb")
        pickle.dump(self.pham, f)
        f.close()
        if save_json:
            json_file = pickle_file.replace(".pickle", ".json")
            self.pham.export_json(json_file)

        args = ["-n", self.pham_no, "-f", pickle_file, '-m', "text"]
        self.make_file(args)

    def merge_report(self):
        merger = PyPDF2.PdfWriter()
        with open(os.path.join(self.output_dir, "OnePham%sGraph.pdf" % self.pham_no), "rb") as graph:
            reader = PyPDF2.PdfReader(graph)
            for page in reader.pages:
                merger.add_page(page)
        with open(os.path.join(self.output_dir, 'Pham%sText.pdf' % self.pham_no), 'rb') as text:
            reader = PyPDF2.PdfReader(text)
            for page in reader.pages:
                merger.add_page(page)
        file_path = os.path.join(self.final_dir, "Pham%sReport.pdf" % self.pham_no)
        merger.write(open(file_path, 'wb'))
        return file_path, "Pham%sReport.pdf" % self.pham_no