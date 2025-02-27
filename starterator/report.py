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

class Report(object):
    def __init__(self, name=None):
        self.output_dir = utils.INTERMEDIATE_DIR
        self.final_dir = utils.FINAL_DIR
        self.base_name = name

    def make_file(self, specifics, whole_phage=False):
        if whole_phage:
            specifics = ["-a", 'All'] + specifics
        args = ['python', utils.MAKING_FILES, "-d", utils.INTERMEDIATE_DIR] + specifics
        sargs = " ".join(args)
        print(sargs)
        subprocess.check_call(args)



class PhageReport(Report):
    def __init__(self, name, gui=None, event=None):
        Report.__init__(self, name)
        self.name = name
        self.gui = gui
        self.event = event
    
    def final_report(self):
        self.get_phams()
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
        self.phage_genes = {}
        pham_counter = 0
        total_no = len(list(self._phams.keys()))
        for pham_no in list(self._phams.keys()):
            if pham_no is not None:
                pham_counter += 1
                self.update_gui("Aligning and Graphing Pham %s \n Pham %d of %d"
                    % (pham_no, pham_counter, total_no), (float(pham_counter)/total_no))
                gene_report = GeneReport(self.base_name, self.phage_)
                pham_no = gene_report.get_pham(pham_no)
                pham = gene_report.make_report(whole=True)
                print(pham_no)
                genes = self._phams[pham_no]
                for gene in genes:
                    print(gene.gene_id)
                    phage_gene = pham.genes[gene.gene_id]
                    gene_no = phamgene.get_gene_number(gene.gene_id)
                    self.phage_genes[gene_no] = {'suggested_start': phage_gene.suggested_start["most_called"],
                                            'gene': gene,
                                            'pham_no': pham_no}

    
    def make_phage_pages(self):
        pickle_file = os.path.join(self.output_dir, "%s.pickle" % (self.name))
        pickle.dump(self.phage_genes, open(pickle_file, "wb"))
        args = ["-p", self.name, "-f", pickle_file, "-l", str(self.seq_length), "-m", "genome"]
        self.make_file(args, True)

    # def make_suggested_starts_page(self):
    #   pickle_file = 
    #   cPickle.dump(self.phage_genes, open(pickle_file, "wb"))

    def merge_reports(self):
        merger = PyPDF2.PdfFileMerger()
        phage_starts = open(os.path.join(self.output_dir, "%sSuggestedStarts.pdf" % (self.name)))
        phage_genome = open("%s/%sPhamsGraph.pdf" %(self.output_dir, self.name))
        merger.append(fileobj=phage_genome)
        merger.append(fileobj=phage_starts)
        phams_added = []
        for gene_no in sorted(self.phage_genes.keys()):
            phage_gene = self.phage_genes[gene_no]
            pham = phage_gene["pham_no"]
            if pham not in phams_added and pham != None:
                graph = open(os.path.join(self.output_dir,"%sAllPham%sGraph.pdf" %(self.name, pham)))
                text = open("%s/%sAllPham%sText.pdf" %(self.output_dir, self.name, pham))              
                merger.append(fileobj=graph)
                merger.append(fileobj=text)
                phams_added.append(pham)
        merger.write(open(os.path.join(self.final_dir, "%sReport.pdf" %(self.name)), 'wb'))
        return (os.path.join(self.final_dir,"%sReport.pdf" % (self.name)), "%sReport.pdf" % self.name)

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
        self.make_reports()
        self.make_phage_pages()
        return self.merge_reports()

    def get_sequence(self):
        if not self.sequence:
            try:
                with open(self.fasta, "rb") as fasta_file:
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
            if self.profile == None:
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
                    with open(self.profile, "rbU") as profile:
                        print(self.profile, "has been opened!")
                        csv_reader = csv.reader(profile)
                        line = next(csv_reader)
                        print(line)
                        next(csv_reader)
                        for row in csv_reader:
                            print(line)
                            feature_type = row[7].strip()
                            print(feature_type)
                            if feature_type == "ORF":
                                number = row[1].replace('"', "")
                                orientation = row[2]
                                start = int(row[4])
                                stop = int(row[5])
                                print(number, start, stop, orientation, self.name)
                                gene = phamgene.UnPhamGene(number, start, stop, orientation, self.name, 
                                    sequence)
                                genes.append(gene)

                                pham_no = gene.blast()
                                if pham_no not in self._phams:
                                    self._phams[pham_no] = []
                                self._phams[pham_no].append(gene)
                except:
                    raise StarteratorError("The profile file (%s) could not be read correctly! Please make sure it is correct." % self.profile)
        return self._phams

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
                    self.phage_genes[gene.gene_id] = {'pham_no': None,
                     "gene": gene,
                     "suggested_start": None}
            pham_counter += 1



class GeneReport(Report):
    def __init__(self, phage_name, number=None, whole_phage=False, fasta_file=None):
        Report.__init__(self)
        self.phage_name = phage_name
        self.number = number
        self.all = whole_phage
        self.fasta = fasta_file

    def get_pham(self, pham_no=None, genes=None):
        if pham_no and genes:
            self.pham = phams.Pham(pham_no, genes)
        elif pham_no:
            self.pham = phams.Pham(pham_no)
        else:
            print(self.phage_name, self.number)
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
        gene = phamgene.UnPhamGene(self.number, start, stop, orientation, self.phage_name, 
                           self.sequence)
        self.pham_no = gene.blast()
        self.get_pham(self.pham_no, [gene])

    def make_report(self, whole=False):
        self.pham.align()
        self.pham.find_most_common_start()
        pickle_file = os.path.join(self.output_dir, "%s%s.pickle" % (self.pham.file, self.pham.pham_no))#TODO: Figure out base name things
        print(pickle_file, self.output_dir)
        f = open(pickle_file, "wb")
        pickle.dump(self.pham, f)
        print(pickle_file)
        sleep(2)
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
        merger = PyPDF2.PdfFileMerger()
        graph = open(os.path.join(self.output_dir, "%sOnePham%sGraph.pdf" % (self.phage_name, self.pham.pham_no)), "rb")
        text = open(os.path.join(self.output_dir, '%sOnePham%sText.pdf' % (self.phage_name, self.pham.pham_no)), 'rb')
        merger.append(fileobj=graph)
        merger.append(fileobj=text)
        file_path = os.path.join(self.final_dir, "%sPham%sReport.pdf" % (self.phage_name, self.pham.pham_no))
        merger.write(open(file_path, 'wb'))
        return file_path, "%sPham%sReport.pdf" % (self.phage_name, self.pham.pham_no)
        

class UnPhamGeneReport(GeneReport):
    def __init__(self, base_name, phage_name, number, start, stop, orientation, fasta_file):
        GeneReport.__init__(self, base_name, phage_name, number)
        self.start = start
        self.stop = stop
        self.orientation = orientation
        self.fasta = fasta

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
            gene = phamgene.UnPhamGene(number, start, stop, orientation, self.name, 
                                sequence)
        except:
            raise StarteratorError("The gene could not be made! Check to make sure coordinates are correct. Start: %s, Stop: %s, Orientation: %s" % 
                (start, stop, orientation))
        self.pham_no = gene.blast()
        self.get_pham(self.pham_no, gene)


    def get_pham(self, pham_no, gene):
        self.pham = pham.Pham(pham_no, [gene])
        pham.add(gene)
        return pham_no

class PhamReport(Report):
    def __init__(self, pham_no):
        Report.__init__(self)
        self.pham_no = pham_no
    
    def final_report(self):
        self.make_report()
        return self.merge_report()

    def make_report(self):
        self.pham = phams.Pham(self.pham_no)
        self.pham.align()
        self.pham.find_most_common_start()
        pickle_file = os.path.join(self.output_dir, "%s.pickle" % (self.pham.pham_no)) #TODO: Figure out base name things
        f = open(pickle_file, "wb")
        pickle.dump(self.pham, f)
        f.close()
        args = ["-n", self.pham_no, "-f", pickle_file, '-m', "text"]
        self.make_file(args)

    def merge_report(self):
        merger = PyPDF2.PdfFileMerger()
        graph = open(os.path.join(self.output_dir, "OnePham%sGraph.pdf" % ( self.pham_no)), "rb")
        text = open(os.path.join(self.output_dir,'Pham%sText.pdf' % (self.pham_no)), 'rb')
        merger.append(fileobj=graph)
        merger.append(fileobj=text)
        file_path = os.path.join(self.final_dir, "Pham%sReport.pdf" % (self.pham_no))
        merger.write(open(file_path, 'wb'))
        return file_path, "Pham%sReport.pdf" % ( self.pham_no)
        
