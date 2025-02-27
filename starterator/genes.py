from .utils import *
from .phams import *
from reportlab.lib.pagesizes import letter, A4
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
import math
import PyPDF2
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastpCommandline as Blastp
from Bio.Blast.Applications import BlastallCommandline
import MySQLdb

class GeneReport(object):
    def __init__(self, number, phage, is_phamerated=True, fasta_file=None, is_all=False):
        self.phage = phage
        self.number = number
        self.name = phage + '_' + str(number)
        self.is_phamerated = is_phamerated
        if self.is_phamerated:
            self.pham_no = get_pham_no(self.db, self.phage, self.number)
        self.is_whole_phage = is_all
        self.fasta = fasta_file

    def blast_gene(self):
        """
            Function that performs BLAST querying the unphamed gene's protein sequence on
            the Mycobacteriophage Proteins database
            Uses an E-value of 10E-30 as in paper for Phamily formation
            Finds the first (best hit) result of the query and returns the line
        """
        # Create a fasta file containing the gene
        SeqIO.write(self.protein, '%s%s.fasta' % (self.output_dir, self.name), 'fasta')
        # Run blast program with the following parameters
        e_value = math.pow(10, -30)
        print(self.legacy_blast) 
        if not self.legacy_blast:
            blast_command = Blastp(
                            query='%s%s.fasta' % (self.output_dir, self.name),
                            db="\"%s\"" % self.protein_db, evalue=e_value, outfmt=5, 
                            out='%s%s.xml' % (self.output_dir, self.name))

            # 'blastp 
            # -out "/home/marissa/.starterator/Intermediate Files/Nilo_1.xml" 
            # -outfmt 5 
            # -query "/home/marissa/.starterator/Intermediate Files/Nilo_1.fasta"
            #  -db /home/marissa/.starterator/Protein Database/Proteins 
            #  -evalue 1e-30
            # blast_command = [self.blast_dir + 'blastp', 
            #     '-out', '\"%s%s.xml\"' % (self.output_dir, self.name),
            #     "-outfmt", '5',
            #     "-query", "\"%s%s.fasta\"" % (self.output_dir, self.name),
            #      '-db', "\"%s\"" % self.protein_db, 
            #      '-evalue', '1e-30' ]
            # for f in blast_command: print f,
            # subprocess.check_call(blast_command)
            stdout, stderr = blast_command()
        else:
            # blastall -p blastn -d nr -i QUERY -o out.QUERY
            # 'blastall -d /home/pbi/Dropbox/StarteratorUI/proteindb/ProteinsDB 
            # -i /home/pbi/Dropbox/StarteratorUI/intermediate_files/Nilo_1.fasta -e 1e-50 
            # -o /home/pbi/Dropbox/StarteratorUI/intermediate_files/Nilo_1.xml
            # blast_command = BlastallCommandline(program=self.blast_dir+ 'blastall',
            #                 database=self.protein_db, infile='%s%s.fasta' % (self.output_dir, self.name),
            #                 outfile='%s%s.xml' % (self.output_dir, self.name), expectation=e_value)
            blast_command = [self.blast_dir + 'blastall', '-p', 'blastp', '-d', "%s" % self.protein_db, '-e', '1e-50',
                '-i', '%s%s.fasta' % (self.output_dir, self.name),
                '-o', '%s%s.xml' % (self.output_dir, self.name), '-m', '7' ]
            # stout, stderr = blast_command()
            subprocess.check_call(blast_command)

        try:
            result_handle = open('%s%s.xml' % (self.output_dir, self.name))
            blast_record = NCBIXML.read(result_handle)
        except:
            result_handle.close()
            result_handle = open('%s%s.xml' % (self.output_dir, self.name))
            blast_records = NCBIXML.parse(result_handle)
            blast_record = next(blast_records)
        # If there is at least one blast result, return the first result, otherwise returns None
        # E-value of 10 exp -30 or better are grouped into the same family
        if len(blast_record.descriptions) > 0:
            first_result = blast_record.descriptions[0].title
            temp = first_result.split(',')
            first_result_name = temp[0].split(' ')[1]
            phage_name = first_result_name.split('_')[0]
            gene_number = first_result_name.split('_')[-1]
            print(phage_name, gene_number)
            self.pham_no =get_pham_no(self.db, phage_name, gene_number)
        else:
            self.pham_no = None
        return self.pham_no

    def make_unphamerated_gene(self, start, stop, orientation):
        print('fasta', self.fasta)
        seq_file = open(self.fasta, 'r')
        first_line = next(seq_file) 
        sequence = ""
        for line in seq_file:
            sequence += line.strip()
        seq_file.close()

        self.gene_record = make_gene(self.name, start-1, stop, orientation, sequence, self.phage)
        self.protein = SeqRecord(self.gene_record.annotations['translation'], id=self.name)

    def load_unphamerated_gene(self, gene_record, protein):
        """
        """
        self.gene_record = gene_record
        self.protein = protein

    def make_report(self, genes=None):
        one_or_all = 'All' if self.is_whole_phage else 'One'
        file_name = None if self.is_phamerated else "%sPham%s" % (self.phage, self.pham_no)
        pham = Pham(self.pham_no, self.phage, root_file_name=file_name)
        if  not self.is_phamerated and not self.is_whole_phage:
            pham.add_gene(self.gene_record)
        if genes != None:
            for gene in genes:
                pham.add_gene(gene.gene_record)
        pham.align(self.output_dir, one_or_all, already_aligned=False)
        pham.put_similar_genes_together()
        start_stats = pham.most_common_start()
        pickle_file = utils.INTERMEDIATE_DIR + '%sPham%s.p' %((self.phage + one_or_all), self.pham_no)
        cPickle.dump(pham, open(pickle_file, 'wb'))

        subprocess.check_call(['python', utils.MAKING_FILES,
            '-p', self.phage, 
            '-n', self.pham_no,
            '-a', one_or_all,
            '-f', "\"%s\"" % pickle_file,
            '-o', "\"%s\"" % utils.INTERMEDIATE_DIR,
            '-d', "\"%s\"" % utils.FINAL_DIR,
            '-m', 'text'])
        subprocess.check_call(['python', utils.MAKING_FILES,
            '-p', self.phage, 
            '-n', self.pham_no,
            '-a', one_or_all,
            '-f', "\"%s\"" % pickle_file,
            '-o', "\"%s\"" % utils.INTERMEDIATE_DIR,
            '-d', "\"%s\"" % utils.FINAL_DIR,
            '-m', 'graph'])
        return pham


    def one_pham_report(self):
        """
            Merges the graphical and text output of one pham into one PDF file
            {Phage}Pham{Number}Report.pdf
        """
        merger = PyPDF2.PdfFileMerger()
        graph = open("%sPham%sGraph.pdf" % (utils.INTERMEDIATE_DIR + self.phage + 'One', self.pham_no), "rb")
        text = open('%s%sPham%sText.pdf' % (utils.INTERMEDIATE_DIR, self.phage +'One', self.pham_no), 'rb')
        merger.append(fileobj=graph)
        merger.append(fileobj=text)
        file_path = "%s%sPham%sReport.pdf" % (utils.FINAL_DIR, self.phage, self.pham_no)
        merger.write(open(file_path, 'wb'))
        return file_path, "%sPham%sReport.pdf" % (self.phage, self.pham_no)