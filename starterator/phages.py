from . import genes
from . import utils
from .genes import *
from .phams import *
from .utils import *
from reportlab.lib.pagesizes import letter, A4
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
import PyPDF2
from reportlab.lib import colors
from .phage import new_phage

class PhageReport(object):
    def __init__(self, name, is_phamerated, fasta_file=None, profile_file=None, gui=None, event=None):
        self.name = name
        self.phage = new_phage(name=name)
        self.genes = {}
        self.is_phamerated = is_phamerated
        self.intermediate_dir = utils.INTERMEDIATE_DIR
        self.final_dir = utils.FINAL_DIR
        self.fasta = fasta_file
        self.profile = profile_file
        self.unphamerated_genes = {}
        self.config = config
        self.gui = gui
        self.current_pham = [0, 0]
        self.total_phams = None
        self.event = event

    def make_report(self):
        if self.is_phamerated:
            print('before phams found')
            phams = self.phage.get_phams()
            seq_length = self.phage.length()
            print('after phams found')
            for gene_name, pham_no in phams:
                if self.gui:
                    if self.event.is_set():
                        print('stopped in make_report')
                        sys.exit(0)
                gene_number = utils.get_gene_number(gene_name)
                # gene_id = self.name +'_' + str(gene_number)

                gene = genes.Gene(gene_number, self.name, self.is_phamerated, self.config, self.db, is_all=True)
                self.genes[gene_number] = (gene, pham_no)
            return self.phage_report(seq_length)
        else:
            
            if self.gui:
                print('in unphamerated pham')
                if self.event.is_set():
                    sys.exit(0)
                self.gui.update('Finished Making Unphamerated Genes.', .02)
            for gene_number in self.genes:
                if self.gui:
                    if self.event.is_set():
                        sys.exit(0)
                gene = self.genes[gene_number][0]
                pham_no = self.genes[gene_number][1]
                if pham_no == None:
                    pass
                else:
                    if pham_no not in self.unphamerated_genes:
                        self.unphamerated_genes[pham_no] = []
                    self.unphamerated_genes[pham_no].append(gene)
            return self.phage_report(seq_length)

    def phage_report(self, seq_length):
        phams_done = []
        pham_genes = {}
        self.total_phams = len(set(self.genes.values()))
        for gene_number in self.genes:
            if self.gui and self.event.is_set():
                sys.exit(0)
            self.current_pham[0] += 1
            gene = self.genes[gene_number][0]
            pham_no = self.genes[gene_number][1]
            self.current_pham[1] = pham_no
            if pham_no not in phams_done and pham_no is not None:
                if not self.is_phamerated:
                    unphamerated_genes = self.unphamerated_genes[pham_no]
                else:
                    unphamerated_genes = None
                if self.gui:
                    if self.event.is_set():
                        sys.exit(0)
                    self.gui.update('Aligning and Graphing Pham %s (Gene %s) \n' % (pham_no, gene_number) + 
                        'Pham %d of %d' % (self.current_pham[0], self.total_phams)
                        , (float(self.current_pham[0])/ self.total_phams) - .05 )
                pham = gene.make_report(genes=unphamerated_genes)

                for gene_record in pham.phage_gene:
                    gene_num = utils.get_gene_number(gene_record.id)
                    suggested_start = pham.suggest_start_site_coord(gene_record.id)
                    gene_record.annotations['suggested_start'] = suggested_start
                    pham_genes[gene_num] = gene_record, gene.pham_no
                phams_done.append(pham_no)
        if self.gui:
            if self.event.is_set():
                sys.exit(0)
            self.gui.update('Making Final Report', .97)
        self.make_pham_genome(pham_genes, seq_length)
        self.make_suggested_starts_page(pham_genes)
        return self.make_phage_report(pham_genes, seq_length)
    
    def phamerated_phage_report(self):
        phage_phams, seq_length = find_phams_of_a_phage(self.db, self.name)
        for gene_name, pham_no in phage_phams:
            gene_number = utils.get_gene_number(gene_name)
            gene_id = self.name +'_' + str(gene_number)
            gene = GeneReport(gene_number, self.name, self.is_phamerated, is_all=True)
            self.genes[gene_number] = (gene, pham_no)
        return self.phage_report(seq_length)

    
    def unphamerated_phage_report(self):
        seq_length = self.load_unphamerated_genes()
        if self.gui:
            if self.event.is_set():
                sys.exit()
            self.gui.update('Finished Making Unphamerated Genes.', .02)
        for gene_number in self.genes:
            gene = self.genes[gene_number][0]
            pham_no = self.genes[gene_number][1]
            if pham_no == None:
                pass
            else:
                if pham_no not in self.unphamerated_genes:
                    self.unphamerated_genes[pham_no] = []
                self.unphamerated_genes[pham_no].append(gene)
        return self.phage_report(seq_length)

    def load_unphamerated_genes(self):
        """
        Function that takes in the fasta file of a unPhamerated phage and
        it's Annotation Profile comma separated (CSV) file from DNAMaster : 
        Using information from these files, it returns a list of UnphamedGenes of the phage
        """
        try:
            if self.gui:
                if self.event.is_set():
                    sys.exit(0)
                self.gui.update('Making Unphamerated Genes', .01)
            seq_file = open(self.fasta, 'r')
            first_line = next(seq_file) 
            seq_info = first_line.split(',')[1]
            seq_length = int(seq_info.strip().split(' ')[0])
            sequence = ""
            for line in seq_file:
                sequence += line.strip()
            seq_file.close()

            genes_file = open(self.profile, 'rb')
            next(genes_file)
            next(genes_file)
            for line in genes_file:
                temp = line.split(',')
                feature_type = temp[7].strip()
                if feature_type == 'ORF':
                    #  Only want possible protein encoding genes
                    number = temp[1].replace('"', '')
                    gene_id = self.name +"_" + str(number)
                    orientation = temp[2]
                    start = int(temp[4]) -1
                    stop = int(temp[5])
                    # Orientation: if F, start and stop are as said
                    # If reverse, start is stop, and stop is start
                    pham_gene = make_gene(gene_id, start, stop, orientation, sequence, self.name)
                    ahead_of_start = pham_gene.annotations['number_before_start']
                    protein = SeqRecord(pham_gene.seq[ahead_of_start:].translate(), id=gene_id)
                    gene = Gene(number, self.name, self.is_phamerated, self.config, self.db, is_all=True)
                    gene.load_unphamerated_gene(pham_gene, protein)
                    pham_no = gene.blast_gene()
                    self.genes[number] = (gene, pham_no)
            genes_file.close()
            return seq_length
        except:
            raise #Exception("Could not read the DNAMaster Profile")



    def make_pham_genome(self, pham_genes, seq_length):
        """
            Creates PDF Graph of the phage genome with the pham number and color
            labeling each gene.
            {Phage}PhamsGraph.pdf
        """
        pham_colors = get_pham_colors(self.db)
        gd_diagram = GenomeDiagram.Diagram(self.name)
        gd_track = gd_diagram.new_track(1, name=self.name, greytrack=1)
        gd_pham_set = gd_track.new_set()
        print("making genome page")
        for gene_num in pham_genes:
            pham = pham_genes[gene_num][1]
            gene = pham_genes[gene_num][0]
            print(pham, gene.id)
            if pham == None:
                pham_no = "None"
                pham_color = 'Black'
            else:
                pham_no = pham
                pham_color = colors.HexColor(pham_colors[pham_no])
            if gene.annotations['orientation'] == 'F':
                gene_location = FeatureLocation(gene.annotations['start'], gene.annotations['stop'])
                strand = 1
            else:
                strand = -1
                gene_location = FeatureLocation(gene.annotations['stop'], gene.annotations['start'])
            gene_feature = SeqFeature(gene_location, strand=strand)
            gene_number = get_gene_number(gene.id)
            # label the gene with the gene number
            gd_pham_set.add_feature(gene_feature, name=str(gene_number), label=True, 
                label_size=6, label_angle=75)
            # label gene with pham color and name
            gd_pham_set.add_feature(gene_feature, color=pham_color, name=pham_no, label=True, label_position='middle')
        
        gd_diagram.draw(format='linear', orientation='portrait', pagesize=letter, fragments=8, start=0, end=seq_length)
        gd_diagram.write('%sPhamsGraph.pdf' % (self.intermediate_dir + self.name), "PDF")
        # gd_diagram.write('%sPhamsGraph.png' % (output_dir + phage), "PNG")

    def make_suggested_starts_page(self, pham_genes):
        """
            Creates a PDF page of the suggested starts of a phage
            Genes are list in order
            {Gene Name} is a member of Pham {Number}: {Suggested Start Coordinates}
        """
        doc = SimpleDocTemplate("%sSuggestedStarts.pdf" % (self.intermediate_dir + self.name), pagesize=letter)
        story = []
        print("making suggested starts page")
        styles = getSampleStyleSheet()
        styles.add(ParagraphStyle(name="paragraph"))
        styles.add(ParagraphStyle(name='Center', alignment=TA_CENTER))
        text = '<font size=14> Suggested Start Coordinates</font>'
        story.append(Paragraph(text, styles['Center']))
        story.append(Spacer(1, 12))
        for gene_no in sorted(pham_genes.keys()):
            gene = pham_genes[gene_no][0]
            pham = pham_genes[gene_no][1]
            if pham == None:
                text = '<font size=12> %s is not a member of an existing Pham </font>' % (gene.id)
            else:
                suggested_start = gene.annotations['suggested_start']
                text = '<font size=12> %s is a member of Pham %s:  %s </font>' % (gene.id, pham, suggested_start)
            story.append(Paragraph(text, styles['Normal']))
        doc.build(story)

    def make_phage_report(self, pham_genes, seq_length):
        """
            Creates a single PDF file for a phage, combines the pham genome graph,
            the suggested starts page, and graph and text output for each pham in the phage
            {Phage}Report.pdf
        """
        one_or_all = 'All'
        print("making report")
        merger = PyPDF2.PdfFileMerger()
        phage_starts = open("%sSuggestedStarts.pdf" % (self.intermediate_dir + self.name), 'rb')
        phage_genome = open('%sPhamsGraph.pdf' % (self.intermediate_dir + self.name), 'rb')
        merger.append(fileobj=phage_genome)
        merger.append(fileobj=phage_starts)
        phams_added = []
        for gene_num in sorted(pham_genes.keys()):
            pham = pham_genes[gene_num][1]
            if pham not in phams_added and pham != None:
                graph = open("%sPham%sGraph.pdf"  % (self.intermediate_dir + self.name + one_or_all, pham), "rb")
                text = open('%sPham%sText.pdf' % (self.intermediate_dir + self.name + one_or_all, pham), 'rb')
                merger.append(fileobj=graph)
                merger.append(fileobj=text)
                phams_added.append(pham)
        merger.write(open("%s%sReport.pdf" % (self.final_dir, self.name), 'wb'))
        return "%s%sReport.pdf" % (self.final_dir, self.name), '%sReport.pdf' % self.name