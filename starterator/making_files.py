#!/usr/bin/env python

# Copyright (c) 2013, 2014 All Right Reserved, Hatfull Lab, University of Pittsburgh
#
# THIS CODE AND INFORMATION ARE PROVIDED "AS IS" WITHOUT WARRANTY OF ANY
# KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A
# PARTICULAR PURPOSE.  USE AT YOUR OWN RISK.
#
# Marissa Pacey
# April 4, 2014
# Functions that create the PDF outputsst

import pickle
import argparse
from Bio.Graphics import GenomeDiagram
from Bio.SeqFeature import SeqFeature, FeatureLocation
from reportlab.lib import colors
from reportlab.lib.pagesizes import letter, A4
from reportlab.pdfgen import canvas
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib import colors
import PyPDF2
from Bio import SeqIO
import math
import io
from . import utils
from . import phams
from . import phamgene
import os
# from phage import



def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--pham_no', default = -1)
    parser.add_argument('-p' , '--phage', default=None)
    parser.add_argument('-d', '--dir', default=None)
    parser.add_argument('-a', '--one_or_all', default='One')
    parser.add_argument('-f', '--pickle_file', default=None)
    parser.add_argument("-l", "--phage_length", type=int, default=-1)
    parser.add_argument('-m', '--make', default=None)

    return parser.parse_args()

def check_file(file_name):
    try: 
        f = open(file_name, "r")
        return True
    except:
        return False

def output_suggested_starts(pham, phage=None, all_genes=False):
    """
        gives a String list of the output of the suggested start site
        If all_genes is true, finds the suggested start site of all genes in a pham
        Otherwise, finds only the suggested start site of the specific phage gene(s)
    """
    output = []
    if all_genes == False:
        for gene in list(pham.genes.values()):
            if phage in gene.gene_id:
                coord_suggestion = gene.suggested_start["most_called"]
                output.append(""+ gene.gene_id + ", " +str(coord_suggestion))
    else:
        for gene in list(pham.genes.values()):
            coord_suggestion = gene.suggested_start["most_called"]
            output.append(gene.gene_id + ", " +str(coord_suggestion))
    return output


def output_start_sites(stats):
        """
            Writes a report of the start sites statistics
            Returns a list of strings with information about
            start sites for this pham  
        """
        most_called_start = stats["most_called_start"]
        total_genes = len(stats["most_called"])+ len(stats["most_not_called"]) + len(stats["no_most_called"])
        output = []
        output.append("Most Called Start: %s (number based on diagram)"
                         % most_called_start)
        percent_with_most_called = (float(len(stats["most_called"])) 
                                    /total_genes *100 )
        output.append("Percent with start called: %10.4f%%" 
                        % percent_with_most_called )
        output.append("Genes that call the most called start:")
        s = '\u2022' + ''
        for gene in stats["most_called"]:
            s += gene+ ", "
        output.append(s)
        output.append("")
        output.append("Genes that have the most called start but do not call it:")
        s = '\u2022' + ''
        for gene in stats["most_not_called"]:
            s += gene + ", "
        output.append(s)
        output.append('')
        output.append("Genes that do not have the most called start:")
        s = '\u2022' + ""
        for gene in stats["no_most_called"]:
            s += gene + ", "
        output.append(s + '')
        output.append('')
        output.append("Other Starts Called:")
        for start, genes in list(stats["called_starts"].items()):
            if len(genes) == 0:
                continue
            if start != most_called_start:
                s = ''
                for gene in genes:
                    s += gene + ", "
                output.append('\u2022' + str(start) + "\t" +s +'')
                percent = float(len(genes)) / total_genes * 100
                output.append("Percent with start called: %10.4f%% \n\t" % percent)
        return output

def add_pham_no_title(args, pham_no, first_graph_path, i=""):
    # print i, type(i)
    # print first_graph_path
    packet = io.StringIO()
    can = canvas.Canvas(packet, pagesize=letter)
    width, height = letter
    # print width, height
    can.drawString(280, 750, 'Pham ' + str(pham_no))
    can.save()

    packet.seek(0)
    new_pdf = PyPDF2.PdfFileReader(packet)
    existing_pdf = PyPDF2.PdfFileReader(file(first_graph_path), 'rb')
    output = PyPDF2.PdfFileWriter()
    print(first_graph_path)
    page = existing_pdf.getPage(0)
    page.mergePage(new_pdf.getPage(0))
    output.addPage(page)
    print(utils.INTERMEDIATE_DIR)
    print("old graph?", os.path.join(args.dir, "%sPham%sGraph%s.pdf" % (args.phage+ args.one_or_all, pham_no, i)))
    outputStream = file(os.path.join(args.dir, "%sPham%sGraph%s.pdf" % (args.phage+ args.one_or_all, pham_no, i)),'wb')
    # print outputStream
    print(outputStream)
    os.remove(first_graph_path)
    output.write(outputStream)
    outputStream.close()

def combine_graphs(args, phage, pham_no, num_pages):
    merger = PyPDF2.PdfFileMerger()
    for i in range(0, num_pages+1):
        print(os.path.join(args.dir, "%sPham%sGraph%d.pdf"  % (phage + args.one_or_all, pham_no, i)))
        graph = open(os.path.join(args.dir, "%sPham%sGraph%d.pdf"  % (phage + args.one_or_all, pham_no, i)), "rb")
        merger.append(fileobj=graph)
    merger.write(open(os.path.join(args.dir, "%sPham%sGraph.pdf"  % (phage + args.one_or_all, pham_no)), "wb"))

def make_gene_track(gd_diagram, pham, gene_group, num_on_diagram, total):
    """"""
    colors = ['purple', 'red', 'green', 'orange', 'yellow', 'brown'] 
    gene = gene_group[0]
    gd_gene_track = gd_diagram.new_track(total - num_on_diagram, name='Track %s'  % (num_on_diagram+1), 
                            label=True, greytrack=1)
    gd_seq_set = gd_gene_track.new_set()
    gd_feature_set = gd_gene_track.new_set()

    start_site = gene.alignment_start_site
    start_site_feature = SeqFeature(FeatureLocation(start_site, start_site +1), 
                            strand=None)
    for feature in gene.alignment.features:
        if feature.type == 'seq':
            gd_seq_set.add_feature(feature, color='pink')
    for site in gene.alignment_candidate_starts:
        site_color = pham.total_possible_starts.index(site) % len(colors) 
        possible_site = SeqFeature(FeatureLocation(site, site ), strand=None)
        gd_feature_set.add_feature(possible_site, color=colors[site_color], 
            name=str(pham.total_possible_starts.index(site)+1), label=True)
    end_gene_feature = SeqFeature(FeatureLocation(len(gene.alignment), 
                        len(gene.alignment)+1), strand=None)
    gd_feature_set.add_feature(start_site_feature, color="blue", label=True)
    gd_feature_set.add_feature(end_gene_feature, color='purple', label=True)

def graph_start_sites(args, pham, file_path):
    """
        graphs the alignment, creates a PDF file called {Phage Name}{One or All}Pham{Pham Number}.pdf
    """

    # genes = sorted(pham.genes_in_pham.values())
    genes = pham.group_similar_genes()
    # for group in genes:
    #     print group, group.id
    if args.phage == None:
        args.phage = ""
    seq_length = len(genes[0][0].sequence.seq)
    if not args.phage:
        graph_path = os.path.join(file_path,"OnePham%sGraph_%%s.pdf" % (pham.pham_no))
    else:
        graph_path = os.path.join(file_path, "%sPham%sGraph_%%s.pdf" % (args.phage+ args.one_or_all, pham.pham_no))
    if len(genes) > 100:
        for i in range(0, int(math.ceil(len(genes)/50.0))):
            gd_diagram = GenomeDiagram.Diagram(pham.pham_no)
            final_graph_path = os.path.join(file_path,"OnePham%sGraph%d.pdf" % (pham.pham_no, i)) if not args.phage else  os.path.join(file_path, "%sPham%sGraph%d.pdf" % (args.phage+args.one_or_all, pham.pham_no, i))
            graph_path = os.path.join(file_path, "OnePham%sGraph_%d.pdf" % ( pham.pham_no, i)) if not args.phage else os.path.join(file_path,"%sPham%sGraph_%d.pdf" % ( args.phage+args.one_or_all, pham.pham_no, i))
            if check_file(final_graph_path):
                continue
            graph_path_svg = os.path.join(file_path, "%sPham%sGraph_%s.svg" % (args.phage+ args.one_or_all, pham.pham_no, i))

            for j in range(0, 50):
                if i*50 + j >= len(genes):
                    print(i * 50, + j, len(genes))
                    gd_gene_track = gd_diagram.new_track(50-j)
                    gd_feature_set = gd_gene_track.new_set()
                    empty_feature = SeqFeature(FeatureLocation(0,1), 
                            strand=None)
                    gd_feature_set.add_feature(empty_feature, color="black", label=True)
                    print('blank track added')
                else:
                    gene = genes[i*50 + j][0]
                    make_gene_track(gd_diagram, pham, genes[i*50 + j], j, 50)
            print(seq_length, i)

            gd_diagram.draw(format="linear", orientation="portrait", pagesize=letter, 
                fragments=1, start=0, end=seq_length)
            gd_diagram.write(graph_path, "PDF")
            gd_diagram.write(graph_path_svg, "SVG")

            add_pham_no_title(args, args.pham_no, graph_path, str(i))

        combine_graphs(args, args.phage, pham.pham_no, i)
    else:  
        final_graph_path = os.path.join(file_path, "Pham%sGraph_.pdf" % (pham.pham_no)) if not args.phage else os.path.join(file_path, "%sPham%sGraph_.pdf" % (args.phage+args.one_or_all, pham.pham_no))        
        # graph_path_svg = "%sPham%sGraph.svg" % (file_path+args.phage+ args.one_or_all, pham.pham_no)
        graph_path = os.path.join(file_path, "Pham%sGraph_.pdf" % (pham.pham_no)) if not args.phage else  os.path.join(file_path,"%sPham%sGraph_.pdf" % (args.phage, pham.pham_no))
        print(graph_path)
        if check_file(final_graph_path):
            pass
        else:
            gd_diagram = GenomeDiagram.Diagram(pham.pham_no)
            i = 0
            for gene_group in genes:
                print('group', i)
                make_gene_track(gd_diagram, pham, gene_group, i, len(genes))
                i += 1
            gd_diagram.draw(format="linear", orientation="portrait", pagesize=letter, 
                fragments=1, start=0, end=len(gene_group[0].alignment))
            gd_diagram.write(graph_path, "PDF")
        # gd_diagram.write(graph_path_svg, "SVG")
            add_pham_no_title(args, pham.pham_no, graph_path)

        # gd_diagram.write("%s.svg" % (file_path+pham.pham_no), "SVG")
        # gd_diagram.write("%s.eps" % (file_path+pham.pham_no), "EPS")
        # gd_diagram.write("%s.png" % (file_path+pham.pham_no), "PNG")

def make_pham_text(args, pham, pham_no, output_dir, only_pham=False):
    """
        Creates a PDF Report for the specific pham.
        From Start sites statisitics
        phage specific
    """
    if not args.phage:
        name = os.path.join(output_dir, "Pham%sText.pdf" % (pham_no))
    else:
        name = os.path.join(output_dir,"%sPham%sText.pdf" % (args.phage + args.one_or_all, pham_no))
    if check_file(name):
        return
    doc = SimpleDocTemplate(name, pagesize=letter)
    story = []
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name="paragraph"))
    styles.add(ParagraphStyle(name='Center', alignment=TA_CENTER))
    text = '<font size=14> Pham %s Report </font>' % pham_no
    story.append(Paragraph(text, styles['Center']))
    story.append(Spacer(1, 12))
    groups = pham.group_similar_genes()
    tracks_info = []
    for index in range(len(groups)):
        text = "Track %s : " % (index+1)
        text += ", ".join(gene.gene_id for gene in groups[index])
        tracks_info.append("<font size=12> "+ '\u2022'+" %s</font>" % text)
    for line in tracks_info:
        story.append(Paragraph(line, styles["Normal"]))
    story.append(Spacer(1, 12))
    if only_pham:
        start_stats = pham.stats["most_common"]
        output = output_start_sites(start_stats)
        for line in output:
            if line == '':
                story.append(Spacer(1, 12))
            text = '<font size=12> %s </font>' % line
            # if 'Genes' not in line or '':
            story.append(Paragraph(text, styles['Normal']))
            # else:
            #     story.append(Paragraph(text, styles['Normal']))
        # story.append()
        story.append(Paragraph("<font size=14>Suggested Starts:</font>", styles["Normal"]))
        suggested_start = output_suggested_starts(pham, all_genes=True)
        story.append(Spacer(1, 12))

        for line in suggested_start:
            text = '<font size=12>%s</font>' % line
            story.append(Paragraph(text, styles["Normal"]))
    else:
        story.append(Paragraph("<font size=14>Suggested Starts: </font>", styles["Normal"]))
        starts = pham.stats["most_common"]
        suggested_start = output_suggested_starts(pham, args.phage)
        story.append(Spacer(1, 12))

        for line in suggested_start:
            text = '<font size=12>%s</font>' % line
            story.append(Paragraph(text, styles["Normal"]))
        story.append(Paragraph("",styles["Normal"]))
        story.append(Paragraph("<font size=12>Gene Information:</font>", styles["Normal"]))
        pham_possible_starts = pham.total_possible_starts
        for gene in list(pham.genes.values()):
            if args.phage in gene.gene_id:
                candidate_starts = []
                for start in gene.alignment_candidate_starts:
                    candidate_starts.append((pham_possible_starts.index(start)+1, gene.alignment_index_to_coord(start)+1))
                if gene.orientation == "R":
                    story.append(Paragraph("<font size = 12> Gene: %s \n Start: %s, Stop: %s </font>" % (gene.gene_id, gene.start, gene.stop+1), styles["Normal"]))
                else:
                    story.append(Paragraph("<font size = 12> Gene: %s \n Start: %s, Stop: %s </font>" % (gene.gene_id, gene.start+1, gene.stop), styles["Normal"]))
                story.append(Paragraph("<font size = 12> Candidate Starts for %s: </font>" % (gene.gene_id), styles["Normal"]))
                story.append(Paragraph("<font size = 12>"+ str(candidate_starts) + "</font>" , styles["Normal"]))

        
    doc.build(story)


def make_pham_genome(phage_genes, phage_name, length, file_path):
    file_name = os.path.join(file_path,'%sPhamsGraph.pdf' % (phage_name))
    if check_file(file_name):
        return
    pham_colors = phams.get_pham_colors()
    gd_diagram = GenomeDiagram.Diagram(phage_name)
    gd_track = gd_diagram.new_track(1, name=phage_name, greytrack=1)
    gd_pham_set = gd_track.new_set()
    print("making genome page")
    for gene_id in sorted(phage_genes.keys()):
        phage_gene = phage_genes[gene_id]
        pham_no = phage_gene["pham_no"]
        gene = phage_gene["gene"]
        print(pham_no, gene.gene_id)
        if pham_no == None:
            pham_no = "None"
            pham_color = 'Black'
        else:
            pham_color = colors.HexColor(pham_colors[str(pham_no)])
        if gene.orientation == 'F':
            gene_location = FeatureLocation(gene.start, gene.stop)
            strand = 1
        else:
            strand = -1
            gene_location = FeatureLocation(gene.stop, gene.start)
        gene_feature = SeqFeature(gene_location, strand=strand)
        gene_number = phamgene.get_gene_number(gene.gene_id)
        # label the gene with the gene number
        gd_pham_set.add_feature(gene_feature, name=str(gene_number), label=True, 
            label_size=6, label_angle=75)
        # label gene with pham color and name
        gd_pham_set.add_feature(gene_feature, color=pham_color, name=str(pham_no), label=True, label_position='middle')
    
    print(type(length), length)
    gd_diagram.draw(format='linear', orientation='portrait', pagesize=letter, fragments=8, start=0, end=length)
    gd_diagram.write(file_name, "PDF")

def make_suggested_starts(phage_genes, phage_name, file_path):
    """
        Creates a PDF page of the suggested starts of a phage
        Genes are list in order
        {Gene Name} is a member of Pham {Number}: {Suggested Start Coordinates}
    """
    file_name = os.path.join(file_path, "%sSuggestedStarts.pdf" % (phage_name))
    if check_file(file_name):
        return
    doc = SimpleDocTemplate(file_name, pagesize=letter)
    story = []
    print("making suggested starts page")
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name="paragraph"))
    styles.add(ParagraphStyle(name='Center', alignment=TA_CENTER))
    text = '<font size=14> Suggested Start Coordinates</font>'
    story.append(Paragraph(text, styles['Center']))
    story.append(Spacer(1, 12))
    for gene_id in sorted(phage_genes.keys()):
        phage_gene = phage_genes[gene_id]
        pham = phage_gene["pham_no"]
        gene = phage_gene["gene"]
        suggested_start = phage_gene["suggested_start"]
        if pham == None:
            text = '<font size=12> %s is not a member of an existing Pham </font>' % (gene.gene_id)
        else:
            text = '<font size=12> %s is a member of Pham %s:  %s </font>' % (gene.gene_id, pham, suggested_start)
        story.append(Paragraph(text, styles['Normal']))
    doc.build(story)

def make_fasta_file(genes, fasta_file):
    count = SeqIO.write(genes, fasta_file, 'fasta')


def main():
    args = parse_arguments()
    print("hi from making files", args.make)
    if 'graph' in args.make:
        print(args.pickle_file)
        pham = pickle.load(open(args.pickle_file.strip('"'), 'rb'))
        graph_start_sites(args, pham, args.dir)

    if 'starts' in args.make:
        phage_genes = pickle.load(open(args.pickle_file.strip('"'), 'rb'))

        make_suggested_starts(phage_genes, args.phage, args.dir)

    if 'genome' in args.make:
        phage = pickle.load(open(args.pickle_file.strip('"'), 'rb'))
        make_pham_genome(phage, args.phage, args.phage_length, args.dir)
        make_suggested_starts(phage, args.phage, args.dir)

    if 'text' in args.make:
        print(args.pickle_file)
        pham = pickle.load(open(args.pickle_file.strip('"'), 'rb'))
        graph_start_sites(args, pham, args.dir)
        print("phage", args.phage)
        if not args.phage:
            print("hello no phage")
            make_pham_text(args, pham, args.pham_no, args.dir, only_pham=True)
        else:
            make_pham_text(args, pham, args.pham_no, args.dir)

    if 'fasta' in args.make:
        pass
        # pickle_file = args.file
        # genes = cPickle.load(open(args.pickle_file.strip('"'), 'rb'))
        # make_fasta_file(genes, (args.dir + '.fasta'))

if __name__ == "__main__":
    main()
