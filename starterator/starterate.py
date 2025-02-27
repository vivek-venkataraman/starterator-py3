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
# Starterate function 

import argparse
from starterator import utils
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import letter, A4
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib import colors
import PyPDF2
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Graphics import GenomeDiagram
from reportlab.lib import colors
from reportlab.lib.pagesizes import letter, A4
# from Bio.Blast import NCBIXML
# from Bio.Blast.Applications import NcbiblastpCommandline as Blastp
# from Bio.Blast.Applications import BlastallCommandline
import MySQLdb
import subprocess
import os, sys
import math
import pickle
# from gi.repository import Gtk, Gdk, GObject
import getpass
from . import report
from . import phamgene
'''
def gui():
    GObject.threads_init()
    Gdk.threads_init()
    win = StarteratorWindow()
    win.connect('delete-event', Gtk.main_quit)
    win.show_all()
    Gdk.threads_enter()
    Gtk.main()
    Gdk.threads_leave()
'''
def get_output_one_pham(pham, pham_no, config):
    """
        Creates a PDF Report for the specific pham.
        From Start sites statisitics
    """
    output_dir = config['intermediate_file_dir']
    doc = SimpleDocTemplate("%s%sPham%sText.pdf" % (output_dir, phage+one_or_all, pham_no), pagesize=letter)
    story = []
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name="paragraph"))
    styles.add(ParagraphStyle(name='Center', alignment=TA_CENTER))
    text = '<font size=14> Pham %s Report </font>' % pham_no
    story.append(Paragraph(text, styles['Center']))
    story.append(Spacer(1, 12))
    output = pham.output_start_sites()
    for line in output:
        if line == '':
            story.append(Spacer(1, 12))
        text = '<font size=12> %s </font>' % line
        story.append(Paragraph(text, styles['Normal']))
    suggested_start = pham.output_suggested_starts()
    story.append(Spacer(1, 12))
    for line in suggested_start:
        text = '<font size=12>%s</font>' % line
        story.append(Paragraph(text, styles["Normal"]))
    doc.build(story)

def get_arguments():
    parser = argparse.ArgumentParser(prog='starterate.py', usage='Phameratored Phage Report: %(prog)s -p {Phage Name}\n'
            +'One Gene of Phameratored Phage Report:  %(prog)s -p {Phage Name} -n {Pham Number}\n'
            + 'Unphameratored Phage Report:  %(prog)s -p {Phage Name} -u True -f {Path to DNAMaster profile file}\n'
            + 'One Gene of Unphameratored Phage Report:  %(prog)s -p {Phage Name} -u True -s {Start of Gene} -t {Stop of Gene} -o {Orientation of Gene} -g {Number of Gene}')
    parser.add_argument('-n', '--pham_no', default = -1, help='Number of the Pham. For case when want report of a phameratored phage gene.')
    parser.add_argument('-p' , '--phage', default=None, help='The Phamerator database Phage Name. Always needed')
    parser.add_argument('-u', '--unphamed', type=bool, default=False,
                 help='Boolean. If phage has been phameratored: False.'
                +' If phage is unphameratored: True. For use when want report with an unphameratored phage')
    parser.add_argument('-s', '--given_start', type=int, default=-1, 
        help= 'The start of a gene. For case when want report of one gene of an unphameratored phage. ')
    parser.add_argument('-t', '--given_stop', type=int, 
        help= 'The stop of a gene. For case when want report of one gene of an unphameratored phage.')
    parser.add_argument('-o', '--given_orientation',
         help='The orientation of a gene. For case when want report of one gene of an unphameratored phage.')
    parser.add_argument('-g', '--gene_number', default=-1,
         help='The number of a gene. For case when want report of one gene of an unphameratored phage.')
    parser.add_argument('-d', '--profile', help='Path to a DNAMaster profile. For case when want whole report of an unphameratored phage')
    parser.add_argument('-f', '--fasta', help='Path to Fasta File')
    return parser.parse_args()



def starterate(info, gui=None, event=None):
    global one_or_all, phage, protein_db, output_dir, final_dir
    one_or_all = 'All' if info['all'] else 'One'
    phage = info['phage']
    # print phage
    protein_db = utils.PROTEIN_DB + 'ProteinsDB'
    output_dir = utils.INTERMEDIATE_DIR
    final_dir = utils.FINAL_DIR
    if info['all'] and info['phamerated']:
        phage = report.PhageReport(phage, gui=gui, event=event)
        final_file, short_final = phage.final_report()
    elif info['all'] and not info['phamerated']:
        # phams.update_protein_db(db, config)
        phage = report.UnPhamPhageReport(phage, fasta_file=info['fasta'], profile_file=info['profile'], gui=gui, event=event)
        final_file, short_final = phage.final_report()
    elif not info['all'] and info['phamerated'] and not info['pham']:
        gene = report.GeneReport( info['phage'], number=info['gene_no'])
        gene.get_pham()
        gene.make_report()
        final_file, s = gene.merge_report()
    elif not info['all'] and not info['phamerated'] and not info['pham']:
        # phams.update_protein_db(db, config)
        gene = report.GeneReport(info['phage'], info["gene_no"], fasta_file=info['fasta'])
        gene.make_unpham_gene(int(info['start']), int(info['stop']), info['orientation'])
        gene.make_report()
        final_file, s = gene.merge_report()
    else:
        # Pham without phage'
        pham = report.PhamReport(info['pham'])
        final_file,s = pham.final_report()
    return final_file

def main():
    args = get_arguments()
    config = utils.get_config()
    print((config["count"]))
    phamgene.check_protein_db(config["count"])
    # --Phamerated and only one gene
    if args.gene_number != -1 and args.phage != None and args.unphamed == False:
        gene = report.GeneReport(args.phage, args.gene_number, True)
        print(gene)
        gene.get_pham()
        gene.make_report()
        final_file, s = gene.merge_report()

    # --Unphameratored Phage with only one gene
    elif args.given_start > -1 and args.phage != None and args.unphamed == True:
        # given start and stop coordinates and orientation
        one_or_all = 'One'
        given_start = args.given_start
        given_stop = args.given_stop
        given_orientation = args.given_orientation
        gene_name = args.phage + '_' + str(args.gene_number)
        gene = report.GeneReport(args.phage, args.gene_number, fasta_file=args.fasta)
        gene.make_unpham_gene(given_start, given_stop, given_orientation)
        print(gene)
        gene.make_report()
        final_file, s = gene.merge_report()

    # --Phameratored or Unphameratored Phages with all genes
    elif args.pham_no == -1 and args.phage != None and args.unphamed == False:
        phage = report.PhageReport(args.phage, gui=None)
        final_file, short_final = phage.final_report()

    elif args.pham_no == -1 and args.phage != None and args.unphamed == True:
        phage = report.UnPhamPhageReport(args.phage, fasta_file=args.fasta, profile_file=args.profile, gui=None)
        final_file, short_final = phage.final_report()
    elif args.phage == None:
        pham = report.PhamReport(args.pham_no)
        final_file, short_final = pham.final_report()
    # clean_up_files()
    # email_final_report(args.email, short_final)


if __name__ == '__main__':
    try:
        # print utils.DIR_PATH ???
        main() 
    except Exception as e:
        # send_error_email(e)
        raise
    finally:
        pass
         # clean_up_files()
