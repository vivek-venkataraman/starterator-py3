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
# Functions that create the PDF outputs

import pickle
import argparse
import time
from Bio.Graphics import GenomeDiagram
from Bio.SeqFeature import SeqFeature, FeatureLocation
import reportlab.lib.pagesizes
from reportlab.pdfgen import canvas
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib import colors
import PyPDF2
from Bio import SeqIO
import math
from io import StringIO, BytesIO
from starterator.utils import *
from starterator.phams import *
from starterator.phamgene import *
import os
# from phage import
# from reportlab.lib import colors


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--pham_no', default=-1)
    parser.add_argument('-p', '--phage', default=None)
    parser.add_argument('-d', '--dir', default=None)
    parser.add_argument('-a', '--one_or_all', default='One')
    parser.add_argument('-f', '--pickle_file', default=None)
    parser.add_argument('-l', '--phage_length', type=int, default=-1)
    parser.add_argument('-m', '--make', default=None)

    return parser.parse_args()


def check_file(file_name):
    try: 
        f = open(file_name, "r")
        f.close()
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
    if not all_genes:
        for gene in pham.genes.values():
            if phage in gene.gene_id:
                coord_suggestion = gene.suggested_start["most_called"]
                output.append("" + gene.gene_id + ", " + str(coord_suggestion))
    else:
        for gene in pham.genes.values():
            coord_suggestion = gene.suggested_start["most_called"]
            output.append(gene.gene_id + ", " + str(coord_suggestion))
    return output


def output_start_sites(stats):
        """
            Writes a report of the start sites statistics
            Returns a list of strings with information about
            start sites for this pham  
        """
        # most_called_start = stats["most_called_start"]
        most_annotated_start = stats["most_annotated_start"]
        annotated_count = stats["annotCount"]
        # noinspection PyListCreation
        output = []
        output.append("")

        # start section on summary of annotations:
        output.append('<b>Summary of Final Annotations (See graph section above for start numbers):</b>')

        if annotated_count > 0:
            annotated_with_most_annotated_called = \
                [g.full_name for g in stats["annot_list"] if g.full_name in stats["called_starts"][most_annotated_start]]
            # annotated_with_most_predicted_called = \
            #    [g.gene_id for g in stats["draft_list"] if g.gene_id in stats["called_starts"][most_called_start]]

            s = 'The start number called the most often in the published annotations is %s,' % str(most_annotated_start)
            s += ' it was called in %s of the %s non-draft genes in the pham.' % \
                (str(len(annotated_with_most_annotated_called)), str(annotated_count))
            output.append(s)

            output.append('')
            # output.append('Called info: "Most Called" Start is %s, called in %s of all %s genes in the pham.'
            #               % (str(stats["most_called_start"]), str(calledCount), str(total_genes)))
            # percent_with_most_annotated = (float(len(stats["most_called"]))
            #                             /total_genes *100 )
            #
            # output.append('Percent of genes that begin at the "Most Annotated" start: %10.1f%%'
            #                 % percent_with_most_called )

            called_most_annotated = [gene for gene in stats["called_starts"][most_annotated_start]]
            called_most_annotated.sort()
            output.append('Genes that call this "Most Annotated" start:')
            s = '\u2022' + ' '
            for gene in called_most_annotated:
                s += gene + ", "
            output.append(s)
            output.append('')

            have_most_annotated = [gene for gene in stats["most_not_annotated"]]
            have_most_annotated.sort()
            output.append('Genes that have the "Most Annotated" start but do not call it:')
            s = '\u2022' + ' '
            for gene in have_most_annotated:
                s += gene + ", "
            output.append(s)
            output.append('')

            has_not_most_annotated = [gene for gene in stats["no_most_annot"]]
            has_not_most_annotated.sort()
            output.append('Genes that do not have the "Most Annotated" start:')
            s = '\u2022' + ' '
            for gene in has_not_most_annotated:
                s += gene + ", "
            output.append(s + '')
            output.append('')
        else:
            output.append("This pham is comprised of all draft annotations. There are no annotations to summarize.")

        # start section summary of start sites by number
        output.append("<b>Summary by start number:</b>")
        output.append('')

        for start, genes in stats["called_starts"].items():
            if len(genes) == 0:
                continue
            output.append("Start %s:" % str(start))
            presence = len(stats["possible"][start])
            percent_present = float(presence) / stats['phamCount'] * 100
            output.append('\u2022' + " Found in %s of %s (%10.1f%% ) of genes in pham\t" %
                          (str(presence), str(stats['phamCount']), percent_present))

            if start in stats['annot_counts'].keys():
                output.append('\u2022' + " Manual Annotations of this start: %s of %s" %
                              (str(stats['annot_counts'][start]), str(annotated_count)))
            else:
                output.append('\u2022' + " No Manual Annotations of this start. " )

            percent_called = float(len(genes)) / presence * 100
            output.append('\u2022' + " Called %10.1f%% of time when present \n\t" % percent_called)

            cluster_dict = {}
            for p_gene in stats['annot_list']:
                cluster_dict[p_gene.full_name] = p_gene.subcluster

            for p_gene in stats['draft_list']:
                cluster_dict[p_gene.full_name] = p_gene.subcluster

            genes.sort()
            s = ''
            cluster_start = []
            for gene in genes:
                cluster_start.append(cluster_dict[gene])
                s += gene + " (" + cluster_dict[gene] + "), "
            output.append('\u2022' + " Phage (with cluster) where this start called:\t" + s + '')
            output.append('')

        output.append("<b>Summary by clusters:</b>")
        output.append('')

        if int(len(stats['clusters_present'])) == 1:
            cluster_string = str(stats['clusters_present'])
            s = "There is one cluster represented in this pham: %s" % str(stats['clusters_present'])[2:-2]
        else:
            cluster_string = ""
            for s in stats['clusters_present']:
                cluster_string += s + ", "
            s = "There are %s clusters represented in this pham: %s" % (len(stats['clusters_present']), cluster_string)

        output.append(s)
        output.append('')

        annotated_genes = [g.gene_id for g in stats['annot_list']]

        for cluster in sorted(stats['clusters_present']):
            if cluster == 'singleton':
                continue

            count_MA = 0
            for geneid in annotated_genes:
                if geneid in stats['genes_by_cluster'][cluster]:
                    count_MA += 1
            if count_MA > 0:
                output.append("Info for manual annotations of cluster %s:" % cluster)
                annotated_cluster_starts = [ph.alignment_start_num_called for ph in stats['annot_list'] if ph.subcluster == cluster]
                start_counts = dict([(x,annotated_cluster_starts.count(x)) for x in set(annotated_cluster_starts)])
                starts_present = sorted(start_counts.keys())
                for start in starts_present:
                    count = start_counts[start]
                    if count > 1:
                        s = "Start number %s was manually annotated %s times for cluster %s." % (start, count, cluster)
                    else:
                        s = "Start number %s was manually annotated 1 time for cluster %s." % (start, cluster)
                    output.append('\u2022' + s)

            output.append('')

        return output


# noinspection PyListCreation
def output_start_sites_by_phage(stats, genelist):
    """
        Writes a report of the start sites statistics
        Returns a list of strings with information about
        start sites for this pham
    """
    # get phage name from the 1st entry in the gene list
    # phage = genelist[0].phage_name

    # Older variables which could be used:
    # sort genelist on gene_id
    # most_called_start = stats["most_called_start"]
    # total_genes = stats["phamCount"]
    # draftCount = stats["draftCount"]
    # calledCount = len(stats["called_starts"][most_called_start])

    # If cluster not known or a singleton, use the standard pham report not this one
    if genelist[0].cluster == 'Unassigned':
        return output_start_sites(stats)


    genelist.sort(key=lambda x: x.start)
    most_annotated_start = stats["most_annotated_start"]
    annotated_count = stats["annotCount"]
    output = []
    output.append("")

    # start section on summary of annotations:
    output.append('Summary of Final Annotations (See graph section above for start numbers):')

    if annotated_count > 0:
        annotated_with_most_annotated_called = \
            [g.full_name for g in stats["annot_list"] if g.full_name in stats["called_starts"][most_annotated_start]]
        # annotated_with_most_predicted_called = \
        #    [g.gene_id for g in stats["draft_list"] if g.gene_id in stats["called_starts"][most_called_start]]

        s = 'The start number called the most often in the published annotations is %s,' % str(most_annotated_start)
        s += ' it was called in %s of the %s non-draft genes in the pham.' % \
             (str(len(annotated_with_most_annotated_called)), str(annotated_count))
        output.append(s)
        output.append('')

        called_most_annotated = [gene for gene in stats["called_starts"][most_annotated_start]]
        have_most_annotated = [gene for gene in stats["most_not_annotated"]]
        has_not_most_annotated = [gene for gene in stats["no_most_annot"]]

        for gene in genelist:
            if gene.gene_id in called_most_annotated:
                output.append('<b>%s did call</b> the "Most Annotated" start (%s).' %
                              (gene.full_name, str(gene.alignment_start_num_called)))
            elif gene.gene_id in have_most_annotated:
                output.append('<b>%s has but does not call</b> the "Most Annotated" start, calling start %s instead.' %
                              (gene.full_name, str(gene.alignment_start_num_called)))
            elif gene.gene_id in has_not_most_annotated:
                output.append('<b>%s does not have</b> the "Most Annotated" start, calling start %s instead. ' %
                              (gene.full_name, str(gene.alignment_start_num_called)))

        output.append('')

    else:
        output.append("This pham is comprised of all draft annotations. There are no annotations to summarize.")
        output.append('')

    # start section summary of start sites by number
    output.append("<b>Summary by start number:</b>")
    for gene in genelist:
        if len(gene.alignment_annot_start_nums) > 1:
            output.append("%s has %d starts with manual annotations (numbers: %s; called: %s times)." %
                          (gene.full_name, len(gene.alignment_annot_start_nums),
                           str(gene.alignment_annot_start_nums)[1:-1], str(gene.alignment_annot_start_counts)[1:-1]))
        elif len(gene.alignment_annot_start_nums) == 1:
            output.append("%s has 1 start with manual annotations (number %s; called: %s times)." %
                          (gene.full_name, str(gene.alignment_annot_start_nums)[1:-1], str(gene.alignment_annot_start_counts)[1:-1]))
        else:
            output.append("%s has no starts with manual annotations in other genes." % gene.full_name)

    output.append('')
    if len(genelist) > 1:
        all_annotated_starts = []
        for gene in genelist:
            all_annotated_starts += gene.alignment_annot_start_nums

        all_annotated_starts = list(set(all_annotated_starts))
    else:
        all_annotated_starts = genelist[0].alignment_annot_start_nums

    for start in all_annotated_starts:
        output.append("Start %s:" % str(start))
        presence = len(stats["possible"][start])
        percent_present = float(presence) / stats['phamCount'] * 100
        output.append('\u2022' + " Found in %s of %s (%10.1f%% ) of genes in pham\t" %
                      (str(presence), str(stats['phamCount']), percent_present))

        percent_called = float(len(stats['called_starts'][start])) / presence * 100
        output.append('\u2022' + " Called %10.1f%% of time when present \n\t" % percent_called)

        genes = stats["called_starts"][start]

        cluster_dict = {}
        for p_gene in stats['annot_list']:
            cluster_dict[p_gene.full_name] = p_gene.subcluster

        for p_gene in stats['draft_list']:
            cluster_dict[p_gene.full_name] = p_gene.subcluster

        genes.sort()
        s = ''
        cluster_start = []
        for gene in genes:
            cluster_start.append(cluster_dict[gene])
            s += gene + " (" + cluster_dict[gene] + "), "
        output.append('\u2022' + " Phage (with cluster) where this start called:\t" + s + '')
        output.append('')

    output.append("<b>Summary by clusters:</b>")
    output.append('')

    if int(len(stats['clusters_present'])) == 1:
        cluster_string = str(stats['clusters_present'])
        s = "There is one cluster represented in this pham: %s" % str(stats['clusters_present'])[2:-2]
    else:
        cluster_string = ""
        for s in stats['clusters_present']:
            cluster_string += s + ", "
        s = "There are %s clusters represented in this pham: %s" % (len(stats['clusters_present']), cluster_string)

    output.append(s)
    output.append('')

    annotated_genes = [g.gene_id for g in stats['annot_list']]

    cluster = genelist[0].subcluster

    if cluster == 'singleton':
        output.append("This phage is a singleton, no other cluster members to compare.")
    else:
        count_MA = 0
        for geneid in annotated_genes:
            if geneid in stats['genes_by_cluster'][cluster]:
                count_MA += 1
        if count_MA > 0:
            output.append("Info for manual annotations of cluster %s:" % cluster)
            annotated_cluster_starts = [ph.alignment_start_num_called for ph in stats['annot_list'] if ph.subcluster == cluster]
            start_counts = dict([(x,annotated_cluster_starts.count(x)) for x in set(annotated_cluster_starts)])
            starts_present = sorted(start_counts.keys())
            for start in starts_present:
                count = start_counts[start]
                if count > 1:
                    s = "Start number %s was manually annotated %s times for cluster %s." % (start, count, cluster)
                else:
                    s = "Start number %s was manually annotated 1 time for cluster %s." % (start, cluster)
                output.append('\u2022' + s)

    output.append('')

    return output


def add_pham_no_title(args, pham_no, first_graph_path, i="", zoom=False):
    # print i, type(i)
    # print first_graph_path
    packet = BytesIO()
    can = canvas.Canvas(packet, pagesize=reportlab.lib.pagesizes.letter)
    # width, height = reportlab.lib.pagesizes.letter
    # print width, height
    if zoom:
        can.drawString(250, 750, 'Zoomed Pham ' + str(pham_no))
    else:
        can.drawString(280, 750, 'Pham ' + str(pham_no))

    can.save()

    packet.seek(0)
    new_pdf = PyPDF2.PdfReader(packet)
    existing_pdf = PyPDF2.PdfReader(open(first_graph_path, 'rb'))
    output = PyPDF2.PdfWriter()
    # print first_graph_path
    page = existing_pdf.pages[0]
    page.merge_page(new_pdf.pages[0])
    output.add_page(page)
    # print utils.INTERMEDIATE_DIR
    # print "old graph?", os.path.join(args.dir, "%sPham%sGraph%s.pdf" % (args.phage + args.one_or_all, pham_no, i))
    output_strm = open(os.path.join(args.dir, "%sPham%sGraph%s.pdf" % (args.phage + args.one_or_all, pham_no, i)), 'wb')
    # print outputStream
    # print output_strm
    os.remove(first_graph_path)
    output.write(output_strm)
    output_strm.close()


def combine_graphs(args, phage, pham_no, num_pages):
    writer = PyPDF2.PdfWriter()
    for j in range(0, num_pages + 1):
        # print os.path.join(args.dir, "%sPham%sGraph%d.pdf" % (phage + args.one_or_all, pham_no, j))
        with open(os.path.join(args.dir, "%sPham%sGraph%d.pdf" % (phage + args.one_or_all, pham_no, j)), "rb") as graph:
            reader = PyPDF2.PdfReader(graph)
            for page in reader.pages:
                writer.add_page(page)
    with open(os.path.join(args.dir, "%sPham%sGraph.pdf" % (phage + args.one_or_all, pham_no)), "wb") as output_file:
        writer.write(output_file)


def make_gene_track(gd_diagram, pham, gene_group, num_on_diagram, total, seqColor):
    """

    :param gd_diagram:
    :param pham:
    :param gene_group:
    :param num_on_diagram:
    :param total:
    :return:
    """

    start_bar_colors = ['purple', 'red', 'lightblue', 'orange', 'tan', 'brown']
    gene = gene_group[0]

    # change track_name to name of fist gene in list
    track_name = str(num_on_diagram+1)
    track_name += ": "
    track_name += gene.full_name
    if len(gene_group) > 1:
        track_name += " + "
        track_name += str(len(gene_group)-1)
    num_on_diagram = num_on_diagram%50 # change from label to positional info
    gd_gene_track = gd_diagram.new_track(total - num_on_diagram, label=True,
                                         name=track_name, greytrack=1, greytrack_labels=1)
    gd_seq_set = gd_gene_track.new_set()
    gd_feature_set = gd_gene_track.new_set()

    start_site = gene.alignment_start_site
    start_site_feature = SeqFeature(FeatureLocation(start_site, start_site + 1))
    for feature in gene.alignment.features:
        if feature.type == 'seq':
            if seqColor % 2 == 0:
                trackColor = (253, 191, 203)
            else:
                trackColor = (237, 205, 223)
            gd_seq_set.add_feature(feature, color=trackColor)
    for site in gene.alignment_candidate_starts:
        site_color = pham.total_possible_starts.index(site) % len(start_bar_colors)
        possible_site = SeqFeature(FeatureLocation(site, site))
        gd_feature_set.add_feature(possible_site, color=start_bar_colors[site_color],
                                   name=str(pham.total_possible_starts.index(site) + 1), label=True)
    end_gene_feature = SeqFeature(FeatureLocation(len(gene.alignment), 
                                  len(gene.alignment)+1))

    # draw blue called start only if non-draft gene in gene group, if all draft use yellow

    all_draft_status = True
    for gene in gene_group:
        all_draft_status = all_draft_status and gene.draftStatus

    if all_draft_status:
        startcolor = "yellow"
    else:
        startcolor = "green"

    gd_feature_set.add_feature(start_site_feature, color=startcolor, label=True)
    gd_feature_set.add_feature(end_gene_feature, color='purple', label=True)


def graph_start_sites(args, pham, file_path):
    """
        graphs the alignment, creates a PDF file called {Phage Name}{One or All}Pham{Pham Number}.pdf
    """

    # genes = sorted(pham.genes_in_pham.values())
    if args.phage is None:
        args.phage = ""
        order_by = None
    else:
        order_by = args.phage

    genes = pham.group_similar_genes(order_by)

    # check for especially long upstream sequences, and define genome diagram left boundary if found
    min_start_coord = min(pham.total_possible_starts)
    max_start_coord = max(pham.total_possible_starts)

    annotated_start_nums = set()
    for gene_list in genes:
        annotated_start_nums.add(gene_list[0].alignment_start_num_called)

    min_annot_num = min(annotated_start_nums)
    max_annot_num = max(annotated_start_nums)
    min_annot_coord = pham.total_possible_starts[min_annot_num - 1]
    max_annot_coord = pham.total_possible_starts[max_annot_num - 1]
    # range of interest (roi) will be range of annots plus one more on each side
    if min_annot_num == 1:
        min_roi_num = 1
    else:
        min_roi_num = min_annot_num - 1

    if max_annot_num == len(pham.total_possible_starts):
        max_roi_num = len(pham.total_possible_starts)
    else:
        max_roi_num = max_annot_num + 1

    roi_size = pham.total_possible_starts[max_roi_num - 1] - pham.total_possible_starts[min_roi_num - 1]
    if roi_size < 3:
        roi_size = 3
    annots_in_roi = range(min_roi_num, max_roi_num + 1)
    max_right_boundary = len(genes[0][0].alignment)

    # test if zoom in needed (i.e. if there are a large number of different starts in a small fraction of a track)
    all_starts_range = (min([max_right_boundary, max_start_coord + 30])) - (max([0, min_start_coord - 30]))
    default_fraction_of_track_with_annots = float(roi_size)/all_starts_range

    max_possible_in_annotation_range = 0
    starts_in_annot_range = set(annots_in_roi)
    for gene_list in genes:
        gene = gene_list[0]
        count = len(set(gene.alignment_candidate_start_nums) & starts_in_annot_range)
        if count > max_possible_in_annotation_range:
            max_possible_in_annotation_range = count

    # if possible set zoom so that there is about 1 start per <scale>% of track
    scale = 1.1

    # test is there are a large number of annots, if so just zoom in to region around them
    if max_possible_in_annotation_range > (100 / scale) - 2:
        should_zoom = True
        right_draw_boundary =  min([max_right_boundary, max_annot_coord + 100])
        left_draw_boundary = max([0, min_annot_coord - 100])

        # threshold for deciding 'too packed' is more than an average of 1 start per 1%
    elif max_possible_in_annotation_range > default_fraction_of_track_with_annots * 100:
        should_zoom = True

        # we have len(annots_in_roi) number of annotations of intest [aoi] that span roi_size number of bases
        # we want this number of bases to take up <len(annots_in_roi) * scale>% of the track
        fraction_with_aoi = float(max_possible_in_annotation_range * scale) / 100
        # so the total should span this number of bases
        total_width = int(roi_size / fraction_with_aoi)
        flank = total_width - roi_size
        if int(flank / 2) + max_annot_coord > max_right_boundary:
            right_draw_boundary = max_right_boundary
            left_draw_boundary = max([0, min_annot_coord - (flank - (max_right_boundary - max_annot_coord))])
        elif min_annot_coord - int(flank / 2) < 1:
            left_draw_boundary = 0
            right_draw_boundary = min([max_right_boundary, max_annot_coord + (flank - min_annot_coord) ])
        else:
            right_draw_boundary = min([max_right_boundary, max_annot_coord + int(flank / 2)])
            left_draw_boundary = max([0, min_annot_coord - int(flank / 2)])

    else:
        # use default track boundies are +/- 30 bases surounding all possible starts
        should_zoom = False
        left_draw_boundary = max([0, min_start_coord - 30])
        right_draw_boundary = min([len(genes[0][0].alignment), max_start_coord + 30])

    if len(genes) > 50:
        seqColor = 0
        for i in range(0, int(math.ceil(len(genes)/50.0))):
            gd_diagram = GenomeDiagram.Diagram(pham.pham_no)
            if not args.phage:
                final_graph_path = os.path.join(file_path, "OnePham%sGraph%d.pdf" % (pham.pham_no, i))
                graph_path = os.path.join(file_path, "OnePham%sGraph_%d.pdf" % (pham.pham_no, i))

            else:
                final_graph_path = os.path.join(file_path,
                                                "%sPham%sGraph%d.pdf" % (args.phage + args.one_or_all, pham.pham_no, i))
                graph_path = os.path.join(file_path,
                                          "%sPham%sGraph_%d.pdf" % (args.phage + args.one_or_all, pham.pham_no, i))

            if check_file(final_graph_path):
                continue

            for j in range(0, 50):
                if i*50 + j >= len(genes):
                    # print i * 50, + j, len(genes)
                    gd_gene_track = gd_diagram.new_track(50-j)
                    gd_feature_set = gd_gene_track.new_set()
                    empty_feature = SeqFeature(FeatureLocation(0, 1))
                    gd_feature_set.add_feature(empty_feature, color="black", label=True)
                else:
                    if i + j > 0: # i.e. not the first track
                        if genes[i*50 + j][0].subcluster != genes[i*50 +j - 1][0].subcluster:
                            seqColor += 1
                    gene = genes[i*50 + j][0]
                    make_gene_track(gd_diagram, pham, genes[i*50 + j], i*50 + j, 50, seqColor)

            gd_diagram.draw(format="linear", orientation="portrait", pagesize=reportlab.lib.pagesizes.letter,
                            fragments=1, start=left_draw_boundary, end=right_draw_boundary)
            gd_diagram.write(graph_path, "PDF")
            # gd_diagram.write(graph_path_svg, "SVG")

            add_pham_no_title(args, args.pham_no, graph_path, str(i), should_zoom)

        combine_graphs(args, args.phage, pham.pham_no, i)
    else:
        if not args.phage:
            final_graph_path = os.path.join(file_path, "Pham%sGraph_.pdf" % pham.pham_no)
        else:
            final_graph_path = os.path.join(file_path, "%sPham%sGraph_.pdf" % (args.phage+args.one_or_all, pham.pham_no))
        # graph_path_svg = "%sPham%sGraph.svg" % (file_path+args.phage+ args.one_or_all, pham.pham_no)
        if not args.phage:
            graph_path = os.path.join(file_path, "Pham%sGraph_.pdf" % pham.pham_no)
        else:
            graph_path = os.path.join(file_path, "%sPham%sGraph_.pdf" % (args.phage, pham.pham_no))
        # print "making_files.graph_start_sites: path to graph is " + str(graph_path)

        if check_file(final_graph_path):
            pass
        else:
            gd_diagram = GenomeDiagram.Diagram(pham.pham_no)
            i = 0
            seqColor = 0

            for gene_group in genes:
                if i > 0:
                    if genes[i][0].subcluster != genes[i-1][0].subcluster:
                        seqColor += 1
                # print 'making_files.graph_start_sites: adding group ' + str(i)
                make_gene_track(gd_diagram, pham, gene_group, i, len(genes), seqColor)
                i += 1
            gd_diagram.draw(format="linear", orientation="portrait", pagesize=reportlab.lib.pagesizes.letter,
                            fragments=1, start=left_draw_boundary, end=right_draw_boundary)
            gd_diagram.write(graph_path, "PDF")
        # gd_diagram.write(graph_path_svg, "SVG")
            add_pham_no_title(args, pham.pham_no, graph_path, zoom=should_zoom)

        # gd_diagram.write("%s.svg" % (file_path+pham.pham_no), "SVG")
        # gd_diagram.write("%s.eps" % (file_path+pham.pham_no), "EPS")
        # gd_diagram.write("%s.png" % (file_path+pham.pham_no), "PNG")


def make_pham_text(args, pham, pham_no, output_dir, only_pham=False):
    """
        Creates a PDF Report for the specific pham.
        From Start sites statistics
        phage specific
    """
    if not args.phage:
        order_by = None
        name = os.path.join(output_dir, "Pham%sText.pdf" % pham_no)
    else:
        order_by = args.phage
        name = os.path.join(output_dir, "%sPham%sText.pdf" % (args.phage + args.one_or_all, pham_no))
    if check_file(name):
        return
    doc = SimpleDocTemplate(name, pagesize=reportlab.lib.pagesizes.letter)
    story = []
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name="paragraph"))
    styles.add(ParagraphStyle(name='Center', alignment=TA_CENTER))

    # increase leading a bit
    styles["Normal"].leading = 14

    note = '<font size=12>Note: Tracks are now grouped by subcluster and scaled. Switching in subcluster is ' \
           'indicated by changes in track color. Track scale is now set by default to display the region 30 ' \
           'bp upstream of start 1 to 30 bp downstream of the last possible start. If this default region ' \
           'is judged to be packed too tightly with annotated starts, the track will be further scaled ' \
           'to only show that region of the ORF with annotated starts. This action will be indicated by ' \
           'adding "Zoomed" to the title. For starts, yellow indicates the location of called starts ' \
           'comprised solely of Glimmer/GeneMark auto-annotations, green indicates the location of called ' \
           'starts with at least 1 manual gene annotation. </font>'

    story.append(Paragraph(note, styles["Normal"]))
    story.append(Spacer(1, 12))

    text = '<font size=14> Pham %s Report </font>' % pham_no
    story.append(Paragraph(text, styles['Center']))
    story.append(Spacer(1, 12))
    current_date = time.strftime("%x")
    db_version = get_version()
    run_date = '<font size=12>This analysis was run %s on database version %s. </font>' % (current_date, db_version)
    story.append(Paragraph(run_date, styles["Normal"]))
    story.append(Spacer(1, 12))

    if pham.count != len(pham.genes):
        warning_text = '<font size=12>WARNING: Pham size does not match number of genes in report. Either ' \
                       'unphamerated genes have been added (by you) ' \
                       'or starterator has removed genes due to invalid start codon. </font>'
        story.append(Paragraph(warning_text, styles["Normal"]))
        story.append(Spacer(1, 12))

    pham_count = len(pham.genes.values())
    draft_count = sum(1 for g in pham.genes.values() if g.draftStatus)
    annot_count = pham_count - draft_count
    summary_text = "<font size=12>Pham number %s has %s members, %s are drafts.</font>" % \
                   (pham_no, len(pham.genes), draft_count)
    story.append(Paragraph(summary_text, styles["Normal"]))
    story.append(Spacer(1, 12))

    story.append(Paragraph('<font size=12>Phages represented in each track:</font>', styles["Normal"]))

    groups = pham.group_similar_genes(order_by)
    tracks_info = []
    for index in range(len(groups)):
        text = "Track %s : " % (index + 1)
        text += ", ".join(gene.full_name for gene in groups[index])
        tracks_info.append("<font size=12> " + '\u2022' + " %s</font>" % text)
    for line in tracks_info:
        story.append(Paragraph(line, styles["Normal"]))
    story.append(Spacer(1, 12))

    start_stats = pham.stats["most_common"]
    start_stats["phamCount"] = pham_count
    start_stats["annotCount"] = annot_count
    start_stats["draftCount"] = draft_count

    if only_pham:
        output = output_start_sites(start_stats)
    else:
        all_genes = pham.genes.values()
        genes_in_phage = []
        for gene in all_genes:
            if args.phage.lower() == gene.phage_name.lower():
                genes_in_phage.append(gene)

        output = output_start_sites_by_phage(start_stats, genes_in_phage)

    for i, line in enumerate(output):
        if line == '':
            story.append(Spacer(1, 12))
        elif i == 1:
            section_title = '<font size=12>' + line + '</font>'
            story.append(Paragraph(section_title, styles["h4"]))
            story.append(Spacer(1, 12))
        else:
            text = '<font size=12> %s </font>' % line
            story.append(Paragraph(text, styles['Normal']))

    story.append(Spacer(1, 12))

    story.append(Paragraph("", styles["Normal"]))
    story.append(Paragraph("<font size=12>Gene Information:</font>", styles["h4"]))
    story.append(Paragraph("", styles["Normal"]))

    pham_possible_starts = pham.total_possible_starts

    text_style = styles["Normal"]
    text_style.fontName = 'Helvetica'
    text_style.leading = 12

    if only_pham:  # The text if working on single pham, no particular phage
        gene_list = list(pham.genes.values())
        gene_list.sort(key=lambda x: x.phage_name)
        for gene in gene_list:
            candidate_starts = ""
            for start_num, start in zip(gene.alignment_candidate_start_nums, gene.alignment_candidate_starts):
                if start_num in gene.alignment_annot_start_nums:
                    count = gene.alignment_annot_counts_by_start[start_num]
                    message = 'MA: <b>' + str(count) + '</b>'
                    candidate_starts += '(Start: ' + str(start_num) + ' @' + str(gene.alignment_index_to_coord(start)) + \
                                        ' has ' + str(count) + " MA's), "
                else:
                    candidate_starts += '(' + str(start_num) + ', ' + str(gene.alignment_index_to_coord(start)) + '), '

            story.append(Paragraph(" Gene: %s \n Start: %s, Stop: %s, Start Num: %s " % (gene.full_name,
                                   gene.start_codon_location, gene.stop_codon_location,
                                   gene.suggested_start["current_start_number"]), text_style))

            story.append(Paragraph(" Candidate Starts for %s: " % gene.full_name, text_style))
            story.append(Paragraph("     " + candidate_starts, text_style))
            story.append(Spacer(1, 12))
    else:   # if working on a pham report for one particular phage then only list starts for that phage
        for gene in genes_in_phage:
            candidate_starts = ""
            for start_num, start in zip(gene.alignment_candidate_start_nums, gene.alignment_candidate_starts):
                if start_num in gene.alignment_annot_start_nums:
                    count = gene.alignment_annot_counts_by_start[start_num]
                    message = 'MA: <b>' + str(count) + '</b>'
                    candidate_starts += '(Start: ' + str(start_num) + ' @' + str(gene.alignment_index_to_coord(start)) + \
                                        ' has ' + str(count) + " MA's), "
                else:
                    candidate_starts += '(' + str(start_num) + ', ' + str(gene.alignment_index_to_coord(start)) + '), '

            story.append(Paragraph(" Gene: %s \n Start: %s, Stop: %s, Start Num: %s " % (gene.full_name,
                                   gene.start_codon_location, gene.stop_codon_location,
                                   gene.suggested_start["current_start_number"]), text_style))

            story.append(Paragraph(" Candidate Starts for %s: " % gene.full_name, text_style))
            story.append(Paragraph("     " + candidate_starts, text_style))
            story.append(Spacer(1, 12))
    doc.build(story)


def make_pham_genome(phage_genes, phage_name, length, file_path):
    file_name = os.path.join(file_path, '%sPhamsGraph.pdf' % phage_name)
    if check_file(file_name):
        return
    phams_in_genome = [g['pham_no'] for g in phage_genes.values()]
    pham_colors = get_pham_colors(phams_in_genome)
    gd_diagram = GenomeDiagram.Diagram(phage_name)
    gd_track = gd_diagram.new_track(1, name=phage_name, greytrack=1)
    gd_pham_set = gd_track.new_set()
    # print "making genome page"
    for gene_id in sorted(phage_genes.keys()):
        phage_gene = phage_genes[gene_id]
        pham_no = phage_gene["pham_no"]
        gene = phage_gene["gene"]
        # print pham_no, gene.gene_id
        if pham_no is None:
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
        gene_feature = SeqFeature(gene_location)
        gene_number = get_gene_number(gene.gene_id)
        # label the gene with the gene number
        gd_pham_set.add_feature(gene_feature, name=str(gene_number), label=True, label_size=6, label_angle=75)
        # label gene with pham color and name
        gd_pham_set.add_feature(gene_feature, color=pham_color, name=str(pham_no), label=True, label_position='middle')
    
    # print type(length), length
    gd_diagram.draw(format='linear', orientation='portrait', pagesize=reportlab.lib.pagesizes.letter, fragments=8,
                    start=0, end=length)
    gd_diagram.write(file_name, "PDF")


def make_suggested_starts(phage_genes, phage_name, file_path):
    """
        Creates a PDF page of the suggested starts of a phage
        Genes are list in order
        {Gene Name} is in Pham {Number}: {Suggested Start Coordinates}
    """
    file_name = os.path.join(file_path, "%sSuggestedStarts.pdf" % phage_name)
    text_file_name = os.path.join(file_path, "%sSuggestedStarts.txt" % phage_name)
    if check_file(file_name):
        return
    doc = SimpleDocTemplate(file_name, pagesize=reportlab.lib.pagesizes.letter)
    story = []
    just_text = []
    # print "making suggested starts page"
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name="paragraph"))
    styles.add(ParagraphStyle(name='Center', alignment=TA_CENTER))
    text = '<font size=14> Suggested Start Coordinates</font>'
    just_text.append(' Suggested Start Coordinates')
    story.append(Paragraph(text, styles['Center']))
    story.append(Spacer(1, 12))

    info = 'The suggested starts list below uses the original algorithm that considers both manual and computational ' \
           'annotations. This algorithm may or may not produce a suggested start. If there is insufficient ' \
           'consensus to suggest a start, all starts in the gene are listed by start number and position. '
    story.append(Paragraph(info, styles['paragraph']))
    story.append(Spacer(1, 12))

    # items to build table
    summary_data = list()

    flagged_startnum_rows = []


    headers = ["Gene", "Pham\nNum", "Pham\nsize", "Start\nNum","Start\nCoord", "Inform\nAnnots", "Agree vs.\ntop Alt"]
    summary_data.append(headers)

    for gene_id in sorted(phage_genes.keys()):
        phage_gene = phage_genes[gene_id]
        pham = phage_gene["pham_no"]
        gene = phage_gene["gene"]
        suggested_start = phage_gene["suggested_start"]
        if pham is None:
            text = '<font size=12> %s is not a member of an existing Pham </font>' % gene.full_name
            simple_text = '%s is not a member of an existing Pham ' % gene.full_name
        else:
            text = '<font size=12> %s is in Pham %s:  %s </font>' % (gene.full_name, pham, suggested_start)
            simple_text = '%s is in Pham %s:  %s ' % (gene.full_name, pham, suggested_start)
        story.append(Paragraph(text, styles['Normal']))
        just_text.append(simple_text)

        # building summary table colm 0 geneID
        gene = phage_genes[gene_id]['gene']
        if gene.pham_no is None:
            continue
        gene_summary = list()
        gene_summary.append(gene.full_name)

        # Colm 1 Pham num
        gene_summary.append(gene.pham_no)

        # Colm 2 Pham size
        gene_summary.append(gene.pham_size)

        # Colm 3 start number
        gene_summary.append(gene.alignment_start_num_called)

        # Colm 4 start coordinate
        start = gene.start_codon_location
        gene_summary.append(start)

        # Colm 5 Number of informative Manual annots
        num_informative = sum(gene.alignment_annot_start_counts)
        gene_summary.append(str(num_informative))

        # Colm 6 supporting annots vs best alternative

        if gene.alignment_start_num_called not in gene.alignment_annot_start_nums:
            support_num=0
            if len(gene.alignment_annot_start_counts) > 0:
                alternative = max(gene.alignment_annot_start_counts)
            else:
                alternative = 0
        else:
            support_num=gene.alignment_annot_counts_by_start[gene.alignment_start_num_called]

            alternative=0
            for start_num, count in gene.alignment_annot_counts_by_start.items():
                if start_num == gene.alignment_start_num_called:
                    continue
                else:
                    alternative = max(count, alternative)
        gene_summary.append(str(support_num) + ':' + str(alternative))

        summary_data.append(gene_summary)

        # Flag this row if the called start is one of the gene's bad adjacent starts
        if getattr(gene, "called_start_is_bad", False):
            flagged_startnum_rows.append(len(summary_data) - 1)


    story.append(Spacer(1, 12))

    text = 'The following table summarizes annotation results with a focus on only those manual ' \
           'annotations which match a start codon available in the listed gene. The total number of this type of ' \
           '"informative" annotation is given in the "Inform Annots" column. ' \
           'The "Agree vs. top Alt" column gives ' \
           'the total number of manual annotations that agree with the ' \
           'currently annotated start (as shown in the "Start Num" column) versus the number of manual annotations ' \
           'for the most manually annotated alternative start found in the listed gene. ' \
           'Background color is based on the difference in the two numbers ' \
           '(green if greater than +3, red if less than -3, and yellow if between +3 and -3).'

    story.append(Paragraph(text, styles['paragraph']))
    story.append(Spacer(1, 12))

    table = Table(summary_data)

    full_grid_style = [('INNERGRID', (0, 0), (-1, -1), 0.5, colors.black), ('BOX', (0, 0), (-1, -1), 0.5, colors.black)]
    align_styles = [('ALIGN', (0, 0), (-1, -1), 'CENTER')]
    color_styles = list()

    # colorize table for now just use absolute difference
    for row, gene_data in enumerate(summary_data[1:]):
        row += 1      # need to add 1 since skipping 1st row in for loop
        score_column = 6
        call_ratios = gene_data[score_column]

        counts = call_ratios.split(':')
        agree_count = int(counts[0])
        alter_count = int(counts[1])

        if alter_count - agree_count > 3:
            color_styles.append(('BACKGROUND', (score_column, row), (score_column, row), colors.red))
        elif agree_count - alter_count > 3:
            color_styles.append(('BACKGROUND', (score_column, row), (score_column, row), colors.green))
        else:
            color_styles.append(('BACKGROUND', (score_column, row), (score_column, row), colors.yellow))

    # Coloring in the "Start Num" column if any of the adjacent start flagged genes are present in the phage.
    start_num_column = 3  # columns: 0 Gene, 1 Pham Num, 2 Pham size, 3 Start Num, 4 Start Coord, 5 Inform Annots, 6 Agree vs top Alt
    for row in flagged_startnum_rows:
        color_styles.append(('BACKGROUND', (start_num_column, row), (start_num_column, row), colors.red))
        color_styles.append(('TEXTCOLOR', (start_num_column, row), (start_num_column, row), colors.white))


    table.setStyle(TableStyle(full_grid_style))
    table.setStyle(TableStyle(align_styles))
    table.setStyle(TableStyle(color_styles))
    story.append(table)

    doc.build(story)
    # print "writing text file"
    with open(text_file_name, 'w') as outfile:
        outfile.write("\n".join(just_text))


def make_fasta_file(genes, fasta_file):
    count = SeqIO.write(genes, fasta_file, 'fasta')
    # print "%s Fasta files written" % count


def main():
    args = parse_arguments()
    # print "making_files:main(); args.make is ", args.make
    if 'graph' in args.make:
        # print "making_files.main() make 'graph': args.pickle_file " + args.pickle_file
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
        # print "making_files.main(): Loading pickle file " + str(args.pickle_file)
        pham = pickle.load(open(args.pickle_file.strip('"'), 'rb'))
        graph_start_sites(args, pham, args.dir)
        # print "making_files.main(): 'text' phage is ", args.phage
        if not args.phage:
            # print "making_files.main() 'text': no phage"
            make_pham_text(args, pham, args.pham_no, args.dir, only_pham=True)
        else:
            make_pham_text(args, pham, args.pham_no, args.dir)

    if 'fasta' in args.make:
        pass
        # pickle_file = args.file
        # genes = pickle.load(open(args.pickle_file.strip('"'), 'rb'))
        # make_fasta_file(genes, (args.dir + '.fasta'))


if __name__ == "__main__":
    main()
