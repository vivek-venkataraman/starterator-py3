# Copyright (c) 2013, 2014 All Right Reserved, Hatfull Lab, University of Pittsburgh
#
# THIS CODE AND INFORMATION ARE PROVIDED "AS IS" WITHOUT WARRANTY OF ANY
# KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A
# PARTICULAR PURPOSE.  USE AT YOUR OWN RISK.
#
# Marissa Pacey
# April 4, 2014
# Annotation of Phage genomes using Glimmer and Genemark

import requests
from html.parser import HTMLParser
from bs4 import BeautifulSoup
from operator import attrgetter

GLIMMER_URL = "http://www.ncbi.nlm.nih.gov/genomes/MICROBES/glimmer_3.cgi"
GENEMARK_URL = "http://www.ncbi.nlm.nih.gov/genomes/MICROBES/genemark.cgi"

class PredictedGene(object):
	def __init__(self, orf_id, start, end, orientation, score=None, **kwargs):
		self.id = orf_id
		self.start = start # 5' end
		self.stop = end  # 3' end
		self.orientation = orientation
		self.score = score
		self.other = kwargs

	def __str__(self):
		string = "Id: %s\tStart: %s\tStop:%s\tOrientation:%s " %  (self.id, self.start, self.stop, self.orientation)
		return "Id: %s\tStart: %s\tStop:%s\tOrientation:%s "% (self.id, self.start, self.stop, self.orientation)

class GenePredictions(object):
	def __init__(self, name, gene_list):
		self.name = name
		self.genes = gene_list

	def combine(self, other):
		"""

		"""
### 	TODO:
#			does not work if two genes are different orientations
#			need to think more about it still
		self_index = 0
		other_index = 0
		new_index = 1
		new_prediction = GenePredictions("Combine %s with %s" % (self.name, other.name), [])
		while(self_index < len(self.genes) or other_index < len(other.genes)):
			predict_gene = self.genes[self_index]
			other_gene = other.genes[other_index]
			print(new_index)
			print("glimmer ", predict_gene)
			print("genemark", other_gene)
			# if predict_gene.orientation != other_gene.orientation:
			# 	new_predictions.genes.append(predict_gene)
			# 	new_predictions.genes.append(other_gene)
			# 	other_gene += 1
			# if predict_gene.orientation = other_gene.orientation:
			if predict_gene.stop == other_gene.stop:
				predict_gene.other[self.name] = predict_gene.start
				predict_gene.other[other.name] = other_gene.start
				predict_gene.id = new_index
				new_prediction.genes.append(predict_gene)
				other_index += 1
				self_index += 1
				new_index += 1
			elif predict_gene.stop > other_gene.stop:
				other_gene.id = new_index
				other_gene.other[self.name] = None
				other_gene.other[other.name] = other_gene.start
				new_prediction.genes.append(other_gene)
				other_index +=1
				new_index +=1
			else: # predict_gene.stop < other_gene.stop:
				predict_gene.id = new_index
				new_prediction.genes.append(predict_gene)
				self_index += 1
				new_index += 1
			
			# else: #the orientations of genes are different
			# 	if predict_gene.orientation == 'F': #other gene is R then

		return new_prediction

	def __str__(self):
		string = "Name: %s\nGenes:\n" % self.name
		for gene in self.genes:
			string += "%s\n" % gene

		return string

	def order_prediction(self):
		new_index = 1
		for gene in self.genes:
			if gene.orientation == 'R':
				start = gene.stop
				stop = gene.start
				gene.start = start
				gene.stop = stop
		self.genes = sorted(self.genes, key=attrgetter('stop'))
		for gene in self.genes:
			gene.id = new_index
			new_index += 1


def genemark_annotate(text):
	genemark_list = []
	index = 1
	lines = text.splitlines()
	predicted_genes = 0
	for line in lines:
		if predicted_genes == 0 and "Predicted genes" in line:
			predicted_genes = 1
			continue
		if predicted_genes == 0:
			continue
		if predicted_genes < 3:
			predicted_genes += 1
			continue
		if "======" in line:
			break
		properties = line.split()
		orientation = properties[1]
		left_end = properties[2]
		right_end = properties[3]
		if orientation == '+':
			gene = PredictedGene(index, int(left_end), int(right_end), 'F')
		else:
			gene = PredictedGene(index, int(right_end), int(left_end), 'R')
		genemark_list.append(gene)
		index += 1
	genemark_list = sorted(genemark_list, key=attrgetter('stop'))
	genemark_prediction = GenePredictions("Genemark", genemark_list)
	return genemark_prediction

def glimmer_annotate(text):
	glimmer_list = []
	index = 1
	lines = text.splitlines()
	for line in lines:
		print(line, index)
		if "orf" in line and "orfID" not in line:
			properties = line.split()
			start = properties[1]
			stop = properties[2]
			orientation = -1 if int(properties[3]) < 1 else 1
			if orientation < 0:
				gene =  PredictedGene(index, int(start), int(stop), 'R')
			else:
				gene = PredictedGene(index, int(start), int(stop), 'F')
			print(gene)
			glimmer_list.append(gene)
			index += 1

	glimmer_list = sorted(glimmer_list, key=attrgetter('stop'))
	glimmer_prediction = GenePredictions("Glimmer", glimmer_list)
	return glimmer_prediction

def post_to_ncbi(end_url, fasta_file):
	sequence = ""
	with open(fasta_file, 'rU') as fasta:
		for line in fasta:
			sequence += line
	body = { "sequence": sequence, "gencode": 11, "topology": 0}
	r = requests.post(end_url, data=body);
	return r.text

def parse_ncbi_response(response):
	soup = BeautifulSoup(response)
	pre_tag = soup.pre #currently, the algorithm's outputs are in the only <ore> tag of the response HTML
						# may change, but for now, it works.
	return pre_tag.contents[0]

def auto_annotate(fasta_file):
	glimmer_response = post_to_ncbi(GLIMMER_URL, fasta_file)
	genemark_response = post_to_ncbi(GENEMARK_URL, fasta_file)
	glimmer = parse_ncbi_response(glimmer_response)
	genemark = parse_ncbi_response(genemark_response)

	glimmer_predict = glimmer_annotate(glimmer)
	print(glimmer_predict)
	genemark_predict = genemark_annotate(genemark)
	print(genemark_predict)
	# print glimmer
	combined = glimmer_predict.combine(genemark_predict)
	combined.order_prediction()
	return combined
