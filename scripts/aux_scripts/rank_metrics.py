#! /usr/bin/env python
import openpyxl
import argparse
import numpy as np

########################### FUNCTIONS #######################
#############################################################

def open_rank_file(rank_file):
	known_ranks = []
	with open(rank_file,"r") as f:
		for line in f:
			line = line.rstrip().split("\t")
			group_name = line[-1]
			candidate_gene = line[0]
			score = line[1]
			percentage_score = line[2]
			known_ranks.append([candidate_gene, score, percentage_score, group_name])
	return known_ranks

def report_stats(data):
  report_stats = []
  report_stats.append(['Elements', len(data)])
  report_stats.append(['Max', round(np.max(data), 5)])
  report_stats.append(['Min', round(np.min(data), 5)])
  report_stats.append(['Average', round(np.mean(data),5)])
  report_stats.append(['Standard_Deviation', round(np.std(data),5)])
  report_stats.append(['Q1', round(np.quantile(data, 0.25),5)])
  report_stats.append(['Median', round(np.quantile(data, 0.5),5)])
  report_stats.append(['Q3', round(np.quantile(data, 0.75),5)])
  return report_stats

def report_ranks(gene_pos_ranks):
	report_ranks = []
	for gene_pos_rank in gene_pos_ranks:
		report_ranks.append(gene_pos_rank)
	return report_ranks

def get_cdf_values(known_ranks):
	known_ranks.sort(key= lambda x: float(x[2]))
	number_known_ranks = len(known_ranks)
	for i, known_rank in enumerate(known_ranks):
			known_rank.insert(3,(i+1)/number_known_ranks)
	return known_ranks

def get_hash2groups(rankings, by_column= 3):
	group2rankings = {}

	for row in rankings:
		group_name = row[by_column]
		if group2rankings.get(group_name) is None:
			group2rankings[group_name] = [row]
		else:
			group2rankings[group_name].append(row)
		
	return group2rankings


########################### OPTPARSE ########################
#############################################################

parser = argparse.ArgumentParser(description='extract metrics from rankings')
parser.add_argument("-r", "--rankings", dest= "rankings", help="The roots to the rankings file")
parser.add_argument("-e", "--execution_mode", dest= "execution_mode", default="stats", help="The mode of execution: stats or ranks")
parser.add_argument("-c", "--by_column", dest= "by_column", help="The column with the factors to separate groups", type= lambda x: int(x))
options = parser.parse_args()

########################### MAIN ############################
#############################################################
all_ranks = open_rank_file(options.rankings)

if options.by_column is not None:
  group2rankings = get_hash2groups(all_ranks, options.by_column)
else:
	group2rankings = {"all_groups": all_ranks}


for group_name, rankings in group2rankings.items():
	known_ranks = get_cdf_values(rankings)

	if known_ranks:
		if options.execution_mode == "stats":
			all_ranks = [float(rank_row[2]) for rank_row in known_ranks]
			for stat in report_stats(all_ranks):
				print(group_name + "\t" + "\t".join([str(col) for col in stat]))

		elif options.execution_mode == "ranks":
			for rank in report_ranks(known_ranks):
				print("\t".join([str(col) for col in rank]))
