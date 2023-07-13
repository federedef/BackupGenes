#!/usr/bin/env python
import argparse

########################### FUNCTIONS #######################
#############################################################

def write_file(output_name, data):
	with open(output_name, "w") as f:
		for row in data:
			f.write("\t".join(str(el) for el in row) + "\n")

def add_tags(file2tag, tags, group_column=0, cases_column=1):
	tagged_file=[]
	for row in file2tag: 
		group = str(row[group_column])
		case_for_group = row[cases_column]
		if tags.get(group) is not None:
			if tags[group][case_for_group] is not None:
				row.append(tags[group][case_for_group])
			tagged_file.append(row)
	return tagged_file

def load_file(file):
	parsed_file = []
	with open(file) as f:
		for line in f:
			fields = line.strip().split("\t")
			parsed_file.append(fields)
	return parsed_file

def tagfile2hash(file):
	tags = {}
	for row in file:
		group = row[0]
		case_for_group = row[1]
		tag = row[2]
		if tags.get(group) is None:
			tags[group] = {}
		tags[group][case_for_group] = int(tag)
	return tags



########################### OPTPARSE ########################
#############################################################

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", default=None, help="The root to the file to add the tags")
parser.add_argument("-g", "--group_column", type=int, default=0, help="The index to column of groups")
parser.add_argument("-c", "--cases_column", type=int, default=1, help="The index to the column of cases")
parser.add_argument("-t", "--tag_file", default=None, help="The root to the file to add the tags")
parser.add_argument("-o", "--output_name", default="added_tags", help="The name of the ranked file")

options = parser.parse_args()

######################### MAIN #########################
########################################################


file2tag = load_file(options.input_file)
tags_file = load_file(options.tag_file)
tags = tagfile2hash(tags_file)
file_with_tags = add_tags(file2tag, tags, options.group_column, options.cases_column)


if file_with_tags is not None:
	write_file(options.output_name,file_with_tags)
