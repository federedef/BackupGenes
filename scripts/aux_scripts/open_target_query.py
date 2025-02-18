#!/usr/bin/env python
import requests
import json
import argparse

query = """query KnownDrugsQuery(
  $efoId: String!
  $cursor: String
  $freeTextQuery: String
  $size: Int = 9000
) {
  disease(efoId: $efoId) {
    id
    knownDrugs(cursor: $cursor, freeTextQuery: $freeTextQuery, size: $size) {
      rows {
        
        target {
          id
          approvedName
          approvedSymbol
        }
      }
    }
  }
}"""

########################### OPTPARSE ########################
#############################################################

parser = argparse.ArgumentParser(description='extract metrics from rankings')
parser.add_argument("-i", "--input_file", dest= "input_file", help="")
options = parser.parse_args()

##############################################################
##############################################################

def get_diseases(file):
	diseases = {}
	with open(file, "r") as f:
		for line in f:
			fields = line.strip().split("\t")
			print(fields)
			diseases[fields[0]] = fields[1]
	return diseases

diseases = get_diseases(options.input_file)
base_url = "https://api.platform.opentargets.org/api/v4/graphql"

disease2target={}
for disease, efo in diseases.items():
  r = requests.post(base_url, json={"query": query, "variables": {"efoId": efo}})
  api_response = json.loads(r.text)["data"]["disease"]
  disease2target[disease] = None
  if api_response is not None:
    disease2target[disease] = []
    for target in api_response["knownDrugs"]["rows"]:
      disease2target[disease].append(target["target"]["approvedSymbol"])
      disease2target[disease] = list(set(disease2target[disease]))

with open("disease2target", "w") as f:
  for disease, targets in disease2target.items():
    if targets: f.write(f"{disease}\t{','.join(targets)}\n")

