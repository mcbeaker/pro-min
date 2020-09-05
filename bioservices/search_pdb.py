#!/usr/bin/env python

# from bioservices import PDB
import requests
if __name__ == '__main__':
    url = 'http://www.rcsb.org/pdb/rest/search'

query = """<?xml version="1.0" encoding="UTF-8"?>
    <orgPdbQuery>
    <version>B0907</version>
    <queryType>org.pdb.query.simple.ResolutionQuery</queryType>
    <description>ResolutionQuery: refine.ls_d_res_high.comparator=between refine.ls_d_res_high.min=0 refine.ls_d_res_high.max=1.5 </description>
    <refine.ls_d_res_high.comparator>between</refine.ls_d_res_high.comparator>
    <refine.ls_d_res_high.min>0</refine.ls_d_res_high.min>
    <refine.ls_d_res_high.max>1.5</refine.ls_d_res_high.max>
    </orgPdbQuery>
    """
print("Query: %s" % query)
print("Querying RCSB PDB REST API...")

header = {'Content-Type': 'application/x-www-form-urlencoded'}

response = requests.post(url, data=query, headers=header)

if response.status_code == 200:
    print("Found %d PDB entries matching query." % len(response.text))
    print("Matches: \n%s" % response.text)
else:
    print("Failed to retrieve results")
