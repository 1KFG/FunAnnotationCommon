#!/usr/bin/env python3

import os
import csv
import xml.etree.ElementTree as ET
from Bio import Entrez
Entrez.email = 'jason.stajich@ucr.edu'
from joblib import Memory

cachedir = os.path.join(os.getcwd(), 'cache')

if not os.path.exists(cachedir):
    os.makedirs(cachedir)
memory = Memory(cachedir, verbose=0)

@memory.cache
def get_bioproject_prefix(BIOPROJECTID):
    LOCUSTAG = ""
    for project in BIOPROJECTID.split(';'):
        # remove umbrella projects
        print(f'Checking {project}')
        if project == "PRJNA533106" or project == "PRJEB40665" or project == "PRJEB43510":
            continue
        try:
            bioproject_handle = Entrez.efetch(db="bioproject",id = project)
            projtree = ET.parse(bioproject_handle)
            projroot = projtree.getroot()
        
            lt = projroot.iter('LocusTagPrefix')
            for locus in lt:
                LOCUSTAG = locus.text
        except Exception as e:
            print(f'Error {e}')
        if len(LOCUSTAG) > 0:
            break
    return LOCUSTAG

accessions = 'ncbi_accessions.csv'
accession_taxonomy = 'ncbi_accessions_taxonomy.csv'

outsamples = 'samples.csv'

accession_dict = {}
fields = ['SPECIES', 'STRAIN', 'BIOPROJECT', 'NCBI_TAXONID', 'BUSCO_LINEAGE']
with open(accessions, 'r') as f:
    reader = csv.reader(f)
    header = next(reader)
    hdr2dict = {header[i]: i for i in range(len(header))}
    for row in reader:
        acc = row[hdr2dict['ACCESSION']]
        species = row[hdr2dict['SPECIES']]
        strain = row[hdr2dict['STRAIN']]
        asm_name = row[hdr2dict['ASM_NAME']]
        asm_base = f'{acc}_{asm_name}'
        accession_dict[asm_base] = [species, strain, 
                                    row[hdr2dict['BIOPROJECT']], 
                                    row[hdr2dict['NCBI_TAXID']] ]
        
with open(accession_taxonomy, 'r') as f:
    reader = csv.reader(f)
    header = next(reader)
    fields.extend(header[4:])
    for row in reader:
        asm_base = row[0]
        
        # add the BUSCO lineage to the dictionary
        buscodb = 'fungi'
        if row[4] == "Ascomycota" or row[4] == "Basidiomycota":
            buscodb = 'dikarya'
        elif row[4] == "Mucoromycota":
            buscodb = 'mucoromycota'
        accession_dict[asm_base].append(buscodb)
        
        # add Taxonomy info string to dictionary
        if asm_base in accession_dict:            
            accession_dict[asm_base].extend(row[4:])
        else:
            print(f'{asm_base} not found in accession_dict')
            
            
print(f'{len(accession_dict)} assemblies processed')

# add LOCUS tag field
fields.append('LOCUSTAG')

# add LOCUSTAG to the dictionary
n = 1
for asm in accession_dict:
    bioproject = accession_dict[asm][2]
    
    locus = get_bioproject_prefix(bioproject)
    print(f'{bioproject} -> {locus}')
    if len(locus) == 0:
        locus = f'FUN{n:03}'   # use strain as prefix
        n += 1
    locus = f'f{locus}'
    accession_dict[asm].append(locus)
    
with open(outsamples, 'wt') as f:
    writer = csv.writer(f)
    writer.writerow(fields)
    for asm in accession_dict:
        writer.writerow([asm] + accession_dict[asm])

