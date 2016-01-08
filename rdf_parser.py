# rdf_parser

import os, sys, rdflib

__usage__ = "rdf_parser.py <input_data_file> <Entrez2Gene.txt>"


def store_rdf(data_lines, genename_dict):

	from rdflib import Graph, Literal, BNode, Namespace, RDF, URIRef

	ppi = Namespace("http://ppi2rdf.org/proteins#")
	pmed = Namespace("http://www.ncbi.nlm.nih.gov/pubmed/")
	# Note: ppi.attribute vs. ppi[item]

	g = Graph()

	sys.stderr.write("Storing RDF")

	for i, line in enumerate(data_lines):

		if (i % 10000 == 0):
			sys.stderr.write(".")
		# Add classes
		interaction_id = line.split('\t')[10]
		pubmed_id = line.split('\t')[5]
		entrez_1 = line.split('\t')[0]
		entrez_2 = line.split('\t')[1]

		interaction = URIRef(ppi[interaction_id])
		reference = URIRef(pmed[pubmed_id])
		gene1 = URIRef(ppi[entrez_1])
		gene2 = URIRef(ppi[entrez_2])

		g.add( (interaction, RDF.type, ppi.interaction) )
		g.add( (reference, RDF.type, ppi.primaryRef) )
		g.add( (gene1, RDF.type, ppi.interactor) )
		g.add( (gene2, RDF.type, ppi.interactor) )

		# Bind to interaction
		method_name = line.split('\t')[2]
		method_id = line.split('\t')[3]
		source = line.split('\t')[9]
		interaction_type = line.split('\t')[8]

		g.add( (interaction, ppi.hasInteractor, gene1) )

		# if self-loop, add type
		if (interaction, ppi.hasInteractor, gene2) in g:
			g.add( (interaction, RDF.type, ppi.selfInteractor) ) 
		else:
			g.add( (interaction, ppi.hasInteractor, gene2) )

		g.add( (interaction, ppi.hasReference, reference) )

		g.add( (interaction, ppi.methodName, Literal(method_name)) )
		g.add( (interaction, ppi.methodID, Literal(method_id)) )
		g.add( (interaction, ppi.source, Literal(source)) )
		g.add( (interaction, ppi.interactionType, Literal(interaction_type)) )

		# Bind to reference
		author = line.split('\t')[4]

		g.add( (reference, ppi.firstAuthor, Literal(author)) )
		g.add( (reference, ppi.pubmed, Literal(pubmed_id)) )

		# Bind to interactor
		taxid_1 = line.split('\t')[6]
		taxid_2 = line.split('\t')[7]

		g.add( (gene1, ppi.taxId, Literal(taxid_1)) )
		g.add( (gene2, ppi.taxId, Literal(taxid_2)) )
		g.add( (gene1, ppi.entrez, Literal(entrez_1)) )
		g.add( (gene2, ppi.entrez, Literal(entrez_2)) )


		if entrez_1 in genename_dict:
			genename_1 = genename_dict[entrez_1]
			g.add( (gene1, ppi.geneName, Literal(genename_1)) )
		if entrez_2 in genename_dict:
			genename_2 = genename_dict[entrez_2]
			g.add( (gene2, ppi.geneName, Literal(genename_2)) )

	return g


def create_schema():

	sys.stderr.write("Creating PPISchema")
	from rdflib.namespace import Namespace, RDF, URIRef, RDFS, OWL
	from rdflib import Graph, Literal, BNode

	gs = Graph() # create new schema

	# Import Namespaces and Schema
	ppi = Namespace("http://ppi2rdf.org/proteins#")
	xmls = Namespace("http://www.w3.org/2001/XMLSchema#")

	# Classes 
	gs.add( (ppi.interaction, RDF.type, OWL.Class) ) 	# the interaction id
	gs.add( (ppi.interactor, RDF.type, OWL.Class) ) 	# the entrez id
	gs.add( (ppi.primaryRef, RDF.type, OWL.Class) ) 	# the pubmed id

	# Subclasses
	gs.add( (ppi.psdGene, RDF.type, OWL.Class) ) # post-synaptic density gene
	gs.add( (ppi.psdGene, RDFS.subClassOf, ppi.interactor) )

	gs.add( (ppi.selfInteractor, RDF.type, OWL.Class) )
	gs.add( (ppi.selfInteractor, RDFS.subClassOf, ppi.interaction) )

	# Interaction Properties
	gs.add( (ppi.hasInteractor, RDF.type, OWL.ObjectProperty) )
	gs.add( (ppi.hasInteractor, RDFS.domain, ppi.interaction) )
	gs.add( (ppi.hasInteractor, RDFS.range, ppi.interactor) )

	gs.add( (ppi.hasReference, RDF.type, OWL.ObjectProperty) )
	gs.add( (ppi.hasReference, RDFS.domain, ppi.interaction) )
	gs.add( (ppi.hasReference, RDFS.range, ppi.primaryRef) )	

	gs.add( (ppi.methodName, RDF.type, OWL.DatatypeProperty) )
	gs.add( (ppi.methodName, RDFS.domain, ppi.interaction) )
	gs.add( (ppi.methodName, RDFS.range, xmls.string) )

	gs.add( (ppi.methodId, RDF.type, OWL.DatatypeProperty) )
	gs.add( (ppi.methodId, RDFS.domain, ppi.interaction) )
	gs.add( (ppi.methodId, RDFS.range, xmls.string) )

	gs.add( (ppi.source, RDF.type, OWL.DatatypeProperty) )
	gs.add( (ppi.source, RDFS.domain, ppi.interaction) )
	gs.add( (ppi.source, RDFS.range, xmls.string) )

	gs.add( (ppi.interactionType, RDF.type, OWL.DatatypeProperty) )
	gs.add( (ppi.interactionType, RDFS.domain, ppi.interaction) )
	gs.add( (ppi.interactionType, RDFS.range, xmls.string) )

	# Interactor Properties

	gs.add( (ppi.entrez, RDF.type, OWL.DatatypeProperty) )
	gs.add( (ppi.entrez, RDFS.domain, ppi.interactor) )
	gs.add( (ppi.entrez, RDFS.range, xmls.string) )

	gs.add( (ppi.taxId, RDF.type, OWL.DatatypeProperty) )
	gs.add( (ppi.taxId, RDFS.domain, ppi.interactor) )
	gs.add( (ppi.taxId, RDFS.range, xmls.string) )	

	gs.add( (ppi.geneName, RDF.type, OWL.DatatypeProperty) )
	gs.add( (ppi.geneName, RDFS.domain, ppi.interactor) )
	gs.add( (ppi.geneName, RDFS.range, xmls.string) )	

	# homologene  (interactor) <--> (interactor)
	gs.add( (ppi.hasHomologene, RDF.type, OWL.ObjectProperty) )
	gs.add( (ppi.hasHomologene, RDFS.domain, ppi.interactor) )
	gs.add( (ppi.hasHomologene, RDFS.range, ppi.interactor) )

	# PrimaryRef Properties

	gs.add( (ppi.pubmed, RDF.type, OWL.DatatypeProperty) )
	gs.add( (ppi.pubmed, RDFS.domain, ppi.primaryRef) )
	gs.add( (ppi.pubmed, RDFS.range, xmls.string) )

	gs.add( (ppi.firstAuthor, RDF.type, OWL.DatatypeProperty) )
	gs.add( (ppi.firstAuthor, RDFS.domain, ppi.primaryRef) )
	gs.add( (ppi.firstAuthor, RDFS.range, xmls.string) )

	return gs

# Loads homologenes from file into 
def load_homologene(homologene_lines):

	HUMAN_ONLY = True # skip groups that don't contain human taxon

	from rdflib.namespace import Namespace, RDF, URIRef, RDFS, OWL
	from rdflib import Graph, Literal, BNode

	gh = Graph() # create new graph

	# Import Namespaces and Schema
	ppi = Namespace("http://ppi2rdf.org/proteins#")
	xmls = Namespace("http://www.w3.org/2001/XMLSchema#")

	homologous_genes = []
	current_id = 0

	for i, line in enumerate(homologene_lines):

		if i > 1000:
			break

		homolid = line.split('\t')[0]
		taxid = line.split('\t')[1]
		entrez = line.split('\t')[2]
		genename = line.split('\t')[3]

		if homolid == current_id:
			homologous_genes.append(entrez + '_' + taxid)
		else:

			# write to graph
			import itertools
			for gene_pair in itertools.permutations(homologous_genes, 2):

				taxid1 = gene_pair[0].split('_')[1]
				taxid2 = gene_pair[1].split('_')[1]


				if (taxid1 != '9606') and (taxid2 != '9606') and HUMAN_ONLY: # only look at humans
					continue
				else:
					gene1 = gene_pair[0].split('_')[0]
					gene2 = gene_pair[1].split('_')[0]
					gh.add( (URIRef(ppi[gene1]), ppi.hasHomologene, URIRef(ppi[gene2])) )

			homologous_genes = [entrez + '_' + taxid]
			current_id = homolid

	return gh


# Loads genenames into dictionary with key: entrezID
# and value: genename and returns it
def load_genenames(genename_lines):	

    genename_dict = {}
    for line in genename_lines:

    	entrezID = line.split('|')[0]
    	genename = line.split('|')[1].strip('\n')

    	genename_dict[entrezID] = genename

    return genename_dict

# Loads Biogrid PSI-MITAB data into lines of 
# "entrez 1 \t entrez 2 \t method name \t method id \t 
# first author \t pubmed id \t taxid 1 \t taxid2 \t 
# interaction type \t source \t interaction id \n"
def load_biogrid(lines):	

    sys.stderr.write("Loading biogrid data")
    new_lines = []
    for i, line in enumerate(lines):

    	if i == 0:
    		continue # Skip header
    	if (i % 10000 == 0):
    		sys.stderr.write(".")
        entrez_1 = line.split('\t')[0].strip("entrez gene/locuslink:")
        entrez_2 = line.split('\t')[1].strip("entrez gene/locuslink:")
        method_name = line.split('\t')[6].split('(')[1].strip(')')
        method_id = line.split('\t')[6].split('(')[0].strip("psi-mi:").strip('\"')
        author = line.split('\t')[7].strip('\"')
        pubmed_id = line.split('\t')[8].strip("pubmed:")
        taxid_1 = line.split('\t')[9].strip("taxid:")
        taxid_2 = line.split('\t')[10].strip("taxid:")
        interaction_type = line.split('\t')[11].split('(')[1].strip(')')
        source = "biogrid"
        interaction_id = line.split('\t')[13].strip("biogrid:")

        new_lines.append("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (entrez_1, entrez_2, method_name, method_id, author, pubmed_id, taxid_1, taxid_2, interaction_type, source, interaction_id))

    return new_lines

# Loads IntAct PSI-MITAB data into lines of 
# "entrez 1 \t entrez 2 \t method name \t method id \t 
# first author \t pubmed id \t taxid 1 \t taxid2 \t 
# interaction type \t source \t interaction id \n"
def load_intact(lines):	

	import fileinput
	sys.stderr.write("Loading intact data")

	#create the mapping dict from uniprot to entrez
	uniprot2entrez = dict()

	#this file was created by using the uniprot mapping tool using the unique uniprot ids from intact
	for line in fileinput.input('uniprot2entrez.txt'):
	    uniprot,entrez_id = line.split("\t")
	    uniprot2entrez[uniprot] = entrez_id

	new_lines = []
	for i, line in enumerate(lines):

		if i == 0:
			continue # Skip header
		if (i % 10000 == 0):
			sys.stderr.write(".")

		proteins = []
		if re.match('^[A-Z][0-9]',fields[0].strip('uniprotkb:')):
			uniprot_id = re.sub('-[0-9]','',fields[0].strip('uniprotkb:'))
			try:
				proteins.append(uniprot2entrez[uniprot_id].strip())
			except:
				pass
		if re.match('^[A-Z][0-9]',fields[1].strip('uniprotkb:')):
			uniprot_id = re.sub('-[0-9]','',fields[1].strip('uniprotkb:'))
			try:
				proteins.append(uniprot2entrez[uniprot_id].strip())
			except:
				pass
		if len(proteins)==2:
			entrez_1 = proteins[0]
			entrez_2 = proteins[1]

			method_name = line.split('\t')[6].split('(')[1].strip(')')
			method_id = line.split('\t')[6].split('(')[0].strip("psi-mi:").strip('\"')
			author = line.split('\t')[7].strip('\"')
			pubmed_id = line.split('\t')[8].strip("pubmed:")
			taxid_1 = line.split('\t')[9].strip("taxid:")
			taxid_2 = line.split('\t')[10].strip("taxid:")
			interaction_type = line.split('\t')[11].split('(')[1].strip(')')
			source = "intact"
			interaction_id = line.split('\t')[13].split('|')[0].strip("intact:")

			new_lines.append("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (entrez_1, entrez_2, method_name, method_id, author, pubmed_id, taxid_1, taxid_2, interaction_type, source, interaction_id))
		else:
			pass

	return new_lines

def load_psd(lines):

	from rdflib.namespace import Namespace, RDF, URIRef, RDFS, OWL
	from rdflib import Graph, Literal, BNode

	gpsd = Graph() # create new graph

	# Import Namespaces and Schema
	ppi = Namespace("http://ppi2rdf.org/proteins#")
	xmls = Namespace("http://www.w3.org/2001/XMLSchema#")

	gene_list = []
	for line in lines:
		gene = line.strip('\n')

		gpsd.add( (URIRef(ppi[gene]), RDF.type, ppi.psdGene) )

	return gpsd


def filter_genes(gene_file, data_lines):
	gene_list = []
	for gene in gene_file:
		gene_list.append(gene.strip('\n'))

	filtered_lines = []
	for line in data_lines:
		entrez_1 = line.split('\t')[0]
		entrez_2 = line.split('\t')[1]

		if (entrez_1 in gene_list) and (entrez_2 in gene_list):
			filtered_lines.append(line)

	exit(0)
	return filtered_lines

def parse_options():

	from rdflib import Graph, Literal, BNode
	import argparse

	parser = argparse.ArgumentParser(description='Provide data file from one or more sources in PSI-MI format. If no data files, only schema will output.')
	parser.add_argument('outfile', type=argparse.FileType('w'), help='provide output filename (i.e. out.nt)')
	parser.add_argument('--biogrid', metavar='DATAFILE', type=file)
	parser.add_argument('--intact', metavar='DATAFILE', type=file, help='uniprot2entrez.txt must exist in path')
	parser.add_argument('--genenames', metavar='Entrez2Gene.txt', type=file)
	parser.add_argument('--homologene', metavar='homologene.data', type=file)
	parser.add_argument('--psd', metavar='psd_genes.txt', type=file, help='adds psd triples to post-synaptic density entrez ids')
	parser.add_argument('--filter', metavar='genes.txt', type=file, help="filters data to entrez ids in newline delimited file")

	args = parser.parse_args(sys.argv[1:])

	data_lines = []

	if args.biogrid != None:
		data_lines.extend(load_biogrid(args.biogrid)) 
		print ".DONE"

	if args.intact != None:
		data_lines.extend(load_intact(args.intact)) 
		print ".DONE"

	if args.genenames != None:
		genename_dict = load_genenames(args.genenames)
	else:
		genename_dict = {}

	if args.homologene != None:
		gh = load_homologene(args.homologene)
	else:
		gh = Graph()

	if args.psd != None:
		gpsd = load_psd(args.psd)
	else:
		gpsd = Graph()

	if args.filter != None:
		data_lines = filter_genes(args.filter, data_lines)

	gs = create_schema()
	print ".DONE"
	if data_lines != []:
		g = store_rdf(data_lines, genename_dict)
		print ".DONE"

	graph = gs + g + gh + gpsd
	graph.serialize(destination=args.outfile, format='nt')

def main():

	parse_options()


if __name__ == "__main__":
    main()











