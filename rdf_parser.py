# rdf_parser

import os, sys, rdflib

__usage__ = "rdf_parser.py <input_data_file> <Entrez2Gene.txt>"


def store_rdf(data_lines, genename_dict):

	from rdflib import Graph, Literal, BNode, Namespace, RDF, URIRef

	ppi = Namespace("http://sample.org/proteins#")
	pmed = Namespace("http://www.ncbi.nlm.nih.gov/pubmed/")
	# Note: ppi.attribute vs. ppi[item]

	g = Graph()

	print "STORING RDF"

	for i, line in enumerate(data_lines):

		if (i % 10000 == 0):
			print i
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
		g.add( (interaction, ppi.hasInteractor, gene2) )
		g.add( (interaction, ppi.hasReference, reference) )

		g.add( (interaction, ppi.methodName, Literal(method_name)) )
		g.add( (interaction, ppi.methodID, Literal(method_id)) )
		g.add( (interaction, ppi.source, Literal(source)) )
		g.add( (interaction, ppi.interactionType, Literal(interaction_type)) )

		# Bind to reference
		author = line.split('\t')[4]

		g.add( (reference, ppi.firstAuthor, Literal(author)) )

		# Bind to interactor
		taxid_1 = line.split('\t')[6]
		taxid_2 = line.split('\t')[7]

		g.add( (gene1, ppi.taxId, Literal(taxid_1)) )
		g.add( (gene2, ppi.taxId, Literal(taxid_2)) )

		if entrez_1 in genename_dict:
			genename_1 = genename_dict[entrez_1]
			g.add( (gene2, ppi.geneName, Literal(genename_1)) )
		if entrez_2 in genename_dict:
			genename_2 = genename_dict[entrez_2]
			g.add( (gene2, ppi.geneName, Literal(genename_2)) )


	return g


def create_schema():

	from rdflib.namespace import Namespace, RDF, URIRef, RDFS, OWL
	from rdflib import Graph, Literal, BNode

	gs = Graph() # create new schema

	# Import Namespaces and Schema
	ppi = Namespace("http://sample.org/proteins#")
	xmls = Namespace("http://www.w3.org/2001/XMLSchema#")

	# Classes 
	gs.add( (ppi.interaction, RDF.type, OWL.Class) ) 	# the interaction id
	gs.add( (ppi.interactor, RDF.type, OWL.Class) ) 	# the entrez id
	gs.add( (ppi.primaryRef, RDF.type, OWL.Class) ) 	# the pubmed id

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

	gs.add( (ppi.taxId, RDF.type, OWL.DatatypeProperty) )
	gs.add( (ppi.taxId, RDFS.domain, ppi.interactor) )
	gs.add( (ppi.taxId, RDFS.range, xmls.string) )	

	gs.add( (ppi.geneName, RDF.type, OWL.DatatypeProperty) )
	gs.add( (ppi.geneName, RDFS.domain, ppi.interactor) )
	gs.add( (ppi.geneName, RDFS.range, xmls.string) )	

	# PrimaryRef Properties

	gs.add( (ppi.firstAuthor, RDF.type, OWL.DatatypeProperty) )
	gs.add( (ppi.firstAuthor, RDFS.domain, ppi.primaryRef) )
	gs.add( (ppi.firstAuthor, RDFS.range, xmls.string) )

	return gs
# Loads genenames into dictionary with key: entrezID
# and value: genename and returns it
def load_genenames(genename_file):	
    f = open(genename_file)
    lines = f.readlines()
    f.close()

    genename_dict = {}
    for line in lines:

    	entrezID = line.split('|')[0]
    	genename = line.split('|')[1].strip('\n')

    	genename_dict[entrezID] = genename

    return genename_dict

# Loads Biogrid PSI-MITAB data into lines of 
# "entrez 1 \t entrez 2 \t method name \t method id \t 
# first author \t pubmed id \t taxid 1 \t taxid2 \t 
# interaction type \t source \t interaction id \n"
def load_biogrid(biogrid_file):	
    f = open(biogrid_file)
    f.readline() # skip header  
    lines = f.readlines()
    f.close()

    new_lines = []
    for i, line in enumerate(lines):

    	if (i % 10000 == 0):
    		print i
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

def main():
    if len(sys.argv) == 3 and (sys.argv[2]) == 'Entrez2Gene.txt':
        data_lines = load_biogrid(sys.argv[1])
        genename_dict = load_genenames(sys.argv[2])
        g = store_rdf(data_lines, genename_dict)
        gs = create_schema()

    	output_filename = (sys.argv[1]).split('.')[0] + ".nt"
        output = open(output_filename,'w')
        output.write(gs.serialize(format='nt'))
        output.write(g.serialize(format='nt'))
        output.close()

    else:
        print __usage__
    


if __name__ == "__main__":
    main()

