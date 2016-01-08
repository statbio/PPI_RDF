# PPI_RDF

Provide data file from one or more sources in PSI-MI format. If no data files, only schema will output.

positional arguments:
  outfile               					provide output filename (i.e. out.nt)

optional arguments:
  --help, -h            				show this help message and exit
  --biogrid DATAFILE				import biogrid MITAB data
  --intact DATAFILE     				uniprot2entrez.txt must exist in path
  --genenames Entrez2Gene.txt		import genenames from entrez ids
  --homologene homologene.data		add homologene relation between interactors
  --psd psd_genes.txt   			adds psd triples to postsynaptic density entrez ids
  --filter genes.txt   				filters data to entrez ids in newline delimited file
  
  
