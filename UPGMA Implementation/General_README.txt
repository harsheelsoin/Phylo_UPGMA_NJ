#To generate genefile.fasta
Usage:
  python makeGeneFile.py

#To generate 'genealign.fasta'
Usage: (Clustal Omega)
clustalo --infile='genefile.fasta' --infmt=fa --outfile='genealign.fasta' --outfmt=fa -seqtype=DNA --iterations=5 --max-guidetree-iterations=3 --max-hmm-iterations=3 --full-iter --full --dealign --force

#To generate Distance Matrix and store in text file 'distanceMatrix.txt'
Usage:
  python makeDistanceMatrix.py




