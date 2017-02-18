import numpy as np
import os

text_file = open("genefile.fasta", "w")

for root, dirs, files in os.walk('Gadiformes'):
	for file in files:
		with open(os.path.join(root, file), "r") as auto:
			seq=auto.readlines()
			#seq=auto.readlines()[1:]
			#print seq
			for line in seq:
				#print line
				text_file.write("%s" % line)
			text_file.write('\n')

text_file.close()