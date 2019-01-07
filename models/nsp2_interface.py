import json
from pprint import pprint

class NSP2():
	def __init__(self):

		f = open('new_fasta.nsp2.json')

		self.data = json.load(f)
		"""
		print(len(self.data))
		print(self.data[0]['desc'])
		print(self.data[0].keys())
		print("bluebleubleh")
		print(self.data[1]['desc'])
		print(len(self.data))
		"""

		self.id_to_index = {}
		for i in range(len(self.data)):
			self.id_to_index[self.data[i]['desc']] = i

	def get_summed_attributes_vector(self, uniprot_id, start, end):
		group = self.data[self.id_to_index[uniprot_id]]
		list_attributes = ['rsa', 'phi', 'psi', 'q3', 'interface', 'asa', 'q8', 'disorder']
		#not included: 'q3_prob', 'q8_prob'
		total_attributes = []
		for attribute in list_attributes:
			vector = group[attribute][start:end]
			if attribute == "q3" or attribute == "q8":
				new_vector = [vector.count("H"), vector.count("S"), vector.count("C"), vector.count("T")]
				for num in new_vector:
					total_attributes.append(num)
			else:
				number = 0
				for thing in vector:
					number = number + thing
				total_attributes.append(number)
		

		return(total_attributes)

		

if __name__ == '__main__':
	a = NSP2().get_summed_attributes_vector("O11312", 5,20)
	print(a)

"""

RSA: relative solvent accessibility
ASA: Absolute solvent accessibility
SS3: 3-class secondary structure
SS8: 8-class secondary structure 
Phi: backbone dihedral phi angle (DSSP software)
Psi: backbone dihedral psi angle (DSSP software)
"""