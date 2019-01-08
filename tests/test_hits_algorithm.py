import unittest
from models import hits_algorithm
import os
import pandas as pd

class QuickTest(unittest.TestCase):
	def check_correct_hits(self, initial_list, the_hits, predicted_hits):
		for num in range(len(initial_list)):
			if num in the_hits:
				self.assertTrue(predicted_hits[num])
			else:
				self.assertFalse(predicted_hits[num])

	def test_benchmark(self):
		[99, 3, 0, 0, 0, 0, 1, 25, 12, 0, 0, 53, 44, 10, 127, 109, 251, 56, 35, 58, 36, 37, 36, 14, 16, 308, 319, 207, 32, 73, 35, 23, 42, 32, 7]
		[0, 0, 0, 0, 32, 56, 17, 10]
		[0, 14, 0, 0, 103, 59]


		list_hits_a = [124, 220, 101, 7, 0, 0, 0, 7, 0, 0, 0, 1, 6, 1000, 0, 0, 0, 0, 0, 1, 18, 395, 256, 0, 0, 0, 0, 0, 0, 0, 10, 14, 0, 0, 0, 0, 10, 0, 0, 0, 0, 0, 0, 0, 1, 48, 6, 3, 0, 1, 0, 10, 51, 6, 0, 0, 1, 0, 1, 0, 0, 2, 37, 109, 95, 26, 21, 32, 5, 6, 0, 0, 0, 0, 1, 1, 1, 1]
		the_hits = [0, 1, 2, 21, 22, 63, 64]
		predicted_hits = hits_algorithm.find_hits(list_hits_a)
		self.check_correct_hits(list_hits_a, the_hits, predicted_hits)

		list_hits_b = [0, 14, 0, 0, 103, 59]
		the_hits = [4,5]
		predicted_hits = hits_algorithm.find_hits(list_hits_b)
		self.check_correct_hits(list_hits_b, the_hits, predicted_hits)

		list_hits_b = [0, 14, 0, 0, 3, 2]
		predicted_hits = hits_algorithm.find_hits(list_hits_b)
		self.assertEqual(len(predicted_hits), 0)

		list_hits_b = [0, 140,]
		predicted_hits = hits_algorithm.find_hits(list_hits_b)
		self.assertEqual(len(predicted_hits), 0)