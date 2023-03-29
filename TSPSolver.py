#!/usr/bin/python3
import numpy as np

from which_pyqt import PYQT_VER
if PYQT_VER == 'PYQT5':
	from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
	from PyQt4.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT6':
	from PyQt6.QtCore import QLineF, QPointF
else:
	raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))


import time
from TSPClasses import *
import heapq
import itertools


class TSPSolver:
	def __init__( self, gui_view ):
		self._scenario = None
		self.pq = []

	def setupWithScenario( self, scenario ):
		self._scenario = scenario


	''' <summary>
		This is the entry point for the default solver
		which just finds a valid random tour.  Note this could be used to find your
		initial BSSF.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of solution,
		time spent to find solution, number of permutations tried during search, the
		solution found, and three null values for fields not used for this
		algorithm</returns>
	'''

	def defaultRandomTour( self, time_allowance=60.0 ):

		results = {}
		cities = self._scenario.getCities()
		ncities = len(cities)
		foundTour = False
		count = 0
		bssf = None
		start_time = time.time()
		while not foundTour and time.time()-start_time < time_allowance:
			# create a random permutation
			perm = np.random.permutation( ncities )
			route = []
			# Now build the route using the random permutation
			for i in range( ncities ):
				route.append( cities[ perm[i] ] )
			bssf = TSPSolution(route)
			count += 1
			if bssf.cost < np.inf:
				# Found a valid route
				foundTour = True
		end_time = time.time()
		results['cost'] = bssf.cost if foundTour else math.inf
		results['time'] = end_time - start_time
		results['count'] = count
		results['soln'] = bssf
		results['max'] = None
		results['total'] = None
		results['pruned'] = None
		return results

	def calculate_lower_bound(self):
		return

	def prune(self):
		return

	''' <summary>
		This is the entry point for the greedy solver, which you must implement for
		the group project (but it is probably a good idea to just do it for the branch-and
		bound project as a way to get your feet wet).  Note this could be used to find your
		initial BSSF.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution,
		time spent to find best solution, total number of solutions found, the best
		solution found, and three null values for fields not used for this
		algorithm</returns>
	'''

	def greedy( self,time_allowance=60.0 ):
		results = {}
		cities = self._scenario.getCities()
		# make a set of unvisited cities (using indices)
		unvisited = set(cities[1:])
		# start with first city
		current_city = cities[0]
		tour = [current_city]

		# run until all cities are visited
		start_time = time.time()
		while unvisited and time.time()-start_time < time_allowance:
			nearest, dist = self.nearest_neighbor(current_city, unvisited)
			# check if there is no complete path
			# if no complete path, break out of loop
			if nearest is None:
				break
			tour.append(nearest)
			unvisited.remove(nearest)
			current_city = nearest

		end_time = time.time()
		bssf = TSPSolution(tour)

		results['cost'] = bssf.cost if not unvisited else np.inf
		results['time'] = end_time - start_time
		results['count'] = 1 if not unvisited else 0
		results['soln'] = bssf
		results['max'] = None
		results['total'] = None
		results['pruned'] = None
		return results

	def nearest_neighbor(self, start_node, unvisited):
		nearest = None
		min_distance = np.inf
		for node in unvisited:
			dist = start_node.costTo(node)
			if dist < min_distance:
				nearest = node
				min_distance = dist
		return nearest, min_distance


	''' <summary>
		This is the entry point for the branch-and-bound algorithm that you will implement
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution,
		time spent to find best solution, total number solutions found during search (does
		not include the initial BSSF), the best solution found, and three more ints:
		max queue size, total number of states created, and number of pruned states.</returns>
	'''

	def reduce_matrix(self, matrix):
		cost_to_reduce = 0.0

		return None

	def initial_partial_path(self):
		pp = PartialPath()
		pp.level = 0
		cities = self._scenario.getCities()
		num_cities = len(cities)
		matrix = [[None for j in range(num_cities)] for i in range(num_cities)]

		for i in range(num_cities):
			for j in range(num_cities):
				if i == j:
					continue
				else:
					matrix[i][j] = cities[i].costTo(cities[j])

		print()

	def branchAndBound( self, time_allowance=60.0 ):
		# initialize priority queue with partial path consisting of only the starting node
		self.bssf = self.greedy()['cost']
		pq = []

		# generate initial partial path
		initial_pp = self.initial_partial_path()

		# while pq not empty
		start_time = time.time()
		while pq and time.time() - start_time < time_allowance:

			# dequeue partial path with lowest cost estimate from pq
			next_partial_path = min(pq, key=lambda x: x.key)

			# if at the end, if solution < BSSF, update BSSF
			# prune all partial paths with lower bound < BSSF

			# expand node into children nodes
			# calculate lower bounds for each child node
				# if lower bound for child node < BSSF, add to pq

		pass



	''' <summary>
		This is the entry point for the algorithm you'll write for your group project.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution,
		time spent to find best solution, total number of solutions found during search, the
		best solution found.  You may use the other three field however you like.
		algorithm</returns>
	'''

	def fancy( self,time_allowance=60.0 ):
		pass
