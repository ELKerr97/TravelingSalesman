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
import copy
from TSPClasses import *
import heapq
import itertools


class TSPSolver:
    def __init__(self, gui_view):
        self._scenario = None
        self.pq = []

    def setupWithScenario(self, scenario):
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

    def defaultRandomTour(self, time_allowance=60.0):

        results = {}
        cities = self._scenario.getCities()
        ncities = len(cities)
        foundTour = False
        count = 0
        bssf = None
        start_time = time.time()
        while not foundTour and time.time() - start_time < time_allowance:
            # create a random permutation
            perm = np.random.permutation(ncities)
            route = []
            # Now build the route using the random permutation
            for i in range(ncities):
                route.append(cities[perm[i]])
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

    def greedy(self, time_allowance=60.0):
        results = {}
        cities = self._scenario.getCities()
        # make a set of unvisited cities (using indices)
        unvisited = set(cities[1:])
        # start with first city
        current_city = cities[0]
        tour = [current_city]

        # run until all cities are visited
        start_time = time.time()
        while unvisited and time.time() - start_time < time_allowance:
            nearest, dist = self.nearest_neighbor(current_city, unvisited)
            # check if there is no complete path
            # if no complete path, break out of loop
            if nearest is None:
                break
            tour.append(nearest)
            unvisited.remove(nearest)
            current_city = nearest

        if tour[-1].costTo(tour[0]) == np.inf:
            unvisited.add(1)
            bssf = None
        else:
            bssf = TSPSolution(tour)

        end_time = time.time()

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

    def reduce_row(self, matrix, row):
        cost_to_reduce = 0.0
        reduced_row = row
        smallest = min(row)
        # if not 0, subtract it from each element in row
        if smallest != 0:
            row = [x - smallest for x in row]
            # add to cost_to_reduce
            cost_to_reduce += smallest
        return reduced_row, cost_to_reduce

    def reduce_matrix(self, matrix):
        cost_to_reduce = 0.0
        # reduce each row
        for i in range(len(matrix)):
            # find min in row
            row = matrix[i]
            smallest = min(row)
            # if not 0, subtract it from each element in row
            if smallest != 0:
                matrix[i] = [x - smallest for x in row]
                # add to cost_to_reduce
                cost_to_reduce += smallest

        # reduce each column
        col_index = 0
        for col in zip(*matrix):
            # find min in col
            # turning tuple into list
            col = list(col)
            smallest = min(col)
            # if not 0, subtract it from each element in col
            if smallest != 0:
                for i in range(len(matrix)):
                    matrix[i][col_index] -= smallest
                cost_to_reduce += smallest
            col_index += 1

        return matrix, cost_to_reduce

    def initial_partial_path(self):
        pp = PartialPath()
        pp.level = 0
        cities = self._scenario.getCities()
        num_cities = len(cities)
        matrix = [[np.inf for j in range(num_cities)] for i in range(num_cities)]

        for i in range(num_cities):
            for j in range(num_cities):
                if i == j:
                    continue
                else:
                    matrix[i][j] = cities[i].costTo(cities[j])

        initial_pp = PartialPath()
        initial_pp.level = 0
        initial_pp.reduced_cost_matrix, initial_pp.lower_bound = self.reduce_matrix(matrix)
        initial_pp.set_key()
        initial_pp.tour.append(0)
        initial_pp.target_cols = [i for i in range(1, len(matrix))]
        return initial_pp


    def branchAndBound(self, time_allowance=60.0):
        # initialize priority queue with partial path consisting of only the starting node
        results = {}
        cities = self._scenario.getCities()

        bssf = None
        # set the best solution cost
        bssf_cost = self.greedy()['cost']
        # number of branched pruned
        pruned = 0
        # number of times bssf was updated
        count = 0
        # number of children created
        states = 0
        # max length of pq
        max_pq = 0
        # initialize priority queue
        initial_pp = self.initial_partial_path()
        pq = [initial_pp]

        # start timer
        start_time = time.time()
        # while pq not empty
        while pq and time.time() - start_time < time_allowance:

            # dequeue partial path with lowest cost estimate from pq
            next_partial_path = min(pq, key=lambda x: x.key)
            pq.remove(next_partial_path)

            # check if at a leaf node.
            if next_partial_path.level == len(next_partial_path.reduced_cost_matrix) - 1:
                # first, check if the last city connects to the first.
                if next_partial_path.reduced_cost_matrix[next_partial_path.tour[-1]][0] == np.inf:
                    continue
                else:
                    next_partial_path.lower_bound += next_partial_path.reduced_cost_matrix[next_partial_path.tour[-1]][0]
                solution = TSPSolution([cities[i] for i in next_partial_path.tour])
                # if at the end and solution < BSSF, update BSSF
                if solution.cost < bssf_cost:
                    bssf = solution
                    bssf_cost = solution.cost
                    count += 1
                # prune all partial paths with lower bound < BSSF
                for pp in pq:
                    if pp.lower_bound > bssf_cost:
                        pq.remove(pp)
                        pruned += 1
                # have pruned function return how many were pruned. Store in results['pruned']

            else:
                # expand children
                # grab the row (the from-node)
                from_node = next_partial_path.get_from_node()
                # iterate through target columns for this particular partial path
                target_cols = next_partial_path.target_cols
                # create a new child for each target column (skip if target column is infinity)
                for target_col in target_cols:

                    # make a new child node (copy from parent)
                    child = copy.deepcopy(next_partial_path)
                    # if the item is infinity, skip (continue)
                    if next_partial_path.reduced_cost_matrix[from_node][target_col] == np.inf:
                        continue

                    # increase the child's depth by 1
                    child.level += 1
                    # add child to tour
                    child.tour.append(target_col)
                    # change target cols for child
                    child.target_cols.remove(target_col)
                    # increment states (number of times child is created)
                    states += 1

                    # grab cost at target, add to total and change to inf
                    target_cost = child.reduced_cost_matrix[from_node][target_col]
                    child.lower_bound += target_cost
                    child.reduced_cost_matrix[from_node][target_col] = np.inf

                    # infinity out the row
                    for i in range(len(child.reduced_cost_matrix)):
                        # if we are in an item not in target column
                        if i != target_col:
                            # if the item is 0
                            if child.reduced_cost_matrix[from_node][i] == 0:
                                # grab column to reduce
                                col_to_reduce = [child.reduced_cost_matrix[row][i] for row in
                                                 range(len(child.reduced_cost_matrix))]
                                # set the number in row to inf
                                col_to_reduce[from_node] = np.inf
                                # compute reduced column and its cost
                                reduced_col, cost = self.reduce_col(col_to_reduce)
                                # add cost to lower bound
                                child.lower_bound += cost
                                # replace reduced column in matrix
                                for row in range(len(reduced_col)):
                                    child.reduced_cost_matrix[row][i] = reduced_col[row]

                            else:
                                # just infinity out the item if not 0
                                child.reduced_cost_matrix[from_node][i] = np.inf

                    # infinity out the column
                    for i in range(len(child.reduced_cost_matrix)):
                        child.reduced_cost_matrix[i][target_col] = np.inf

                    # set the key based on lower bound and depth
                    child.set_key()
                    # if child's lowerbound < BSSF, add to the queue
                    if child.lower_bound < bssf_cost:
                        pq.append(child)
                        # update the max pq
                        if len(pq) > max_pq:
                            max_pq = len(pq)
                    else:
                        pruned += 1

        end_time = time.time()
        if bssf == None:
            results['cost'] = bssf_cost
        else :
            results['cost'] = bssf.cost

        if not pq:
            pruned += len(pq)


        results['time'] = end_time - start_time
        results['count'] = count
        results['soln'] = bssf
        results['max'] = max_pq
        results['total'] = states
        results['pruned'] = pruned
        return results

    def reduce_col(self, col):
        cost_to_reduce = 0.0
        smallest = min(col)
        # if not 0, subtract it from each element in col
        if smallest != 0:
            for i in range(len(col)):
                col[i] -= smallest
            cost_to_reduce += smallest
        return col, cost_to_reduce

    ''' <summary>
		This is the entry point for the algorithm you'll write for your group project.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution,
		time spent to find best solution, total number of solutions found during search, the
		best solution found.  You may use the other three field however you like.
		algorithm</returns>
	'''

    # Ant Colony Optimization
    def fancy(self, time_allowance=60.0):
        # Initialize constants
        # Influence of pheromone on route selection
        self.alpha = 1.0
        # Influence of distance on route selection
        self.beta = 5.0
        # Pheromone evaporation rate
        self.evaporation_rate = 0.02
        # Pheromone deposit amount
        self.q = 1.0
        # Number of ants used in each iteration
        num_ants = 100
        # Best solution found so far
        bssf = self.greedy()['soln']
        # Best distance found so far
        bssf_cost = bssf.cost
        # number of times bssf was updated
        count = 0
        # Initialize pheromone levels
        num_cities = len(self._scenario.getCities())
        self.pheromone_matrix = [[np.inf] * num_cities for _ in range(num_cities)]

        max_tau = 10000 / bssf_cost
        for i in range(num_cities):
            for j in range(num_cities):
                if i == j:
                    self.pheromone_matrix[i][j] = np.inf
                else:
                    self.pheromone_matrix[i][j] = random.uniform(0.1, max_tau)
                    self.pheromone_matrix[j][i] = self.pheromone_matrix[i][j]

        path_improved = True
        start_time = time.time()
        num_iterations = 100
        while num_iterations > 0 and time.time() - start_time < time_allowance:
            path_improved = False
            ants = []
            # For each ant
            for ant in range(num_ants):
                # Start each ant on a random city:
                cities = list(range(num_cities))
                ant_tour = [np.random.randint(num_cities)]
                cities.remove(ant_tour[0])
                current_city = ant_tour[-1]
                # While there are cities to visit...
                while cities:
                    best_city = self.best_neighbor(current_city, cities)
                    # Stop if ant cannot find best route
                    if best_city == None:
                        # and could not find complete path
                        break
                    ant_tour.append(best_city)
                    cities.remove(best_city)
                    current_city = best_city
                # Update best route if ant has found complete path
                if not cities:
                    ant_solution = TSPSolution([self._scenario.getCities()[i] for i in ant_tour])
                    ants.append(ant_tour)
                    if ant_solution.cost < bssf_cost:
                        bssf = ant_solution
                        bssf_cost = ant_solution.cost
                        path_improved = True
                        count += 1
            num_iterations -= 1

            # Reduce pheromone on all edges
            for i in range(len(self.pheromone_matrix)):
                for j in range(len(self.pheromone_matrix)):
                    self.pheromone_matrix[i][j] = self.pheromone_matrix[i][j]*(1-self.evaporation_rate)
                    self.pheromone_matrix[j][i] = self.pheromone_matrix[j][i]*(1-self.evaporation_rate)

            # Update pheromone matrix based on which ants visited which edge
            reward = 0.001
            for ant_tour in ants:
                cost = TSPSolution([self._scenario.getCities()[i] for i in ant_tour]).cost
                for i in range(len(ant_tour) - 1):
                    from_node = ant_tour[i]
                    to_node = ant_tour[i+1]
                    if cost == bssf_cost:
                        # reward the best path to find local optimum faster
                        self.pheromone_matrix[from_node][to_node] += (10 / cost) + reward
                        self.pheromone_matrix[to_node][from_node] += (10 / cost) + reward
                    else:
                        self.pheromone_matrix[from_node][to_node] += 10 / cost
                        self.pheromone_matrix[to_node][from_node] += 10 / cost


        # Return best route
        end_time = time.time()
        results = {}
        results['count'] = count
        results['cost'] = bssf_cost
        results['time'] = end_time - start_time
        results['soln'] = bssf
        return results

    def best_neighbor(self, from_index, unvisited):
        denominator = 0

        # calculate the denominator in probability function
        for to_index in unvisited:
            if to_index != from_index:
                from_city = self._scenario.getCities()[from_index]
                to_city = self._scenario.getCities()[to_index]
                pheromone_var = pow(self.pheromone_matrix[from_index][to_index],self.alpha)
                distance = from_city.costTo(to_city)
                if distance == np.inf or distance == 0:
                    break
                distance_var = pow(1.0/from_city.costTo(to_city),self.beta)
                if distance_var != np.inf:
                    denominator += pheromone_var * distance_var

        choices = []
        for to_index in unvisited:
            from_city = self._scenario.getCities()[from_index]
            to_city = self._scenario.getCities()[to_index]
            pheromone_var = pow(self.pheromone_matrix[from_index][to_index], self.alpha)
            distance = from_city.costTo(to_city)
            if distance == np.inf or distance == 0:
                break
            distance_var = pow(1.0/from_city.costTo(to_city), self.beta)
            if distance_var != np.inf and denominator != 0.0:
                numerator = pheromone_var*distance_var
                probability = numerator/denominator
                choices.append((to_index, probability))

        # choose a random city with a probability
        if choices:
            city = random.choices(population=[city[0] for city in choices], weights=[city[1] for city in choices], k=1)[0]
            return city
        else:
            return None
