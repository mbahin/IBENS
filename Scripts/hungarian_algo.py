#!/usr/bin/python

__author__ = 'bahin'
''' Script to apply the Hungarian algorithm to decide the seats positions in the new offices. '''

import argparse
from munkres import Munkres
import sys
import random
import time

# Getting command-line options back
parser = argparse.ArgumentParser()
parser.add_argument('-i', dest='input', required=True, help='Input points matrix')
options = parser.parse_args()

# Reading the input matrix
people = []
points = []
with open(options.input, 'r') as input_file:
    # Indexing the seats line
    line = input_file.readline()
    seats = line.rstrip().split(',')[1:]
    # Indexing the points
    for line in input_file:
        people.append(line.split(',')[0])
        points.append(line.rstrip().split(',')[1:])

# Shuffling input rows
shuffled = range(len(people))
random.shuffle(shuffled)
random_people = []
random_points = []
for i in shuffled:
    random_people.append(people[i])
    random_points.append(points[i])

# Transforming the profit matrix into a cost one
cost_matrix = []
for row in random_points:
    cost_row = []
    for col in row:
        cost_row += [sys.maxsize - int(col)]
    cost_matrix += [cost_row]

# Computing the assignment
print '### Results ###'
m = Munkres()
indexes = m.compute(cost_matrix)
total = 0
for row, column in sorted(indexes, key=lambda tup: tup[1]):
    value = int(random_points[row][column])
    total += value
    raw_input('Next seat?')
    raw_input('Seat ' + seats[column] + ' is attributed to...')
    for l in random_people[row]:
        sys.stdout.write(l)
        sys.stdout.flush()
        # Suspense for people that have the same name first letter
        if (l == 'S') or (l == 'F') or (l == 'A'):
            time.sleep(2)
        time.sleep(0.3)
    print '\n'

# Showing data
print '### Data ###'
print 'Seats\t\t', '\t'.join(seats), '\n'
for i in range(len(random_people)):
    if len(random_people[i]) >= 8:
        print random_people[i] + '\t' + '\t'.join(random_points[i])
    else:
        print random_people[i] + '\t\t' + '\t'.join(random_points[i])