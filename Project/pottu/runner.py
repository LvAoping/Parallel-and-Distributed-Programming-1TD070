#!/usr/bin/python3

import sys
import random
import os

NUMRUNS = 100

LB = -1000000
UB =  1000000
SIZES = range(1, 400)
NUMPE = range(1, 17)

def genFile(n):
  with open("input/{}x{}.txt".format(str(n), str(n)), "w") as f:
    f.write(str(n) + '\n')
    
    nums = random.choices(range(LB, UB), k=n*n)

    for row in range(n):
      for col in range(n):
        f.write(str(nums.pop()) + ' ')
      f.write('\n')

for i in range(1, NUMRUNS):
  p = random.choice(NUMPE)
  n = random.choice(SIZES) * p
  genFile(n)
  c = "mpirun -n {} ./shearsort -cs input/{}x{}.txt".format(p, n, n)
  print("Running: {}".format(c))
  os.system(c)

