#!/usr/bin/python3

import sys
import random



try:
  n = int(sys.argv[1])
except:
  print("Usage: ./gen.py <side length>")
  quit()

with open("input/{}x{}.txt".format(str(n), str(n)), "w") as f:
  f.write(str(n) + '\n')
  
  nums = random.choices(range(-n, n), k=n*n)

  for row in range(n):
    for col in range(n):
      f.write(str(nums.pop()) + ' ')
    f.write('\n')