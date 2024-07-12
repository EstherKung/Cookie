#cooking
import numpy as np
from bisect import bisect_left, bisect_right

l = [0.04, 0.08, 0.09, 0.30, 0.90, 0.2, 0.5, 0.7, 1.0]

check = 0.2
pos = bisect_left(l, check)
print(pos)
if pos == 0:
    print(l[0])
if pos == len(l):
    print(l[-1])

before = l[pos-1]
print(before)
if check == l[pos]:
    after = l[pos + 1]
else: after = l[pos]
print(after)
print(l[pos])
if after - check < check - before:
    print(after)
else: print(before)

# def nearest(des: str, coords, xloc):
#     pos = bisect_left(coords, xloc)
#     if pos == 0:
#         return(coords[0])
#     if pos == len(coords):
#         return(coords[-1])
#     before = coords[pos - 1]
#     after = coords[pos]
    
#     match des:
#         case 'prev': 
#             return before
#         case 'next':
#             return after

# print("test this: " , nearest('next', l, 0.3))