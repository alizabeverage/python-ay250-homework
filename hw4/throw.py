from math import sqrt
from random import uniform

# Simulate dart throwing. 
# For each dart, find random position for it to fall
# Test if it falls within the circle by finding distance from the origin to the dart.
# Darts within 0.5 of the origin are within the circle
    
def throw(number_of_darts):
    number_of_darts_in_circle = 0
    for n in range(int(number_of_darts)):
        x, y = uniform(0,1), uniform(0,1)
        if sqrt((x - 0.5)**2 + (y - 0.5)**2) <= 0.5:
            number_of_darts_in_circle +=  1
            
    return number_of_darts_in_circle

def throw_dask(number_of_darts):
    import dask.array as da
    x = da.random.normal(size=(number_of_darts))
    y = da.random.normal(size=(number_of_darts))
    number_of_darts_in_circle = da.sum(da.sqrt((x - 0.5)**2 + (y - 0.5)**2) <= 0.5)

    return number_of_darts_in_circle.compute()
