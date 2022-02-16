import math
from unittest import TestCase


def square_me(x):
    return x**2


def test_1():
    # test the simple square_me function
    assert square_me(2.) == 4.


def test_2():
    # testing without an absolute equality
    assert abs(square_me(2.34) - 5.4756) < 0.001

    
## can also make a class which is a subclass of TestCase...
class TestOneNumber(TestCase):
    def test_floats(self):
        for num in [1617161771.7650001, math.pi, math.pi**100,
                    math.pi**-100, 3.1]:
            self.assertEqual(square_me(num), pow(num,2))
