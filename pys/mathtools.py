from math import sqrt, floor
import itertools

def is_prime(n):
	if n == 0 or n == 1 or n % 2 == 0:
		return False
	elif n == 2:
		return True
	else:
		sqrt_ = floor(sqrt(n)) + 1
		for i in range(3, sqrt_, 2):
			if n % i == 0:
				return False
		return True

def hcomb(n, deg):
	t = (0, )*n
	yield t
	for i in range(1, deg):
		s = tuple(i for i in range(n))
		for p in itertools.combinations_with_replacement(s, i):
			yield tuple(p.count(v) for v in s)