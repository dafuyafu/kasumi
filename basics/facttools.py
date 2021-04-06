from basics.pytools import validate_type
from basics.basictools import dp, DP
from basics.polytools import poly, Poly

def gcd(f, g):
	validate_type(f, Poly)
	validate_type(g, Poly)
	q = f % g
	if q == 0:
		return g
	else:
		return sffgcd(g, q) * (-1)