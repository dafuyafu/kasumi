from pys.pytools import validate_type
from basics.basictools import dp, DP
from basics.polytools import poly, Poly

def gcd(f, g):
	validate_type(f, Poly)
	validate_type(g, Poly)
	q = f % g
	if q == 0:
		return g
	else:
		return gcd(g, q) * (-1)

def factor(f):
	"""
	returns irreducible component over f.coeff_dom with Cantor-Zassenhaus algorithm
	"""
	pass

def factor_equal_degree(f, d):
	validate_type(f, Poly)
	r = f.degree() // d
	q = f.dom.number()
	F = [f]
	var = f.indet_vars[0]
	while len(F) < r:
		g = f.dom.random_poly(var, r * d)
		g = poly((g ** ((q ** d - 1) // 2) - 1), dom=f.dom)
		F_1 = []
		while len(F) > 0:
			h = F.pop()
			z = gcd(h, g)
			if z == 1 or z == h:
				F_1.append(h)
			else:
				F_1.append(z)
				F_1.append(h // z)
		F = F_1
	return F