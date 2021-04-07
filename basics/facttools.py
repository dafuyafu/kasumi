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
		return sffgcd(g, q) * (-1)

def factor_equal_degree(f, d):
	validate_type(f, Poly)
	r = f.degree() // d
	q = f.number_of_domain_elements()
	F = [f]
	var = f.indet_vars[0]
	while len(F) < r:
		g = SFFQuotientPoly.random_poly(var, f.dom, r * d, f.rep)
		g = sffpoly((g ** ((q ** d - 1) // 2) - 1).rep, f.dom)
		F_1 = []
		while len(F) > 0:
			h = F.pop(-1)
			z = sffgcd(h, g)
			if z == 1 or z == h:
				F_1.append(h)
			else:
				F_1.append(z)
				F_1.append(h // z)
		F = F_1
	return F