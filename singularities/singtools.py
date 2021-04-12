from basics.basictools import it_symbols
from basics.polytools import poly, Poly, diff, solve
from geometries.geotools import point

def sing(f, as_point=True):

	"""

	return singular locus of f over algebraic closure of f.coeff_dom

	Example:
	>>> f = poly(x**2 - y**3, mod = 3)
	>>> sing(f)
	[(0, 0)]

	"""

	validate_type(f, Poly)
	mod, k = f.get_modulus(), f.get_coeff_dom()
	pd_dict, sol, it_vars = dict(), list(), it_symbols("a")
	for v in f.indet_vars:
		pd_dict[v] = diff(f, v)
		sol_ = [p for p in k.it_points(f.indet_vars) if f.subs(p) == 0]
		if sol_:
			sol.append(sol_)
		else:
			pass
			
	return sol_f, sff_.as_SFF()