from pys.pytools import validate_type
from basics.basictools import it_symbols
from basics.polytools import poly, Poly, diff, solve, polynomialring
from geometries.geotools import point
import itertools
import math

def sing(f, as_point=True):

	"""
		# sing()
		return singular locus of f over algebraic closure of f.coeff_dom


		## Example:
		>>> f = poly(x**2 - y**3, mod = 3)
		>>> sing(f)
		[(0, 0)]
	"""

	validate_type(f, Poly)
	mod, k = f.get_modulus(), f.get_coeff_dom()
	pd_dict, sol, it_vars, var_list = dict(), list(), it_symbols("a"), list()
	for v in f.indet_vars:
		pd_dict[v] = diff(f, v)

		"""
			* Caution!!
			This function only support the case that pd_dict[v] is univariate for any v in indet_vars.
		"""

		if not pd_dict[v].is_univariate():
			raise NotImplementedError("the case that one of partial differentials of f is multivariate has not been implemented yet")

		sol_ = [p for p in k.it_elements() if f.subs({pd_dict[v].get_univariate(): p}) == 0]
		if sol_:
			sol.append(sol_)
		else:
			var = next(it_vars)
			var_list.append(var)
			if k.is_prime():
				rel = pd_dict[v].subs({v: var}).as_dp().sort_vars(var_list)
			else:
				rel = factor(poly(pd_dict[v].rep, coeff_dom = k), deg = pd_dict[v].degree(v) // math.gcd(deg, k.ex_degree()))[0].subs({v: var}).as_dp().sort_vars(var_list)
			k = k.extend(var, rel=rel)
			print(repr(k))
			sol_ = [k.reduce(var ** k.mod ** i) for i in range(rel.degree())]
			sol.append(sol_)
	sing_ = list()
	f_ = poly(f.rep, coeff_dom=k, quo=f.dom.quo)
	for p in itertools.product(*sol):
		d = dict()
		for i in range(len(f_.indet_vars)):
			d[f_.indet_vars[i]] = p[i]
		if f_.subs(d) == 0:
			sing_.append(d)
	return sing_
