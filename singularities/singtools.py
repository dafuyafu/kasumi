from pys.pytools import validate_type, tuple_union
from basics.basictools import it_symbols
from basics.polytools import poly, Poly, diff, solve, polynomialring, factor
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
	mod = f.get_modulus()
	k = kp = f.get_coeff_dom()
	pd_dict, sol, it_vars, var_list = dict(), list(), it_symbols("a"), list()
	for v in f.indet_vars:
		pd_dict[v] = diff(f, v)
		if pd_dict[v] == 0:
			sol.append([p for p in k.it_elements()])
			continue
		elif pd_dict[v].is_constant:
			return list()

		"""
			* Caution!!
			This function only support the case that pd_dict[v] is univariate for any v in indet_vars.
		"""

		if not pd_dict[v].is_univariate():
			raise NotImplementedError("the case that one of partial differentials of f is multivariate has not been implemented yet")

		sol_ = solve(pd_dict[v])
		if sol_:
			sol.append(sol_)
		else:
			var = next(it_vars)
			deg = pd_dict[v].degree(v)
			if k.is_prime():
				var_list.insert(0, var)
				rel = pd_dict[v].subs({v: var}).as_dp().sort_vars(var_list)
			else:
				rep = pd_dict[v].rep.sort_vars(tuple_union((v, ), tuple(var_list)))
				poly_ = poly(rep, v, coeff_dom=k)
				var_list.insert(0, var)
				rel = factor(poly_, deg = deg // math.gcd(deg, k.ex_degree()))[0].as_dp().subs({v: var}).sort_vars(var_list)
			k = k.extend(var, rel=rel)
			sol_ = [poly(var.as_dp(), *f.indet_vars, coeff_dom=k) ** k.mod ** i for i in range(deg)]
			print("sol_%s = %s" % (str(v), str(sol_)))
			sol.append(sol_)
	sing_ = list()
	for p in itertools.product(*sol):
		d = dict(zip(f.indet_vars, p))
		sing_.append(d)
		# if f_.subs(d) == 0:
		# 	sing_.append(d)
	return sing_, k
