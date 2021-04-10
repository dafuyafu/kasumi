from pys.pytools import validate_type, tuple_union, tuple_minus
from basics.polytools import poly, Poly, Symbol

def diff(f, *var):
	validate_type(f, Poly)
	if len(var) > 1:
		raise ValueError("the number of variable argument must be 0 or 1, not %s" % str(len(var)))
	elif len(var) == 1:
		var_ = var[0]
		sorted_vars = tuple_union(var, tuple_minus(f.inner_vars, var))
		rep_ = f.rep.sort_vars(sorted_vars)
	else:
		var_ = f.indet_vars[0]
		rep_ = f.rep

	return poly(rep_.diff(var_), *f.indet_vars, dom=f.dom)