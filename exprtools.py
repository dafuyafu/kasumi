from sympy.core.numbers import Integer
from sympy.core.expr import Expr
from sympy.core.symbol import symbols

def reduce_coefficients(f, mod):
	codict = f.as_coefficients_dict()
	colist = [c % mod for c in codict.values()]
	print(colist)
	g, i = 0, 0
	for var in codict:
		g += colist[i] * var
		i += 1
	return g