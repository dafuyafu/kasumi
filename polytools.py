from basics import dp, DP
from pytools import tuple_minus, tuple_union

class Poly:
	"""
	represents an element of polynomial ring

	
	"""

	def __new__(cls, rep, mod=0, quo=0, rel=0):
		if isinstance(rep, int):
			self = super().__new__(Integer)
		elif isinstance(rep, DP):
			self = super().__new__(cls)
		elif isinstance(rep, Poly):
			pass
		else:
			raise TypeError("rep must be int, DP, Poly, Contant or Integer object, not %s" % rep.__class__.__name__)
		return self

	def __init__(self, rep, *var, mod=0, quo=0, rel=0):
		self.var, self.const_var = self._set_var(rep, var)
		self.rep = rep.sort_vars(self.var + self.const_var)
		self.inner_vars = self.rep.inner_vars
		self.deg = rep.deg
		self.mod = mod
		self.quo = quo
		self.rel = rel

	def _set_var(self, rep, var):
		if len(var) == 0:
			return rep.inner_vars
		else:
			if set(rep.inner_vars) >= set(var):
				return (var, tuple_minus(rep.inner_vars, var))
			else:
				raise ValueError("variable must be of its rep")

	def __repr__(self):
		repr_ = "%s(%s, %s" % (self.__class__.__name__, str(self.rep), self.var)
		if modulus > 0:
			repl_ += ", mod = %s" % self.mod
		if quo != 0:
			repl_ += quo

	def __str__(self):
		return str(self.rep)

	"""
	Binary operators:
	add, sub, mul, floordiv, mod, pow
	"""

	def __add__(self):
		pass

class Constant(Poly):
	pass

class Integer(Constant):
	pass

def poly(f):
	return Poly(f)