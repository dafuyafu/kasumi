from pys.pytools import tuple_intersection, tuple_minus, tuple_union, tuple_or_object, validate_type
from basics.basictools import dp, DP, Symbol, int_to_dp
from basics.domains import ring, Ring, polynomialring, PolynomialRing, Relation
import itertools as it

class Poly:
	"""
	represents an element of polynomial ring

	Example:
	>>> x = symbols("x")
	>>> Poly(x**2 - 2*x + 1)
	Poly(x**2 - 2*x + 1, x)
	>>> Poly(x ** 2 - 2*x + 1, mod=3)
	Poly(x**2 - 2*x + 1, x, mod=3)
	>>> from basics.domains import ring
	>>> ff = ring(mod=3, a**2 - 2)
	>>> Poly(x**2 - 2*x + 1, dom=ff)
	Poly(x**2 - 2*x + 1, x, mod=3, rel=a**2 - 2)

	"""

	def __new__(cls, rep, *var, **options):
		if "dom" in options:
			if isinstance(rep, (DP, Poly, Symbol)):
				if isinstance(rep, DP):
					rep = options["dom"].reduce(rep)
				elif isinstance(rep, Poly):
					rep = rep.reduce()
				else:
					rep = options["dom"].reduce(rep.as_dp())
				if rep.is_constant():
					rep = rep.get_constant()
		if isinstance(rep, int):
			self = super().__new__(Integer)
		elif isinstance(rep, Symbol):
			if rep in var:
				self = super().__new__(Poly)
			else:
				self = super().__new__(Constant)
		elif isinstance(rep, DP):
			if not rep.degree(*var):
				self = super().__new__(Constant)
			else:
				self = super().__new__(Poly)
		else:
			raise TypeError("rep must be int, Symbol, DP, Poly, Contant or Integer object, not %s" % rep.__class__.__name__)
		return self

	def __init__(self, rep, *var, **options):
		if isinstance(rep, Poly):
			self.rep = rep.rep
			self.indet_vars = rep.indet_vars
			self.const_vars = rep.const_vars
			self.inner_vars = rep.inner_vars
			self.dom = rep.dom
		else:
			if isinstance(rep, DP):
				pass
			elif isinstance(rep, int):
				rep = int_to_dp(rep, *var)
			elif isinstance(rep, Symbol):
				rep = rep.as_dp()
			else:
				raise TypeError("rep must be Poly or DP, not %s" % rep.__class__.__name__)
			self.indet_vars, self.const_vars = self._set_var(rep, var)
			self.inner_vars = self.rep.inner_vars
			self.deg = rep.deg

			"""

			* Set domain options

			Given "dom" option, initialize with it.
			Not given, construct PolynomialRing object and set to "dom" option

			"""

			if "dom" in options:
				validate_type(options["dom"], PolynomialRing)
				self.dom = options["dom"]
			else:
				if "mod" in options:
					validate_type(options["mod"], int)
					mod = options["mod"]
				else:
					mod = 0

				if "rel" in options:
					validate_type(options["rel"], int, tuple, DP, Relation)
					rel = options["rel"]
				else:
					rel = 0

				self.coeff_dom = ring(*self.const_vars, mod=mod, rel=rel)

				if "quo" in options:
					validate_type(options["rel"], int, tuple, DP, Relation)
					quo = options["rel"]
				else:
					quo = 0

				self.dom = polynomialring(*self.indet_vars, coeff_dom=self.coeff_dom, quo=quo)
			self.rep = self.reduce()

	def _set_var(self, rep, var):
		if len(var) == 0:
			self.rep = rep
			return (rep.inner_vars, tuple())
		else:
			self.rep = rep.sort_vars(tuple_union(var, rep.inner_vars))
			return var, tuple_minus(rep.inner_vars, var)

	def __repr__(self):
		repr_ = "%s(%s, %s" % (self.__class__.__name__, str(self), tuple_or_object(self.indet_vars))
		if self.dom.mod > 0:
			repr_ += ", mod = %s" % self.dom.mod
		if self.dom.rel != 0:
			repr_ += ", rel = %s" % str(self.dom.rel)
		if self.dom.quo != 0:
			repr_ += ", quo = %s" % str(self.dom.quo)
		repr_ += ")"
		return repr_

	def __str__(self):
		return self.as_dist()

	"""
	Binary operators:
	add, sub, mul, floordiv, mod, pow
	"""

	def reduce(self):
		return self.dom.reduce(self.rep)

	def __add__(f, g):
		f_, g_ = binary_uniform(f, g)
		if f.dom.mod:
			add_ = (f_ + g_) % f.dom.mod
		else:
			add_ = f_ + g_
		return poly(add_, *f.indet_vars, dom=f.dom)

	def __radd__(f, g):
		return f + g

	def __sub__(f, g):
		return f + (- g)

	def __rsub__(f, g):
		return - f + g

	def __mul__(f, g):
		f_, g_ = binary_uniform(f, g)
		if f.dom.mod:
			mul_ = (f_ * g_) % f.dom.mod
		else:
			mul_ = f_ * g_
		return poly(mul_, *f.indet_vars, dom=f.dom)

	def __rmul__(f, g):
		return f * g

	def __neg__(f):
		return f * (-1)

	def __truediv__(f, g):
		if isinstance(g, int):
			if f.mod > 0:
				return poly(f.rep.div(g, f.dom.mod), *f.indet_vars, dom=f.dom)
			else:
				raise TypeError("can divide only on field")
		else:
			raise TypeError("unsupported operand type(s) for '/'")

	def __floordiv__(f, g):
		if isinstance(g, int):
			if f.mod > 0:
				return poly(f.rep.div(g, f.dom.mod), *f.indet_vars, dom=f.dom)
			else:
				raise TypeError("can divide only on fields")

		"""
		* Validation part
		validate f (and g) is univariate and f has the same indeterminate as g
		"""

		if not f.is_univariate():
			raise TypeError("operands must be univariate polynomial")
		if not f.indet_vars == g.indet_vars:
			raise TypeError("operands must have the same variables")

		if f.degree() < g.degree():
			return poly(0, *f.indet_vars, dom=f.dom)

		"""
		* Calculation part
		"""

		q, r, v = poly(0, *f.indet_vars, dom=f.dom), f, poly(f.indet_vars[0], *f.indet_vars, dom=f.dom)
		while r.degree() >= g.degree():
			t = v ** (r.degree() - g.degree()) * (r.LC() / g.LC())
			r = r - t * g
			q = q + t
		return q

	def __mod__(f, g):
		return f - (f // g) * g

	def __pow__(f, e):
		validate_type(e, int)
		if e < 0:
			raise ValueError("exponent must be positive")
		if e == 0:
			return poly(1, *f.indet_vars, dom=f.dom)
		num_ = bin(e).replace('0b', '')
		len_ = len(num_)
		list_ = [len_ - d - 1 for d in range(len_) if num_[d] == '1']
		pow_ = 1
		for l in list_:
			pow_ *= pow_self(f, l)
		return pow_

	def __eq__(f, g):
		if isinstance(g, Poly):
			if f.rep == g.rep:
				return True
			else:
				return False
		elif isinstance(g, DP) or isinstance(g, int):
			if f.rep == g:
				return True
			else:
				return False
		elif isinstance(g, Symbol):
			if f.rep == as_dp(g):
				return True
			else:
				return False
		else:
			raise TypeError("unsupported operand type(s) for '=='")
	
	"""
	* Iterable magic methods
	"""

	def __getitem__(self, item):
		return self.rep[item]

	def __setitem__(self, key, value):
		raise NotImplementedError()

	def __iter__(self):
		return iter(self.rep.coeffs)

	def __len__(self):
		return len(self.rep)

	"""
	* Iterator methods
	"""

	def it_dist(self, termorder="lex", with_index=False):
		"""
		generate coefficients of each index of exponents,
		return int or DP object with index tuple if with_index is True.
		"""

		validate_type(termorder, str)
		if termorder == "lex":
			iters = [range(self.degree(v, non_negative=True), -1, -1) for v in self.indet_vars]
			for p in it.product(*iters):
				c = self.rep.get(p, as_int=False)
				if c == 0:
					continue
				if with_index:
					yield (p, c)
				else:
					yield c
		elif termorder == "grevlex":
			raise NotImplementedError()
		else:
			raise TypeError("given unsupported termorder '%s'" % termorder)

	def it_reversed_dist(self):
		pass

	"""
	* Reprisentation methods
	"""

	def as_dist(self, termorder="lex", total=False):
		if total:
			return self.rep.as_dist(termorder)
		else:
			rep_, lt, flag_one, flag_paren = str(), True, False, False
			if self.is_zero():
				return "0"
			for data in self.it_dist(termorder=termorder, with_index=True):
				if data[1] == 0:
					continue

				# plus minus part
				if lt:
					lt = False
					if data[1] < 0:
						rep_ += "- "
				else:
					if data[1] > 0:
						rep_ += " + "
					else:
						rep_ += " - "

				# DP bracket part
				if isinstance(data[1], DP):
					if data[1].is_monomial() or not any(data[0]):
						pass
					else:
						flag_paren = True
						rep_ += "("

				# coeffcients part
				if data[1] == 1 or data[1] == -1:
					flag_one = True
				else:
					flag_one = False
					rep_ += str(abs(data[1]))
					if flag_paren:
						flag_paren = False
						rep_ += ")"
					if any(data[0]):
						rep_ += "*"

				# variables part
				for i in range(len(self.indet_vars)):
					if data[0][i] > 1:
						rep_ += str(self.indet_vars[i]) + "**" + str(data[0][i])
					elif data[0][i] == 1:
						rep_ += str(self.indet_vars[i])
					else:
						continue
					if i < len(self.indet_vars) - 1 and any(data[0][(i + 1):]):
						rep_ +="*"

				# last one modification
				if flag_one and not any(data[0]):
					rep_ += "1"

			return rep_

	def as_list(self):
		return self.rep.as_list()

	def as_dp(self):
		return self.rep

	"""
	* Degree methods
	"""

	def degree(self, *var, total=False, as_dict=False, any_vars=False, non_negative=False):
		if total:
			d, i = 0, len(self.indet_vars)
			for v in self.indet_vars:
				d_ = self.rep.degree(v)
				if d_ == -1:
					i -= 1
				else:
					d += d_
			if i == 0:
				if non_negative:
					return 0
				else:
					return -1
			else:
				return d
		else:
			return self.rep.degree(*var, total=total, as_dict=as_dict, any_vars=any_vars, non_negative=non_negative)

	"""
	* Validators
	These functions return bool.
	"""

	def is_univariate(self):
		if len(self.indet_vars) == 1:
			return True
		else:
			return False

	def is_zero(self):
		return self.rep.is_zero()

	"""
	* Other methods
	"""

	def LC(self, termorder="lex"):
		for d in self.it_dist(termorder=termorder):
			if d != 0:
				return poly(d, *self.indet_vars, dom=self.dom)
			else:
				continue
		return poly(0, *self.indet_vars, dom=self.dom)

def poly(f, *var, **options):
	return Poly(f, *var, **options)

class Constant(Poly):
	"""
	* Binary operations
	"""

	def __truediv__(f, g):
		if f.mod == 0:
			raise TypeError("cannot divide not on field")
		if isinstance(g, Constant):
			return f * g ** (f.mod - 2)
		elif isinstance(g, int):
			return f / g

	def __floordiv__(f, g):
		pass

	"""
	* Validators
	These functions return bool.
	"""

	def is_univariate(self):
		if super().is_univariate():
			return True
		else:
			if len(self.const_vars) == 1:
				return True
			else:
				return False

	"""
	* Degree methods
	"""

	def degree(self, *var, total=False, as_dict=False, any_vars=False, non_negative=False):
		return 0

class Integer(Constant):

	"""
	* Validators
	These functions return bool.
	"""

	def is_univariate(self):
		return False

	"""
	* Degree methods
	"""

	def degree(self, *var, total=False, as_dict=False, any_vars=False, non_negative=False):
		if non_negative:
			return 0
		else:
			if self.rep == 0:
				return -1
			else:
				return 0

def binary_uniform(f, g):
	if isinstance(g, int) or isinstance(g, DP):
		f_, g_ = f.rep, g
	elif isinstance(g, Poly):
		f_, g_ = f.rep, g.rep
	else:
		raise TypeError("unsupported operand type(s) for '+'")
	return f_, g_

def pow_self(f, e):
	for i in range(e):
		f = f * f
	return f