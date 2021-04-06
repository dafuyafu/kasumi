from basics.basictools import dp, DP, Symbol, int_to_dp
from basics.pytools import tuple_intersection, tuple_minus, tuple_union, validate_type
from math import floor, sqrt
import itertools as it

class Poly:
	"""
	represents an element of polynomial ring

	"""

	def __new__(cls, rep, *var, **options):
		if isinstance(rep, int):
			self = super().__new__(Integer)
		elif isinstance(rep, Symbol):
			if rep in var:
				self = super().__new__(cls)
			else:
				self = super().__new__(Constant)
		elif isinstance(rep, DP):
			if set(rep.inner_vars) >= set(var):
				self = super().__new__(cls)
			elif not tuple_intersection(rep.inner_vars, var):
				self = super().__new__(Constant)
			else:
				raise ValueError("variable argument must be some of rep's or others")
		elif isinstance(rep, Poly):
			if set(rep.var) >= set(var):
				self = super().__new__(cls)
			elif tuple_intersection(rep.var, var):
				self = super().__new__(Constant)
			else:
				raise ValueError("variable argument must be some of rep's or others")
		else:
			raise TypeError("rep must be int, Symbol, DP, Poly, Contant or Integer object, not %s" % rep.__class__.__name__)
		return self

	def __init__(self, rep, *var, **options):
		if isinstance(rep, Poly):
			rep = rep.rep
		elif isinstance(rep, DP):
			pass
		elif isinstance(rep, int):
			rep = int_to_dp(rep, *var)
		elif isinstance(rep, Symbol):
			rep = rep.as_dp()
		else:
			raise TypeError("rep must be Poly or DP, not %s" % rep.__class__.__name__)
		self.rep = rep
		self.inner_vars = rep.inner_vars
		self.indet_vars, self.const_vars = self._set_var(rep, var)
		self.deg = rep.deg

		if "mod" in options:
			validate_type(options["mod"], int)
			if options["mod"] == 0:
				self.mod = 0
			elif is_prime(options["mod"]):
				self.mod = options["mod"]
			else:
				raise ValueError("mod must be a prime number or 0")
		else:
			self.mod = 0

		if "rel" in options:
			if isinstance(options["rel"], Relation):
				self.rel = options["rel"]
			elif isinstance(options["rel"], DP) or isinstance(options["rel"], list) or isinstance(options["rel"], tuple):
				self.rel = Relation(options["rel"], self.mod)
			elif isinstance(options["rel"], int):
				self.rel = 0
			else:
				raise TypeError("rel option must be 0(int), DP or Relation, not %s" % options["rel"].__class__.__name__)
		else:
			self.rel = 0

		if "quo" in options:
			if isinstance(options["quo"], Relation):
				self.quo = options["quo"]
			elif isinstance(options["quo"], DP) or isinstance(options["quo"], list) or isinstance(options["quo"], tuple):
				self.quo = Relation(options["quo"], self.mod)
			elif isinstance(options["quo"], int):
				self.quo = 0
			else:
				raise TypeError("quo option must be 0(int), DP or Relation, not %s" % options["quo"].__class__.__name__)
		else:
			self.quo = 0

		self.rep = self._reduce(self.rep)
		self.options = {"mod": self.mod, "rel": self.rel, "quo": self.quo}

	def _set_var(self, rep, var):
		if len(var) == 0:
			return (rep.inner_vars, tuple())
		else:
			if set(rep.inner_vars) >= set(var):
				return var, tuple_minus(rep.inner_vars, var)
			else:
				union_ = tuple_union(var, rep.inner_vars)
				self.rep = rep.sort_vars(union_)
				return var, tuple_minus(rep.inner_vars, var)

	def get_var(self):
		if len(self.indet_vars) == 1:
			return self.indet_vars[0]
		else:
			return self.indet_vars

	def __repr__(self):
		repr_ = "%s(%s, %s" % (self.__class__.__name__, str(self), self.get_var())
		if self.mod > 0:
			repr_ += ", mod = %s" % self.mod
		if self.rel != 0:
			repr_ += ", rel = %s" % str(self.rel)
		if self.quo != 0:
			repr_ += ", quo = %s" % str(self.quo)
		repr_ += ")"
		return repr_

	def __str__(self):
		return str(self.rep)

	"""
	Binary operators:
	add, sub, mul, floordiv, mod, pow
	"""

	def _reduce(self, rep):
		if self.rel:
			rep = reduction(rep, self.rel)
		if self.quo:
			rep = reduction(rep, self.quo)
		if self.mod:
			rep = rep % self.mod
		return rep

	def __add__(f, g):
		f_, g_ = binary_uniform(f, g)
		if f.mod:
			add_ = (f_ + g_) % f.mod
		else:
			add_ = f_ + g_
		return poly(f._reduce(add_), *f.indet_vars, **f.options)

	def __radd__(f, g):
		return f + g

	def __sub__(f, g):
		return f + (- g)

	def __rsub__(f, g):
		return - f + g

	def __mul__(f, g):
		f_, g_ = binary_uniform(f, g)
		if f.mod:
			mul_ = (f_ * g_) % f.mod
		else:
			mul_ = f_ * g_
		return poly(f._reduce(mul_), *f.indet_vars, **f.options)

	def __rmul__(f, g):
		return f * g

	def __neg__(f):
		return f * (-1)

	def __truediv__(f, g):
		if isinstance(g, int):
			if f.mod > 0:
				return poly(f.rep.div(g, f.mod), *f.indet_vars, **f.options)
			else:
				raise TypeError("can divide only on field")
		else:
			raise TypeError("unsupported operand type(s) for '/'")

	def __floordiv__(f, g):
		if isinstance(g, int):
			if f.mod > 0:
				return poly(f.rep.div(g, f.mod), *f.indet_vars, **f.options)
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
			return poly(0, *f.indet_vars, **f.options)

		q, r, v = poly(0, *f.indet_vars, **f.options), f, poly(f.indet_vars[0], *f.indet_vars, **f.options)
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
			return poly(1, *f.indet_vars, **f.options)
		num_ = bin(e).replace('0b', '')
		len_ = len(num_)
		list_ = [len_ - d - 1 for d in range(len_) if num_[d] == '1']
		pow_ = 1
		for l in list_:
			pow_ *= f._pow_self(e)
		return pow_

	def _pow_self(f, e):
		for i in range(e):
			f = f * f
		return f

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
				c = self.rep.get(p)
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

	def as_dist_rep(self):
		return self.rep.as_dist_rep()

	def as_list(self):
		return self.rep.as_list()

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

	"""
	* Other methods
	"""

	def LC(self, termorder="lex"):
		for d in self.it_dist(termorder=termorder):
			if d != 0:
				return poly(d, *self.indet_vars, **self.options)
			else:
				continue
		return poly(0, *self.indet_vars, **self.options)

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

class Relation:
	"""
	Represent quotient relations.
	"""
	def __init__(self, reps, mod=0):
		if isinstance(reps, tuple):
			pass
		elif isinstance(reps, DP):
			reps = (reps, )
		elif isinstance(reps, list):
			reps = tuple(reps)
		else:
			raise TypeError("reps must be tuple, list or a DP, not '%s'" % reps.__class__.__name__)
		rel_list = list()
		for r in reps:
			if r == 0:
				continue
			lc_, lt_ = r.LC(), r.LT(as_data=True)
			if mod > 0:
				if lc_ != 1:
					# make rep monic
					r = (r * pow(lc_, mod - 2, mod)) % mod
				else:
					pass
			else:
				pass
			var = r.inner_vars[next(i for i in range(len(lt_[0])) if lt_[0][i] != 0)]
			rel_list.append({'var': var,
							 'rep': r,
							 'deg': r.deg})
		self.rel = tuple(rel_list)

	"""
	* Representation magic methods
	repr, str
	"""

	def __repr__(self):
		return "Relation(%s)" % str(self.rel)

	def __str__(self):
		if len(self) == 1:
			return "%s" % str(self.rel[0]["rep"])
		else:
			return "%s" % tuple([str(r["rep"]) for r in self.rel])
	"""
	* Iterable magic methods
	"""
	def __iter__(self):
		return iter(self.rel)

	def __len__(self):
		return len(self.rel)

	def __getitem__(self, key):
		return self.rel[key]

	"""
	* Binary Operations
	"""

	def __bool__(self):
		if len(self) == 0:
			return False
		elif len(self) == 1 and self[0]['rep'] == 0:
			return False
		else:
			return True

	def __eq__(r, s):
		if isinstance(s, int):
			if len(r) == 1 and r[0]['rep'] == s:
				return True
			else:
				return False
		elif isinstance(s, Relation):
			if r.rel == s.rel:
				return True
			else:
				return False
		else:
			raise TypeError("unsupported operand type(s) for '=='")

def relation(reps, mod=0):
	return Relation(reps, mod)

def reduction(f, relation):
	validate_type(f, DP)
	for r in relation:
		if not r['var'] in f.inner_vars:
			continue
		else:
			f = _simple_red(f, r)
	return f

def _simple_red(f, r):
	lm = r['var'] ** r['deg']
	if f.degree(r['var']) < r['rep'].degree(r['var']):
		return f
	else:
		return _simple_red_rec(f, lm, lm - r['rep'], r['var'])

def _simple_red_rec(f, lm, sub, var):
	if f.degree(var) == lm.degree(var):
		return f.subs({lm: sub})
	else:
		return _simple_red_rec(f, lm * var, sub * var, var).subs({lm: sub})

def binary_uniform(f, g):
	if isinstance(g, int) or isinstance(g, DP):
		f_, g_ = f.rep, g
	elif isinstance(g, Poly):
		f_, g_ = f.rep, g.rep
	else:
		raise TypeError("unsupported operand type(s) for '+'")
	return f_, g_

def is_prime(n):
	if n == 2:
		return True
	if n % 2 == 0:
		return False
	sqrt_ = floor(sqrt(n)) + 1
	for i in range(3, sqrt_, 2):
		if n % i == 0:
			return False
	return True