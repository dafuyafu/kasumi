from basics import dp, DP
from pytools import tuple_minus, tuple_union

class Poly:
	"""
	represents an element of polynomial ring

	"""

	def __new__(cls, rep, *var, **options):
		if isinstance(rep, int):
			self = super().__new__(Integer)
		elif isinstance(rep, DP):
			if set(rep.inner_vars) > set(var):
				self = super().__new__(cls)
			elif tuple_intersection(rep.inner_vars, var):
				self = super().__new__(Constant)
			else:
				raise ValueError("variable argument must be some of rep's or others")
		elif isinstance(rep, Poly):
			if set(rep.var) > set(var):
				self = super().__new__(cls)
			elif tuple_intersection(rep.var, var):
				self = super().__new__(Constant)
			else:
				raise ValueError("variable argument must be some of rep's or others")
		else:
			raise TypeError("rep must be int, DP, Poly, Contant or Integer object, not %s" % rep.__class__.__name__)
		return self

	def __init__(self, rep, *var, **options):
		if isinstance(rep, Poly):
			rep = rep.rep
		elif isinstance(rep, DP):
			pass
		elif isinstance(rep, int):
			if "mod" in options:
				rep = dp(var[0], rep, mod=options["mod"])
			else:
				rep = dp(var[0], rep)
		else:
			raise TypeError("rep must be Poly or DP, not %s" % rep.__class__.__name__)
		self.rep = rep
		self.inner_vars = rep.inner_vars
		self.indet_vars, self.const_vars = self._set_var(rep, var)
		self.deg = rep.deg

		if "mod" in options:
			self.mod = options["mod"]
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
		self.options = {"mod": self.mod, "rel": self.rel, "quo": self.quo}

	def _set_var(self, rep, var):
		if len(var) == 0:
			return (rep.inner_vars, tuple())
		else:
			if set(rep.inner_vars) >= set(var):
				return (var, tuple_minus(rep.inner_vars, var))
			else:
				raise ValueError("variable must be of its rep")

	def __repr__(self):
		repr_ = "%s(%s, %s" % (self.__class__.__name__, str(self), self.var)
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

	def __add__(f, g):
		if isinstance(g, int):
			return poly(f.rep + g, *f.var, **self.options)
		elif isinstance(g, Poly):
			return poly(f.rep + g.rep, *f.var, **self.options)

	"""
	Reprisentation Functions
	"""

	def as_dist_rep(self):
		return self.rep.as_dist_rep()



def poly(f, *var, **options):
	return Poly(f, *var, **options)

class Constant(Poly):
	pass

class Integer(Constant):
	pass

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
	* Iterable magic methods
	"""
	def __iter__(self):
		return iter(self.rel)

	def __len__(self):
		return len(self.rel)

def reduction(f, rel):
	validate_type(f, DP)
	for r in rel:
		if not r['var'] in r.inner_vars:
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