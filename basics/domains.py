from pys.mathtools import is_prime
from pys.pytools import tuple_or_object, validate_type
from basics.basictools import symbols, dp, DP
from abc import ABCMeta, abstractmethod

"""
* Relation class and functions
We construct Ring and PolynomialRing object with some Relation object
"""

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
		self.rel_list = tuple(rel_list)
		self.var_list = tuple([rel["var"] for rel in self.rel_list])

	"""
	* Representation magic methods
	repr, str
	"""

	def __repr__(self):
		return "Relation(%s)" % str(self.rel_list)

	def __str__(self):
		if len(self) == 1:
			return str(self[0]["rep"])
		else:
			str_, flag_first = str(), True
			for rel in self:
				if flag_first:
					flag_first = False
				else:
					repr_ += ", "
				str_ += str(rel["rep"])
			return str_
	"""
	* Iterable magic methods
	"""
	def __iter__(self):
		return iter(self.rel_list)

	def __len__(self):
		return len(self.rel_list)

	def __getitem__(self, key):
		return self.rel_list[key]

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
			if r.rel_list == s.rel_list:
				return True
			else:
				return False
		else:
			raise TypeError("unsupported operand type(s) for '=='")

def relation(reps, mod=0):
	return Relation(reps, mod)

"""
* Reduction methods
Reduce variables with variable relations.
"""

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

class Ring(metaclass=ABCMeta):
	"""
	Represent coefficient rings of Poly
	Ring and some objects that inherit it have only relational data and generators of their elements.
	Note that they are nothing to do with calculations of themselves.

	Example:
	>>> Ring()
	ZZ
	>>> Ring(complex=True)
	CC
	>>> Ring(rational=True)
	QQ
	>>> Ring(mod=4) # ZZ/(4)
	FiniteRing(mod: 4)
	>>> Ring(mod=3) # ZZ/(3)
	FiniteField(mod: 3)
	>>> Ring(mod=3, rel=Relation(a**2 - 2))
	FiniteField(mod: 3, rel: a**2 - 2, var: (a, ))

	"""

	is_domain = True
	is_finite = False
	is_field = False
	is_alg_closed = False

	def __new__(cls, *var, **options):
		if not options:
			self = super().__new__(ZZ)
		else:
			if "complex" in options:
				if options["complex"]:
					self = super().__new__(CC)
			elif "rational" in options:
				if options["rational"]:
					self = super().__new__(QQ)
			elif "mod" in options:
				validate_type(options["mod"], int)
				if options["mod"] < 0:
					raise ValueError("modulus option must be 0 or positive integer")
				elif options["mod"] > 0:
					""" positive characteristic """
					if is_prime(options["mod"]):
						self = super().__new__(FiniteField)
					else:
						self = super().__new__(FiniteRing)
				else:
					self = super().__new__(ZZ)
			else:
				raise ValueError("unsupported option '%s'" % next(iter(options)))
		return self
	
	def __init__(self, *var, **options):
		self.options = dict()
		if "mod" in options:
			self.mod = self.options["mod"] = options["mod"]
		else:
			self.mod = 0

		if "rel" in options:
			validate_type(options["rel"], Relation)
			self.rel = self.options["rel"] = options["rel"]
		else:
			self.rel = 0

		if len(var) > 0:
			self.var = self.options["var"] = var
		else:
			self.var = tuple()

	def __repr__(self):
		repr_ = "%s" % self.__class__.__name__
		if self.options:
			repr_ += "("
			flag_first = True
			for k in self.options:
				if flag_first:
					flag_first = False
				else:
					repr_ +=", "
				repr_ += str(k) + ": " + str(self.options[k])
			repr_ += ")"
		else:
			pass
		return repr_

	def __str__(self):
		return repr(self)

	@abstractmethod
	def number(self):
		raise NotImplementedError("This class is abstract class")

	@abstractmethod
	def has_relation(self):
		raise NotImplementedError("This class is abstract class")


def ring(**options):
	return Ring(**options)

class ZZ(Ring):
	"""

	* ZZ

	represent the ring of rational integers
	given no options, the constructor generates ZZ object

	"""

	def number(self):
		raise TypeError("%s is infinte" % self.__class__.__name__)

	"""
	* Validators
	"""

	def has_relation(self):
		return False

	def has_variable(self):
		return bool(self.var)

class Field(Ring):
	is_field = True

class QQ(ZZ):
	pass

class CC(ZZ):
	is_alg_closed = True

class FiniteRing(Ring):
	"""

	* FiniteRing

	represent finite rings such as ZZ/(6) or ZZ/(10)[a]/(a**2 - 3)

	"""
	is_finite = True

	def number(self):
		pass

	def has_relation(self):
		return "rel" in self.options

class FiniteField(Field, FiniteRing):
	"""

	* FiniteField

	represent finite fields such as ZZ/(5) or ZZ(5)[a]/(a**2 - 2)

	"""
	pass

class ACFiniteField(Field):
	is_alg_closed = True

class PolynomialRing:
	pass