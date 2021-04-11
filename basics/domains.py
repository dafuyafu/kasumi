from pys.mathtools import is_prime
from pys.pytools import tuple_or_object, validate_type
from basics.basictools import symbols, dp, DP, 
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
		if isinstance(reps, Relation):
			self.rel_list = reps.rel_list
			self.var_list = reps.var_list
		else:
			if isinstance(reps, tuple):
				pass
			elif isinstance(reps, DP):
				reps = (reps, )
			elif isinstance(reps, list):
				reps = tuple(reps)
			else:
				raise TypeError("reps must be tuple, list, DP or Relation, not '%s'" % reps.__class__.__name__)
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
				rel_list.append({'var': var, 'rep': r, 'deg': r.degree(var)})
			self.rel_list = tuple(rel_list)
			self.var_list = tuple([rel["var"] for rel in self.rel_list])

	"""
	* Representation magic methods
	repr, str
	"""

	def __repr__(self):
		repr_ = "Relation("
		for i in range(len(self)):
			repr_ += str(self[i]["rep"])
			if not i == len(self) - 1:
				repr_ += ", "
		return repr_ + ")"

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
* Domain class and functions
"""

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
				self = super().__new__(ZZ)
		return self
	
	def __init__(self, *var, **options):
		self.options = dict()
		if "mod" in options:
			self.mod = self.options["mod"] = options["mod"]
		else:
			self.mod = 0

		if "rel" in options and not isinstance(options["rel"], int):
			self.rel = self.options["rel"] = relation(options["rel"])
		else:
			self.rel = 0

		if len(var) > 0:
			self.const_vars = self.options["const_vars"] = var
		else:
			self.const_vars = tuple()

	def __repr__(self):
		repr_ = "%s" % self.__class__.__name__
		if self.options:
			repr_ += "("
			i = len(self.options)
			for k in self.options:
				repr_ += str(k) + " = " + str(self.options[k])
				i -= 1
				if i == 0:
					repr_ += ")"
				else:
					repr_ += ", "
		return repr_

	def __str__(self):
		return repr(self)

	@abstractmethod
	def number(self):
		raise NotImplementedError("This class is abstract class")

	@abstractmethod
	def has_relation(self):
		raise NotImplementedError("This class is abstract class")

def ring(*var, **options):
	return Ring(*var, **options)

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
	def __init__(self, *var, **options):

		"""
		Arugments:
		* var: indeterminate variables
		* options: allowed keys are "quo" and them of Rings'
		"""

		self.indet_vars = var

		if "coeff_dom" in options:
			validate_type(options["coeff_dom"], Ring)
			self.coeff_dom = options["coeff_dom"]
		else:
			self.coeff_dom = ring(**options)

		self.reduction_step = dict()

		if not self.coeff_dom.mod == 0:
			self.mod = self.coeff_dom.mod
			self.reduction_step["mod"] = True
		else:
			self.mod = 0
			self.reduction_step["mod"] = False

		if not self.coeff_dom.rel == 0:
			self.rel = relation(self.coeff_dom.rel)
			self.reduction_step["rel"] = True
		else:
			self.rel = 0
			self.reduction_step["rel"] = False

		self.const_vars = self.coeff_dom.const_vars

		if "quo" in options and not isinstance(options["quo"], int):
			self.quo = relation(options["quo"])
			self.reduction_step["quo"] = True
		else:
			self.quo = 0
			self.reduction_step["quo"] = False

	def __repr__(self):
		repr_ = "PolynomialRing("
		for v in self.indet_vars:
			repr_ += str(v) + ", "
		repr_ += "coeff: " + str(self.coeff_dom)
		if self.reduction_step["quo"]:
			repr_ += ", quo: " + str(self.quo)
		repr_ += ")"
		return repr_

	def reduce(self, rep):
		validate_type(rep, DP, int)
		if isinstance(rep, int):
			return rep
		else:
			reduce_ = rep
			if self.reduction_step["rel"]:
				reduce_ = _reduce(reduce_, self.rel)
			if self.reduction_step["quo"]:
				reduce_ = _reduce(reduce_, self.quo)
			if self.reduction_step["mod"]:
				reduce_ = reduce_ % self.mod
		return reduce_

def polynomialring(*var, **options):
	return PolynomialRing(*var, **options)

"""
* Reduction methods
Reduce variables with variable relations.
"""

def _reduce(f, relation):
	for r in relation[::-1]:
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