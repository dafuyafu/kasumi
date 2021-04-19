from pys.mathtools import is_prime
from pys.pytools import tuple_union, tuple_or_object, validate_type
from basics.basictools import symbols, Symbol, dp, DP, monomial_from_index
from abc import ABCMeta, abstractmethod
import itertools
import random

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
			str_, flag_first = "(", True
			for rel in self:
				if flag_first:
					flag_first = False
				else:
					str_ += ", "
				str_ += str(rel["rep"])
			return str_ + ")"
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

	def get_variable(self):
		var_ = list()
		for r in self:
			var_.append(r["var"])
		return tuple(var_)

	def as_tuple(self):
		tuple_ = list()
		for r in self:
			tuple_.append(r["rep"])
		return tuple(tuple_)

def relation(reps, mod=0):
	return Relation(reps, mod)

"""
	* Reduction methods
	Reduce variables with variable relations.
"""

def reduce_relation(f, relation):
	validate_type(relation, Relation)
	return _reduce(f, relation)

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
		FiniteField(mod: 3, rel: a**2 - 2, const_vars: (a, ))

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

					""" 
						* Require Updates !! 
						Given variables without relation, must return infinite objects.

					"""
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
			var_tuple_ = tuple(r["var"] for r in self.rel)
			self.const_vars = self.options["const_vars"] = tuple_union(var, var_tuple_)
		else:
			self.rel = 0
			if len(var) == 0:
				self.const_vars = tuple()
			else:
				self.const_vars = self.options["const_vars"] = var

	def __repr__(self):
		repr_ = "%s" % self.__class__.__name__
		if self.options:
			repr_ += "("
			i = len(self.options)
			for k in self.options:
				repr_ += str(k) + ": " + str(self.options[k])
				i -= 1
				if i == 0:
					repr_ += ")"
				else:
					repr_ += ", "
		return repr_

	def __str__(self):
		return repr(self)

	def __len__(self):
		if self.is_finite:
			return self.number()
		else:
			raise TypeError("this is an infinite domain")

	@abstractmethod
	def number(self):
		raise NotImplementedError("This class is abstract class")

	@abstractmethod
	def has_relation(self):
		raise NotImplementedError("This class is abstract class")

	def ex_degree(self):
		raise TypeError("can calculate the extension degree only of finite domains")

	"""
	* Iterators
	"""

	def it_elements(self, infinite=False):
		if infinite:
			pass
		else:
			if self.mod == 0:
				raise TypeError("cannot generate elements of infinite domain")
		if self.const_vars:
			n = self.number()
			for i in itertools.count():
				if i == n:
					break
				yield self.element(i)
		else:
			for i in itertools.count():
				if i == self.mod:
					break
				yield i

	def element(self, *n):
		if len(n) > 1:
			return tuple(self.element(i) for i in n)
		elif len(n) == 1:
			i = n[0]
		else:
			i = random.randrange(self.number())
		if i > len(self) - 1:
			raise ValueError("index argument must be equal to or less than the number of domain elements")
		if self.const_vars:
			coeff_list = list()
			exdeg = self.ex_degree()
			for d in range(exdeg)[::-1]:
				c = i // (self.mod ** d)
				if c > self.mod / 2:
					c -= self.mod
				coeff_list.append(c)
				i %= self.mod ** d
			rep_ = i = 0
			for e in itertools.product(*[range(r["deg"]) for r in self.rel]):
				index = (e, coeff_list[exdeg - i - 1])
				rep_ += monomial_from_index(index, self.const_vars)
				i += 1
			return rep_
		else:
			return i

	def it_points(self, *var, infinite=False):
		if infinite:
			pass
		else:
			if self.mod == 0:
				raise TypeError("cannot generate elements of infinite domain")
		for e in itertools.product(self.it_elements(), repeat=len(var)):
			yield dict(zip(var, e))

	def extend(self, var, rel=0):
		validate_type(var, Symbol)
		if var in self.const_vars:
			raise ValueError("it already has the variable '%s'" % str(var))
		validate_type(rel, int, DP)
		const_vars_ = self.const_vars + (var, )
		if rel == 0:
			return ring(*const_vars_, mod=self.mod, rel=self.rel)
		else:
			if self.rel == 0:
				return ring(*const_vars_, mod=self.mod, rel=rel)
			else:
				dps_ = self.rel.as_tuple() + (rel, )
				return ring(*const_vars_, mod=self.mod, rel=dps_)

	def is_prime(self):

		"""

			* is_prime()
			validate if self is a prime field

		"""

		if self.is_field:
			if self.mod == 0:
				if self.__class__ == QQ:
					return True
				else:
					return False
			else:
				if self.rel == 0:
					return True
				else:
					return False
		else:
			return False

	def reduce(self, rep):
		validate_type(rep, Symbol, DP, int)
		if isinstance(rep, int):
			if self.mod > 0:
				rep %= self.mod
			if rep > self.mod // 2:
				rep -= self.mod
			return rep
		else:
			rep = rep.as_dp()
			if self.rel != 0:
				rep = reduce_relation(rep, self.rel)
			if self.mod > 0:
				rep %= self.mod
		return rep

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
		return self.mod ** self.ex_degree()

	def has_relation(self):
		return "rel" in self.options

	def ex_degree(self):
		if self.has_relation():
			ex_degree_ = 1
			for r in self.rel:
				ex_degree_ *= r['deg']
			return ex_degree_
		else:
			return 1

class FiniteField(Field, FiniteRing):
	"""

	* FiniteField

	represent finite fields such as ZZ/(5) or ZZ(5)[a]/(a**2 - 2)

	"""
	pass

class ACFiniteField(Field):
	is_alg_closed = True