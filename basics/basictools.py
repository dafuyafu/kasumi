from pys.pytools import validate_type, tuple_union, tuple_minus

from multiprocessing import Pool
import os
import itertools
import numpy as np

class Symbol:
	def __init__(self, char):
		validate_type(char, str)
		self.char = char
		self.inner_vars = (self, )

	def __repr__(self):
		return self.char

	def __str__(self):
		return self.char

	def __add__(a, b):
		return as_dp(a) + as_dp(b, a)

	def __radd__(a, b):
		return a + b

	def __sub__(a, b):
		return as_dp(a) - as_dp(b, a)

	def __rsub__(a, b):
		return - a + b

	def __mul__(a, b):
		return as_dp(a) * as_dp(b, a)
	
	def __rmul__(a, b):
		return a * b

	def __pow__(a, b):
		return as_dp(a) ** b

	def __neg__(a):
		return a * -1

	def __eq__(a,b):
		return a.char == b.char

	def __hash__(a):
		return hash(a.char)
 
	def as_dp(self):
		return dp(self, [0, 1])

def symbol(c):
	validate_type(c, str)
	return Symbol(c)

def symbols(c):
	validate_type(c, str)
	symbols_ = c.replace(" ", "").split(",")
	if len(symbols_) == 1:
		return symbol(symbols_[0])
	else:
		return tuple([symbol(s) for s in symbols_])

def it_symbols(c, limit=0):
	for i in itertools.count(1):
		if i == limit:
			break
		c_ = c + "_" + str(i)
		yield symbol(c_)

class DP:
	"""

	Represents a dense polynomial.
	Coefficients need to be integer or DP object.

	"""

	def __init__(self, symbol, coeffs):
		if not isinstance(symbol, Symbol):
			raise TypeError("symbol must be Symbol, not %s" % symbol.__class__.__name__)
		if isinstance(coeffs, list):
			coeffs = tuple(coeffs)
		elif isinstance(coeffs, tuple):
			pass
		else:
			raise TypeError("coefficients must be tuple or list, not %s" % coeffs.__class__.__name__)
		self.var = symbol
		self._set_init(coeffs)
		if self.is_zero():
			self._set_zero(coeffs)

	"""
	* Initialize methods
	Use only to initialize instances.
	"""

	def _set_init(self, coeffs):	
		self.coeffs = coeffs
		self.deg = self._set_deg()
		self.trdeg = self._trailing_deg()
		self.inner_vars = self._inner_vars()

	def _set_deg(self):
		deg_ = -1
		for i in range(len(self)):
			if isinstance(self[i], int):
				if self[i] != 0:
					deg_ = i
			else:
				if not self[i].is_zero():
					deg_ = i
		return deg_

	def _set_zero(self, coeffs):
		self.coeffs = coeffs
		self.deg = -1
		self.trdeg = -1
		self.inner_vars = self._inner_vars()

	def zero_reduce(self):
		self_ = self
		self_._zero_reduce()
		return self_

	def _zero_reduce(self):
		for i in range(len(self)):
			if i == 0:
				continue
			else:
				if isinstance(self[i], DP) and self[i].is_zero:
					self[i] = 0

	def _inner_vars(self):
		vars_ = (self.var, )
		inner_ = tuple()
		for c in self:
			if isinstance(c, DP):
				inner_ = tuple_union(inner_, c.inner_vars)
			else:
				continue
		return vars_ + inner_

	def _trailing_deg(self):
		for d in range(self.deg + 1):
			if isinstance(self[d], int) and self[d] != 0:
				return d
			elif isinstance(self[d], DP) and not self[d].is_zero():
				return d
		return None

	"""
	* Representation magic methods
	"""

	def __repr__(self):
		if len(self.inner_vars) > 1:
			repr_ = "DP(" + self.as_dist() + ", " + str(self.inner_vars)
		else:
			repr_ = "DP(" + self.as_dist() + ", " + str(self.var)
		return repr_ + ")"

	def __str__(self):
		return self.as_dist()

	"""
	Arithmetic operations:
	add, sub, mul, flootdiv, mod, eq and so on.
	"""
	def __add__(f, g):
		return f._add(as_dp(g, f.var))

	def _add(f, g):
		if g.inner_vars == f.inner_vars:
			m = min(g.degree(non_negative=True), f.degree(non_negative=True)) + 1
			coeffs_ = [f[i] + g[i] for i in range(m)]
			if f.deg > g.deg:
				coeffs_ += f[m:]
			else:
				coeffs_ += g[m:]
			return dp(f.var, coeffs_)
		else:
			vars_ = tuple_union(f.inner_vars, g.inner_vars)
			return f.sort_vars(vars_)._add(g.sort_vars(vars_))

	def __radd__(f, g):
		return f + g

	def __sub__(f, g):
		return f + g * (-1)

	def __rsub__(f, g):
	 	return - f + g

	def __mul__(f, g):
		return f._mul(as_dp(g, f.var)) 

	def _mul(f, g):
		if g.inner_vars == f.inner_vars:
			f_d, g_d = f.degree(non_negative=True), g.degree(non_negative=True)
			maxdeg = f_d + g_d
			coeffs_ = []
			for i in range(maxdeg + 1):
				co_ = 0
				for j in range(i + 1):
					if j < f_d + 1 and i - j < g_d + 1:
						co_ += f[j] * g[i - j]
					else:
						pass
				coeffs_.append(co_)
			return dp(f.var, coeffs_)
		else:
			vars_ = tuple_union(f.inner_vars, g.inner_vars)
			return f.sort_vars(vars_)._mul(g.sort_vars(vars_))

	def _mur_ka(f, g):

		"""
		multiple f and g with Karatsuba Algorithm
		"""

		f_d, g_d = f.degree(non_negative=True), g.degree(non_negative=True)
		if f_d > g_d:
			f_, g_ = f.coeffs, g.coeffs + (0, )*(f_d - g_d)
		else:
			f_, g_ = f.coeffs + (0, )*(g_d - f_d), g.coeffs
		list_ = [f_[i] * g_[i] for i in range(max(f_d, g_d) + 1)]
		dict_ = dict()
		for s in range(len(f_)):
			for t in range(s + 1, len(list_)):
				dict_[(s,t)] = (f[s] + f[t]) * (g[s] + g[t])
		coeffs_ = list()
		for i in range(max(f_d, g_d) * 2 + 1):
			if i == 0:
				coeffs_.append(list_[0])
			elif i == f_d + g_d:
				coeffs_.append(list_[-1])
			else:
				co_former = 0
				for s in range(i):
					if s >= i - s:
						break
					try:
						co_former += dict_[(s, i - s)]
					except KeyError:
						pass
				co_latter = 0
				for s in range(i):
					if s >= i - s:
						break
					if s >= len(list_) or i - s >= len(list_):
						continue
					co_latter += list_[s] + list_[i - s]
				if i % 2 == 0:
					coeffs_.append(co_former - co_latter + list_[i // 2])
				else:
					coeffs_.append(co_former - co_latter)
		print(coeffs_)
		return dp(f.var, coeffs_)


	def __rmul__(f, g):
		return f * g

	def __pow__(f, e):
		if not isinstance(e, int):
			raise TypeError("second operand must be int, not %s" % e.__class__.__name__)
		if e < 0:
			raise ValueError("exponent must be positive")
		if e == 0:
			return dp(f.var, [1]).sort_vars(f.inner_vars)
		num_ = bin(e).replace('0b', '')
		pow_ = 1
		for i in num_:
			pow_ = pow_ * pow_
			if i == '1':
				pow_ = f * pow_
		return pow_

	def __truediv__(f, g):
		raise TypeError("unsupported operand type(s) for '/'")

	def __floordiv__(f, g):
		if isinstance(g, int):
			coeffs_ = [c // g for c in f.coeffs]
			return dp(f.var, coeffs_)
		else:
			pass

	def div(f, g, mod):
		if isinstance(g, int):
			coeffs_ = list()
			for c in f:
				if isinstance(c, int):
					coeffs_.append(c * pow(g, mod - 2, mod))
				else:
					coeffs_.append(c.div(g, mod))
			return dp(f.var, coeffs_)
		else:
			raise TypeError("unsupported operand type(s) for div()")

	def __eq__(f, g):
		if isinstance(g, DP):
			if f.inner_vars == g.inner_vars and f.coeffs == g.coeffs:
				return True
			else:
				return False
		elif isinstance(g, int):
			if f.is_constant():
				if f[0] == g:
					return True
				else:
					return False
			else:
				return False
		else:
			raise TypeError("unsupported operand type(s) for '='")

	def __ne__(f, g):
		return not f == g

	def __gt__(f, g):
		if isinstance(g, int):
			return f[-1] > g

	def __lt__(f, g):
		return f != g and not f > g

	def __ge__(f, g):
		return f == g or f > g

	def __le__(f, g):
		return not f > g

	def __mod__(f, g):
		if isinstance(g, int):
			coeffs_ = list()
			for c in f:
				if isinstance(c, int) and c % g > g // 2:
					coeffs_.append(c % g - g)
				else:
					coeffs_.append(c % g)
			while len(coeffs_) > 1:
				if coeffs_[-1] == 0:
					coeffs_.pop()
				else:
					break
			return dp(f.var, coeffs_)
		elif isinstance(g, DP):
			if f.is_univariate() and f.var == g.var:
				return f - (f // g) * g
			else:
				raise TypeError("both of operands must have the same one variate")
		else:
			raise TypeError("unsupported operand type(s) for '%'")

	def __neg__(f):
		return f * -1

	def __abs__(self):
		if self.deg == 0:
			if self[0] > 0:
				return self
			else:
				return - self
		else:
			return self

	def __bool__(self):
		if self.is_zero():
			return False
		else:
			return True

	def __hash__(self):
		return hash(self.as_tuple())

	"""
	* Iterable magic methods
	"""

	def __getitem__(self, key):
		if isinstance(key, (int, slice)):
			try:
				return self.coeffs[key]
			except IndexError:
				return 0
		elif isinstance(key, dict):
			list_ = list()
			for v in self.inner_vars:
				list_.append(key[v])
			while len(list_) < len(self.inner_vars):
				list_.append(0)
			return self[list_]
		elif isinstance(key, (list, tuple)):
			if len(key) == 1:
				return self[key[0]]
			else:
				try:
					return self[key[0]][key[1:]]
				except TypeError:
					if any(key[1:]):
						return 0
					else:
						return self[key[0]]
		else:
			raise TypeError("key must be int, slice, dict, tuple or list, not '%s'" % key.__class__.__name__)

	def get_constant(self):
		return self[(0, )*len(self.inner_vars)]

	def __setitem__(self, key, value):
		if isinstance(key, (int, slice)):
			coeffs_ = list(self.coeffs)
			try:
				coeffs_[key] = value
			except IndexError:
				for i in range(key - len(self)):
					coeffs_.append(0)
				coeffs_.append(value)
			self._set_init(tuple(coeffs_))
		elif isinstance(key, dict):
			key_ = list()
			for v in self.inner_vars:
				key_.append(key[v])
			while len(key_) < len(self.inner_vars):
				key_.append(0)
			self._set(key_, value, self.inner_vars)
		elif isinstance(key, (list, tuple)):
			if len(key) > len(self.inner_vars):
				raise ValueError("length of key tuple(or list) must be less than or equal to the one of self.inner_vars")
			elif len(key) < len(self.inner_vars):
				key = list(key)
				for i in range(len(self.inner_vars) - len(key)):
					key.append(0)
			self._set(key, value, self.inner_vars)
		else:
			raise TypeError("key must be int, slice, dict, list or tuple, not '%s" % key.__class__.__name__)

	def _set(self, var, value, inner_vars):
		if len(var) == 1:
			self[var[0]] = value
		else:
			if isinstance(self[var[0]], DP):
				return self[var[0]]._set(var[1:], value, inner_vars[1:])
			else:
				if not any(var[1:]):
					self[var[0]] = value
				else:
					dp_ = dp(inner_vars[1], [self[var[0]]])
					dp_._set(var[1:], value, inner_vars[1:])
					self[var[0]] = dp_

	def __iter__(self):
		return iter(self.coeffs)

	def __len__(self):
		return len(self.coeffs)	

	"""
	* Iterator methods
	"""

	def it_dist(self, termorder="lex", with_index=False):
		validate_type(termorder, str)
		if termorder == "lex":
			iters = [range(self.degree(v, non_negative=True), -1, -1) for v in self.inner_vars]
			for p in itertools.product(*iters):
				c = self[p]
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

	def it_reversed_dist(self, index=None, with_index=False):
		i = 0
		if index is None:
			index = dict()
		for c in self:
			index_ = index
			if isinstance(c, int):
				if c == 0:
					pass
				else:
					if with_index:
						for v in self.inner_vars:
							if v == self.var:
								index[v] = i
							else:
								index[v] = 0
						yield (index, c)
					else:
						yield c
			else:
				index[self.var] = i
				yield from c.it_reversed_dist(index, with_index)
			index = index_
			i += 1

	"""
	* Representation Methods
	Default method is 'as_dist()'
	"""
	def as_rec_rep(self):
		if self.is_zero():
			return "0"
		rep_ = ""
		deg = self.deg # parse from higher degree
		while deg >= 0:
			c = self[deg]

			"""
			* coefficient part
			translate int or DP object to str
			"""

			if isinstance(c, DP):
				if c.is_zero():
					deg -= 1
					break
				elif c.is_monomial():
					c = str(c)
				else:
					if deg > 0:
						c = "(" + str(c) + ")"
					else:
						c = str(c)
			else:
				if deg == self.deg:
					c = str(c)
				else:
					c = str(abs(c))

			"""
			* variable and exponent part
			connect coeff, variable and exponent index
			"""

			if deg > 1:
				if c == "1":
					pass
				elif c == "-1":
					rep_ += "- "
				else:
					rep_ += c + "*"
				rep_ += str(self.var) + "**" + str(deg)
			elif deg == 1:
				if c == "1":
					pass
				elif c == "-1":
					rep_ += "- "
				else:
					rep_ += c + "*"
				rep_ += str(self.var)
			else:
				rep_ += c

			"""
			* plus or minus part
			output + (repr. -) if the next coeff which is not zero is positive (repr. negative)
			"""

			if deg != self.trdeg:
				while True:
					deg -= 1
					if self[deg] != 0:
						break
				if isinstance(self.coeffs[deg], int):
					if self[deg] > 0:
						rep_ += " + "
					else:
						rep_ += " - "
				else:
					rep_ += " + "
			else:
				break
		return rep_

	def as_dict(self, l=-1, t=tuple(), d=dict()):
		"""

		Return self as dict type.
		Note that this is also dense representation
		because return 0 if given keys which it does not have.

		* Example:

		>>> f = a + b + a*b
		>>> f.as_dict()
		{(1,0): 1, (0,1): 1, (1,1): 1}

		"""
		if l == -1:
			l = len(self.inner_vars)
		for i in range(len(self)):
			t += (i, )
			if isinstance(self[i], int):
				d[t + zero_tuple(l - 1)] = self[i]
			else:
				d.update(self[i].as_dict(l - 1, t, d))
			t = t[:(len(t) - 1)]
		return d

	def as_list(self):
		list_ = list()
		for c in self:
			if isinstance(c, int):
				list_.append(c)
			else:
				list_.append(c.as_list())
		return list_

	def as_tuple(self):
		list_ = list()
		for c in self:
			if isinstance(c, int):
				list_.append(c)
			else:
				list_.append(c.as_tuple())
		return tuple(list_)

	def as_dist_tuple(self, termorder="lex", with_index=False):
		validate_type(termorder, str)
		if termorder == "lex":
			iters, list_ = list(), list()
			for v in self.inner_vars[::-1]:
				iters.append(range(self.degree(v) + 1)[::-1])
			for p in itertools.product(*iters):
				if with_index:
					list_.append((p[::-1], self[p[::-1]]))
				else:
					list_.append(self[p[::-1]])
		elif termorder == "grevlex":
			raise NotImplementedError()
		else:
			raise TypeError("given unsupported termorder '%s'" % termorder)
		return tuple(list_)

	def as_dist(self, termorder="lex"):
		"""
		Return distributed representation.
		"""

		rep_, lt = str(), True
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

			# coeffcients part
			if data[1] == 1 or data[1] == -1:
				flag_ = True
			else:
				flag_ = False
				rep_ += str(abs(data[1]))
				if any(data[0]):
					rep_ += "*"

			# variables part
			for i in range(len(self.inner_vars)):
				if data[0][i] > 1:
					rep_ += str(self.inner_vars[i]) + "**" + str(data[0][i])
				elif data[0][i] == 1:
					rep_ += str(self.inner_vars[i])
				else:
					continue
				if i < len(self.inner_vars) - 1 and any(data[0][(i + 1):]):
					rep_ +="*"

			# last one modification
			if flag_ and not any(data[0]):
				rep_ += "1"

		return rep_

	def as_dp(self):
		return self

	"""
	* Validators
	These functions return bool.
	"""

	def is_constant(self):
		if self.degree(total=True, non_negative=True) == 0: 
			return True
		else:
			return False

	def is_monomial(self):
		if self.is_zero() or (self.is_univariate() and self.deg == self.trdeg):
			return True
		else:
			return False

	def is_univariate(self):
		for c in self:
			if isinstance(c, DP):
				if c.is_constant():
					continue
				return False
		return True

	def is_zero(self):
		for c in self:
			if isinstance(c, int):
				if not c == 0:
					return False
			else:
				if not c.is_zero():
					return False
		return True

	def is_monic(self, termorder="lex"):
		if self.LC(termorder) == 1:
			return True
		else:
			return False

	"""
	* Degree function and its subfunctions
	"""

	def degree(self, *var, total=False, as_dict=False, any_vars=False, non_negative=False):
		"""
		Return the degree of self in the given variable.
		Given multiple variables, returns dict of each degree.
		If opetion "total" is True, returns total degree of it.

		Example:
		>>> a = symbols("a, b, c")
		>>> p = dp(a, [dp(b, [2, 0, dp(c, [1, 1])]), dp(b, [0, 3, dp(c, [1, 2, 1, 2])])])
		>>> p
		((2*c**3 + c**2 + 2*c + 1)*b**2 + 3*b)*a + (c + 1)*b**2 + 2
		>>> p.degree()
		1
		>>> p.degree(a)
		1
		>>> p.degree(a, b)
		(1, 2)
		>>> p.degree(any_var)
		(1, 2, 3)
		>>> p.degree(total=True)
		6
		>>> p.degree(a, as_dict=True)
		{a: 1}
		>>> p.degree(any_var=True, as_dict=True)
		{a: 1, b: 2, c: 3}
		"""
		if total:
			deg_, i = 0, 0
			for v in self.inner_vars:
				d_ = self._degree(v)
				if d_ < 0:
					i += 1
				else:
					deg_ += d_
			if i == len(self.inner_vars):
				if non_negative:
					return 0
				else:
					return -1
			else:
				return deg_
		elif as_dict:
			degs_ = dict()
			if any_vars:
				var = self.inner_vars
			elif len(var) == 0:
				var = (self.var, )
			for v in var:
				degs_[v] = self._degree(v, non_negative)
			return degs_
		else:
			if len(var) == 0:
				if non_negative:
					return max(0, self.deg)
				else:
					return self.deg
			elif len(var) == 1:
				return self._degree(var[0], non_negative)
			else:
				return tuple([self._degree(v, non_negative) for v in var])

	def _degree(self, var, non_negative=False):
		if self.var == var:
			if non_negative:
				return max(0, self.deg)
			else:
				return self.deg
		else:
			deg_ = -1
			for c in self:
				if isinstance(c, DP):
					if c.var == var:
						if c.deg > deg_:
							deg_ = c.deg
					else:
						deg_ = c._degree(var, non_negative)
				else:
					continue
			if non_negative:
				return max(0, deg_)
			else:
				return deg_

	def subs(self, subsdict):
		validate_type(subsdict, dict)
		inner_vars_values = tuple_union(*[v.inner_vars for v in subsdict.values() if isinstance(v, (DP, Symbol))])
		inner_vars_ = tuple_union(self.inner_vars, inner_vars_values)
		f_ = self.sort_vars(inner_vars_)
		for k, v in subsdict.items():
			if isinstance(k, Symbol):
				k = as_dp(k).sort_vars(inner_vars_)
			elif isinstance(k, DP):
				if k.is_monomial():
					if k.is_monic():
						k = k.sort_vars(inner_vars_)
					else:
						raise ValueError("keys of subsdict must be monic")
				else:
					raise ValueError("keys of subsdict must be monomial")
			else:
				raise TypeError("keys of subsdict must be Symbol or DP, not '%s'" % k.__class__.__name__)
			f_ = f_._subs(k, v)
		return f_

	def _subs(self, key, value):
		key_index = np.array(next(key.it_dist(with_index=True))[0])
		f_ = self
		for mon in self.it_dist(with_index=True):
			mon_index, i = np.array(mon[0]), 0
			while True:
				mon_index -= key_index
				if np.all(mon_index >= 0):
					i += 1
				else:
					mon_index += key_index
					break
			if i > 0:
				mon_ = monomial_from_index(mon, self.inner_vars)
				val_ = monomial_from_index((mon_index, mon[1]), self.inner_vars) * value ** i
				f_ = f_ - mon_ + val_
			else:
				continue
		return f_

	"""
	* Variable sorting methods
	"""

	def sort_vars(self, vars_tuple):
		dp_ = dp_from_int(0, vars_tuple)
		for c in self.it_reversed_dist(with_index=True):
			vars_dict = c[0]
			for v in tuple_minus(vars_tuple, self.inner_vars):
				vars_dict[v] = 0
			dp_[vars_dict] = c[1]
		return dp_

	"""
	* Singularity methods
	"""

	def diff(self):
		coeffs_ = list()
		for i in range(len(self) - 1):
			coeffs_.append((i + 1) * self[i + 1])
		return dp(self.var, coeffs_)

	"""
	* Groebner basis methods
	"""

	def LT(self, termorder="lex", as_data=False):
		for d in self.it_dist(termorder=termorder, with_index=True):
			if d[1] != 0:
				if as_data:
					return d
				else:
					return monomial_from_index(d, self.inner_vars)
			else:
				continue
		return monomial_from_index(((0, )*len(self.inner_vars), 0), self.inner_vars)

	def LM(self, termorder="lex"):
		for d in self.it_dist(termorder=termorder, with_index=True):
			if d[1] != 0:
				return monomial_from_index((d[0], 1), self.inner_vars)
			else:
				continue
		return monomial_from_index(((0, )*len(self.inner_vars), 0), self.inner_vars)

	def LC(self, termorder="lex"):
		for d in self.it_dist(termorder=termorder):
			if d != 0:
				return d
			else:
				continue
		return 0

def dp(symbol, coeffs):
	return DP(symbol, coeffs)

def dp_from_dict(var, coeffs):
	validate_type(coeffs, dict)
	if len(var) == 0:
		raise ValueError("needs one variable at least")
	if len(var) == 1:
		coeffs_, i, j = [], 0, 0
		while i < len(coeffs):
			try:
				coeffs_.append(coeffs[j])
				i += 1
			except KeyError:
				coeffs_.append(0)
			finally:
				j += 1
		return dp(var[0], coeffs_)
	else:
		pass

def dp_from_list(var, coeffs):
	validate_type(coeffs, tuple, list)
	if len(var) == 0:
		raise ValueError("needs one variable at least")
	if len(var) == 1:
		return dp(var[0], coeffs)
	else:
		coeffs_ = []
		for r in coeffs:
			if isinstance(r, list):
				coeffs_.append(dp_from_list(var[1:], r))
			elif isinstance(r, int):
				coeffs_.append(r)
			else:
				raise ValueError("elements of coeffs list must be int or dp, not %s" % r.__class__.__name__)
		return dp(var[0], coeffs_)

def monomial_from_index(data, var):
	"""
	Return monomial from data of iterator.
	The argument data must be a tuple of exponents of variables and its coefficient.
	"""
	coeffs_ = [0]*data[0][0]
	if len(data[0]) == 1:
		coeffs_.append(data[1])
	else:
		data_ = (data[0][1:], data[1])
		coeffs_.append(monomial_from_index(data_, var[1:]))
	return dp(var[0], coeffs_)

def dp_from_int(n, var):
	validate_type(n, int)
	validate_type(var, (tuple, list))
	coeffs_ = [n]
	for v in var[1:][::-1]:
		coeffs_ = [dp(v, coeffs_)]
	return dp(var[0], coeffs_)

def as_dp(f, *symbol):
	validate_type(f, int, Symbol, DP)
	if isinstance(f, int):
		if len(symbol) == 0:
			raise ValueError("needs one symbol when f is int")
		else:
			return dp(symbol[0], [f])
	elif isinstance(f, Symbol):
		return f.as_dp()
	else:
		return f