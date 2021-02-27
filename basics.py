class Symbol:
	def __init__(self, char):
		if not isinstance(char, str):
			raise TypeError("can generate Symbol only with str, not %s" % char.__class__.__name__)
		self.char = char

	def __repr__(self):
		return self.char

	def __str__(self):
		return self.char

	def __add__(a, b):
		return as_dup(a) + as_dup(b, a)

	def __radd__(a, b):
		return a + b

	def __sub__(a, b):
		return as_dup(a) - as_dup(b, a)

	def __rsub__(a, b):
		return - a + b

	def __mul__(a, b):
		return as_dup(a) * as_dup(b, a)
	
	def __rmul__(a, b):
		return a * b

	def __neg__(a):
		return a * -1
 
	def as_dup(self):
		return dup(self, [0, 1])

def symbol(c):
	if not isinstance(c, str):
		raise TypeError("symbol of variable must be str, not %s" % c.__class__.__name__)
	return Symbol(c)

def symbols(c):
	if not isinstance(c, str):
		raise TypeError("symbol of variable must be str, not %s" % c.__class__.__name__)
	symbols_ = c.replace(" ", "").split(",")
	return tuple([symbol(s) for s in symbols_])

class DUP:
	"""

	Represents a dense univariate polynomial.
	Coefficients need to be integer or DUP object.

	"""

	def __init__(self, var, coeffs=[]):
		if not isinstance(var, Symbol):
			raise TypeError("symbol must be Symbol, not %s" % symbol.__class__.__name__)
		elif not isinstance(coeffs, list):
			raise TypeError("coefficients must be list, not %s" % coeffs.__class__.__name__)
		else:
			self.var = var
			if coeffs == []:
				self._zero()
			else:
				self._set(coeffs)

	def _set(self, coeffs):
		self.is_zero = False
		self.coeffs = coeffs
		self.deg = len(coeffs) - 1
		self.trdeg = self.trailing_deg()
		self.inner_vars = self._inner_vars()

	def _zero(self):
		self.is_zero = True
		self.coeffs = [0]
		self.deg = "-oo"
		self.trdeg = "-oo"
		self.inner_vars = (self.var, )

	def __repr__(self):
		return self.sparse_rep()

	def __str__(self):
		return self.sparse_rep()

	"""
	Binary operations:
	add, sub, mul, flootdiv, mod, eq
	"""
	def __add__(f, g):
		return f._add(as_dup(g, f.var))

	def _add(f, g):
		if g.var == f.var:
			m = min(g.deg, f.deg) + 1
			coeffs_ = [f.coeffs[i] + g.coeffs[i] for i in range(m)]
			if f.deg > g.deg:
				coeffs_ += f.coeffs[m:]
			else:
				coeffs_ += g.coeffs[m:]
			return dup(f.var, coeffs_)
		else:
			add_ = bivariate_uniform(f, g)
			return add_[0] + add_[1]

	def __radd__(f, g):
		return f + g

	def __sub__(f, g):
		return f + g * (-1)

	def __rsub__(f, g):
	 	return - (f - g)

	def __mul__(f, g):
		return f._mul(as_dup(g, f.var)) 

	def _mul(f, g):
		if g.var == f.var:
			maxdeg = f.deg + g.deg
			coeffs_ = []
			for i in range(maxdeg + 1):
				co_ = 0
				for j in range(i + 1):
					if j < f.deg + 1 and i - j < g.deg + 1:
						co_ += f.coeffs[j] * g.coeffs[j - i]
					else:
						pass
				coeffs_.append(co_)
			return dup(f.var, coeffs_)
		else:
			mul_ = bivariate_uniform(f, g)
			return mul_[0] * mul_[1]

	def __rmul__(f, g):
		return f * g

	def __eq__(f, g):
		if isinstance(g, DUP):
			if f.var == g.var and f.coeffs == g.coeffs:
				return True
			else:
				return False
		elif isinstance(g, int):
			if f.deg == 0 and f.coeffs[0] == g:
				return True
			else:
				return False

	def __gt__(f, g):
		return True

	def __mod__(f, g):
		if isinstance(g, int):
			return dup(f.var, [c % g for c in f.coeffs]) 
		else:
			raise NotImplementedError()

	def __neg__(f):
		return f * -1

	def trailing_deg(self):
		for d in range(self.deg + 1):
			if self.coeffs[d] != 0:
				return d
		return None

	def sparse_rep(self):
		if self.is_zero:
			return "0"
		rep = ""
		deg = self.deg
		while deg >= 0:
			c = self.coeffs[deg]
			if isinstance(c, DUP) and c.deg > 0:
				if deg > 0:
					c = "(" + str(c) + ")"
				else:
					c = str(c)
			else:
				c = str(abs(c))
			if deg > 1:
				if c != "1":
					rep += c + "*"
				rep += str(self.var) + "**" + str(deg)
			elif deg == 1:
				if c != "1":
					rep += c + "*"
				rep += str(self.var)
			else:
				rep += c
			if deg != self.trdeg:
				while True:
					deg -= 1
					if self.coeffs[deg] != 0:
						break
				if self.coeffs[deg] > 0:
					rep += " + "
				else:
					rep += " - "
			else:
				break
		return rep

	def isconstant(self):
		if self.deg == 0:
			return True
		else:
			return False

	def degree(self, *var, total=False, as_dict=False, any_vars=False):
		"""
		Return the degree of self in the given variable.
		Given multiple variables, returns dict of each degree.
		If opetion "total" is True, returns total degree of it.

		Example:
		>>> a = symbols("a, b, c")
		>>> p = dup(a, [dup(b, [2, 0, dup(c, [1, 1])]), dup(b, [0, 3, dup(c, [1, 2, 1, 2])])])
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
				try:
					deg_ += self._degree(v)
				except TypeError:
					i += 1
			if i == len(self.inner_vars):
				return "-oo"
			else:
				return deg_
		elif as_dict:
			degs_ = dict()
			if any_vars:
				var = self.inner_vars
			elif len(var) == 0:
				var = (self.var, )
			for v in var:
				degs_[v] = self._degree(v)
			return degs_
		else:
			if len(var) == 0:
				return self.deg
			elif len(var) == 1:
				return self._degree(var[0])
			else:
				degs_ = tuple()
				for v in var:
					degs_ += (self._degree(v), )
				return degs_

	def _degree(self, var):
		if not var in self.inner_vars:
			raise ValueError("does not have the variable %s" % str(var))
		if var == self.var:
			return self.deg
		deg_ = -1
		for c in self.coeffs:
			if isinstance(c, DUP):
				if c.var == var:
					try:
						if c.deg > deg_:
							deg_ = c.deg
					except TypeError:
						continue
				else:
					deg_ = c._degree(var)
			else:
				continue
		if deg_ == -1:
			return "-oo"
		else:
			return deg_

	def as_list(self):
		coeffs_ = []
		for c in self.coeffs:
			if isinstance(c, DUP):
				coeffs_.append(c.as_list())
			else:
				coeffs_.append(c)
		return coeffs_

	def _inner_vars(self):
		vars_ = (self.var, )
		for c in self.coeffs:
			if isinstance(c, DUP):
				vars_ += c._inner_vars()
				break
			else:
				continue
		return vars_

	def is_univariate(self):
		for c in self.coeffs:
			if isinstance(c, DUP):
				return False
		return True

	def subs(self, subsdict):
		if not isinstance(subsdict, dict):
			raise TypeError("argument must be dict, not %s" % subsdict.__class__.__name__)
		if not self.var in subsdict:
			return self
		else:
			subs_ = 0
			value_ = subsdict[self.var]
			for i in range(self.deg + 1):
				subs_ += self.coeffs[i] * value_ ** i
			return subs_

	def variables_sort(self, var):
		"""
		Sort variables to given order.
		
		Example:
		>>> a, b, c = symbols("a, b, c")
		>>> p = dup(a, [dup(b, [2, 0, dup(c, [1, 1])]), dup(b, [0, 3, dup(c, [1, 2, 1, 2])])])
		>>> p
		((2*c**3 + c**2 + 2*c + 1)*b**2 + 3*b)*a + (c + 1)*b**2 + 2
		>>> p.variables_sort((b, a, c))

		"""
		if len(self.inner_vars) == 1:
			raise TypeError("can not sort the variable of univariate polynomials")
		if not set(self.inner_vars) == set(var):
			raise ValueError("variables must be the same as former one")
		pre, post = self.as_list(), []
		for c in self.coeffs:
			pass

def dup(symbol, coeffs):
	return DUP(symbol, coeffs)

def dup_from_dict(var, rep):
	if not isinstance(var, tuple):
		raise TypeError("first argument must be tuple, not %s" % var.__class__.__name__)
	if len(var) == 0:
		raise ValueError("needs one variable at least")
	if len(var) == 1:
		coeffs_, i, j = [], 0, 0
		while i < len(rep):
			try:
				coeffs_.append(rep[j])
				i += 1
			except KeyError:
				coeffs_.append(0)
			finally:
				j += 1
		return dup(var[0], coeffs_)
	else:
		pass

def dup_from_list(var, rep):
	if not isinstance(var, tuple):
		raise TypeError("first argument must be tuple, not %s" % var.__class__.__name__)
	if len(var) == 0:
		raise ValueError("needs one variable at least")
	if len(var) == 1:
		return dup(var[0], rep)
	else:
		coeffs_ = []
		for r in rep:
			if isinstance(r, list):
				coeffs_.append(dup_from_list(var[1:], r))
			elif isinstance(r, int):
				coeffs_.append(r)
			else:
				raise ValueError("elements of coeffs list must be int or DUP, not %s" % r.__class__.__name__)
		return dup(var[0], coeffs_)

def as_dup(f, *symbol):
	if isinstance(f, int):
		if len(symbol) == 0:
			raise ValueError("needs one symbol when f is int")
		else:
			return dup(symbol[0], [f])
	elif isinstance(f, Symbol):
		return f.as_dup()
	elif isinstance(f, DUP):
		return f
	elif isinstance(f, DMP):
		return f.rep
	else:
		raise TypeError("argument must be int, Symbol, DUP or DMP, not %s" % f.__class__.__name__)

class DMP:
	"""

	Represents a dense multivariate polynomial.

	Member:
	* symbols: tuple of Symbol objects
	* rep: DUP object (recursive polynomial)

	"""

	def __init__(self, symbols, rep):
		if isinstance(symbols, Symbol):
			symbols = (symbols, )
		elif isinstance(symbols, list):
			symbols = tuple(symbols)
		elif isinstance(symbols, tuple):
			pass
		else:
			raise TypeError("symbols argument must be Symbol, list or tuple, not %s" % symbols.__class__.__name__)
		self.symbols = symbols
		if isinstance(rep, DUP):
			self.rep = rep
		elif isinstance(rep, list):
			self.rep = dup_from_list(symbols, rep)
		elif isinstance(rep, dict):
			raise NotImplementedError

	def __repr__(self):
		return self.rep.sparse_rep()

	def __str__(self):
		return self.rep.sparse_rep()

	def __add__(f, g):
		if isinstance(g, DMP):
			if f.same_variables(g):
				if f.symbols == g.symbols:
					return dmp(f.symbols, f.rep + g.rep)
				else:
					# uniform variables with variable_sort()
					raise NotImplementedError()
			else:
				# uniform variables
				pass
		else:
			pass

	def __radd__(f, g):
		return f + g

	def __sub__(f, g):
		pass

	def __rsub__(f, g):
		return - f + g

	def __mul__(f, g):
		if isinstance(g, DMP):
			if f.same_variables(g):
				if f.symbols == g.symbols:
					return dmp(f.symbols, f.rep * g.rep)
				else:
					raise NotImplementedError()
			else:
				# uniform variables
				pass
		else:
			pass

	def variables_sort(self, t):
		if not isinstance(t, tuple):
			raise TypeError("needs one tuple not %s", t.__class__.__name__)
		if not set(self.symbols) == set(t):
			raise ValueError("can sort only if both polynomial have the same variables")
		if self.symbols == t:
			raise ValueError("already sorted")
		return dmp(t, self.rep.variables_sort(t))

	def degree(self, *var):
		pass

	def subs(self, subsdict):
		pass

	def same_variables(self, other):
		if set(self.symbols) == set(other.symbols):
			return True
		else:
			return False

	def as_list(self):
		return self.rep.as_list()

	def as_dup(self):
		return self.rep

def dmp(symbols, rep):
	return DMP(symbols, rep)

def abs(rep):
	if isinstance(rep, int):
		if rep >= 0:
			return rep
		else:
			return - rep
	elif isinstance(rep, DUP):
		if rep.deg > 0:
			return rep
		else:
			if rep.coeffs[0] > 0:
				return rep
			else:
				return - rep

def bivariate_uniform(f, g):
	"""

	Transform two univariate polynomials to have the same two variables.
	Given f(x) and g(y), it returns f'(x,y) and g'(x,y)
	with f'(x, 0) = f(x), f'(0, y) = 0, g'(0, y) = g(y) and g'(x, 0) = 0

	Example:
	>>> x, y = symbols('x, y')
	>>> f = dup(x, [1, 1]) # f = x + 1
	>>> g = dup(y, [1, 1]) # g = y + 1
	>>> bivariate_uniform(f, g)
	(DMP((x, y), DMP(x, [DUP(y, [1]), 1])), DMP((x,y), DUP(x, [DUP(y, [1, 1]), 1])))

	"""
	symbols_ = (f.symbol, g.symbol)
	f_list = [dup(g.symbol, [f.coeffs[0]])] + f.coeffs[1:]
	f_ = dmp(symbols_, dup(f.symbol, f_list))
	g_ = dmp(symbols_, dup(f.symbol, [g]))
	return (f_, g_)