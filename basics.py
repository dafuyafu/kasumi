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

	def __neg__(a):
		return a * -1
 
	def as_dp(self):
		return dp(self, [0, 1])

def symbol(c):
	if not isinstance(c, str):
		raise TypeError("symbol of variable must be str, not %s" % c.__class__.__name__)
	return Symbol(c)

def symbols(c):
	if not isinstance(c, str):
		raise TypeError("symbol of variable must be str, not %s" % c.__class__.__name__)
	symbols_ = c.replace(" ", "").split(",")
	return tuple([symbol(s) for s in symbols_])

class DP:
	"""

	Represents a dense polynomial.
	Coefficients need to be integer or DP object.

	"""

	def __init__(self, symbol, coeffs=[]):
		if not isinstance(symbol, Symbol):
			raise TypeError("symbol must be Symbol, not %s" % symbol.__class__.__name__)
		elif not isinstance(coeffs, list):
			raise TypeError("coefficients must be list, not %s" % coeffs.__class__.__name__)
		else:
			self.var = symbol
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
		return "DP(" + str(self.inner_vars) + ", " + self.sparse_rep() + ")"

	def __str__(self):
		return self.sparse_rep()

	"""
	Binary operations:
	add, sub, mul, flootdiv, mod, eq
	"""
	def __add__(f, g):
		return f._add(as_dp(g, f.var))

	def _add(f, g):
		if g.var == f.var:
			m = min(g.deg, f.deg) + 1
			coeffs_ = [f.coeffs[i] + g.coeffs[i] for i in range(m)]
			if f.deg > g.deg:
				coeffs_ += f.coeffs[m:]
			else:
				coeffs_ += g.coeffs[m:]
			return dp(f.var, coeffs_)
		else:
			add_ = bivariate_uniform(f, g)
			return add_[0] + add_[1]

	def __radd__(f, g):
		return f + g

	def __sub__(f, g):
		return f + g * (-1)

	def __rsub__(f, g):
	 	return - f + g

	def __mul__(f, g):
		return f._mul(as_dp(g, f.var)) 

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
			return dp(f.var, coeffs_)
		else:
			mul_ = bivariate_uniform(f, g)
			return mul_[0] * mul_[1]

	def __rmul__(f, g):
		return f * g

	def __eq__(f, g):
		if isinstance(g, DP):
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
			return dp(f.var, [c % g for c in f.coeffs]) 
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
			if isinstance(c, DP) and not c.is_monomial():
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

	def is_constant(self):
		if self.deg == 0:
			return True
		else:
			return False

	def is_monomial(self):
		if self.is_zero or self.deg == self.trdeg:
			return True
		else:
			return False

	def is_univariate(self):
		try:
			return self._is_univariate
		except AttributeError:
			for c in self.coeffs:
				if isinstance(c, DP):
					if c.is_constant():
						continue
					self._is_univariate = False
					return False
			self._is_univariate = True
			return True

	def degree(self, *var, total=False, as_dict=False, any_vars=False):
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
			if isinstance(c, DP):
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
			if isinstance(c, DP):
				coeffs_.append(c.as_list())
			else:
				coeffs_.append(c)
		return coeffs_

	def _inner_vars(self):
		vars_ = (self.var, )
		inner_ = tuple()
		for c in self.coeffs:
			if isinstance(c, DP):
				if len(c.inner_vars) > len(inner_):
					inner_ = c.inner_vars
			else:
				continue
		return vars_ + inner_

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

	def sort_vars(self, var):
		"""
		Sort variables to given order.
		
		Example:
		>>> a, b, c = symbols("a, b, c")
		>>> p = dp(a, [dp(b, [2, 0, dp(c, [1, 1])]), dp(b, [0, 3, dp(c, [1, 2, 1, 2])])])
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

def dp(symbol, coeffs):
	return DP(symbol, coeffs)

def dp_from_dict(var, rep):
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
		return dp(var[0], coeffs_)
	else:
		pass

def dp_from_list(var, rep):
	if not isinstance(var, tuple):
		raise TypeError("first argument must be tuple, not %s" % var.__class__.__name__)
	if len(var) == 0:
		raise ValueError("needs one variable at least")
	if len(var) == 1:
		return dp(var[0], rep)
	else:
		coeffs_ = []
		for r in rep:
			if isinstance(r, list):
				coeffs_.append(dp_from_list(var[1:], r))
			elif isinstance(r, int):
				coeffs_.append(r)
			else:
				raise ValueError("elements of coeffs list must be int or dp, not %s" % r.__class__.__name__)
		return dp(var[0], coeffs_)

def as_dp(f, *symbol):
	if isinstance(f, int):
		if len(symbol) == 0:
			raise ValueError("needs one symbol when f is int")
		else:
			return dp(symbol[0], [f])
	elif isinstance(f, Symbol):
		return f.as_dp()
	elif isinstance(f, DP):
		return f
	else:
		raise TypeError("argument must be int, Symbol, dp or DMP, not %s" % f.__class__.__name__)

def abs(rep):
	if isinstance(rep, int):
		if rep >= 0:
			return rep
		else:
			return - rep
	elif isinstance(rep, DP):
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
	>>> f = dp(x, [1, 1]) # f = x + 1
	>>> g = dp(y, [1, 1]) # g = y + 1
	>>> bivariate_uniform(f, g)
	(DMP((x, y), DMP(x, [dp(y, [1]), 1])), DMP((x,y), dp(x, [dp(y, [1, 1]), 1])))

	"""
	f_list = [dp(g.symbol, [f.coeffs[0]])] + f.coeffs[1:]
	f_ = dp(f.symbol, f_list)
	g_ = dp(f.symbol, [g])
	return (f_, g_)