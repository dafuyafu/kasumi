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
	Coefficients need to be integers.

	"""

	def __init__(self, symbol, coeffs=[]):
		if not isinstance(symbol, Symbol):
			raise TypeError("symbol must be Symbol, not %s" % symbol.__class__.__name__)
		elif not isinstance(coeffs, list):
			raise TypeError("coefficients must be list, not %s" % coeffs.__class__.__name__)
		else:
			self.symbol = symbol
			if coeffs == []:
				self._zero()
			else:
				self._set(coeffs)

	def _set(self, coeffs):
		self.is_zero = False
		self.coeffs = coeffs
		self.deg = len(coeffs) - 1
		self.trdeg = self.trailing_deg()

	def _zero(self):
		self.is_zero = True
		self.coeffs = [0]
		self.deg = None
		self.trdeg = None

	def __repr__(self):
		return self.sparse_rep()

	def __str__(self):
		return self.sparse_rep()

	"""
	Binary operations:
	add, sub, mul, flootdiv, mod, eq
	"""
	def __add__(f, g):
		return f._add(as_dup(g, f.symbol))

	def _add(f, g):
		if g.symbol == f.symbol:
			m = min(g.deg, f.deg) + 1
			coeffs_ = [f.coeffs[i] + g.coeffs[i] for i in range(m)]
			if f.deg > g.deg:
				coeffs_ += f.coeffs[m:]
			else:
				coeffs_ += g.coeffs[m:]
			return dup(f.symbol, coeffs_)
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
		return f._mul(as_dup(g, f.symbol)) 

	def _mul(f, g):
		if g.symbol == f.symbol:
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
			return dup(f.symbol, coeffs_)
		else:
			mul_ = bivariate_uniform(f, g)
			return mul_[0] * mul_[1]

	def __rmul__(f, g):
		return f * g

	def __eq__(f, g):
		if isinstance(g, DUP):
			if f.symbol == g.symbol and f.coeffs == g.coeffs:
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
			return dup(f.symbol, [c % g for c in f.coeffs]) 
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
				rep += str(self.symbol) + "**" + str(deg)
			elif deg == 1:
				if c != "1":
					rep += c + "*"
				rep += str(self.symbol)
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

	def degree(self):
		return self.deg

	def coeffs_list(self):
		coeffs_ = []
		for c in self.coeffs:
			if isinstance(c, DUP):
				coeffs_.append(c.coeffs_list())
			else:
				coeffs_.append(c)
		return coeffs_

	def inner_vars(self):
		vars_ = (self.symbol, )
		for c in self.coeffs:
			if isinstance(c, DUP):
				vars_ += c.inner_vars()
			else:
				continue
		return vars_

	def subs(self, subsdict):
		if not isinstance(subsdict, dict):
			raise TypeError("argument must be dict, not %s" % subsdict.__class__.__name__)
		if not self.symbol in subsdict:
			return self
		else:
			subs_ = 0
			value_ = subsdict[self.symbol]
			for i in range(self.deg + 1):
				subs_ += self.coeffs[i] * value_ ** i
			return subs_

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

class DMP:

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
		elif True:
			pass

	def __radd__(f, g):
		return f + g

	def variable_sort(self, t):
		if not isinstance(t, tuple):
			raise TypeError("needs one tuple not %s", t.__class__.__name__)
		if not set(self.symbols) == set(t):
			raise ValueError("can sort only if both polynomial have the same variables")
		if self.symbols == t:
			raise ValueError("already sorted")

	def degree(self, *var):
		pass

	def subs(self, subsdict):
		pass

	def same_variables(self, other):
		if set(self.symbols) == set(other.symbols):
			return True
		else:
			return False

def dmp(symbols, rep):
	return DMP(symbols, rep)

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
	symbols_ = (f.symbol, g.symbol)
	f_list = [dup(g.symbol, [f.coeffs[0]])] + f.coeffs[1:]
	f_ = dmp(symbols_, dup(f.symbol, f_list))
	g_ = dmp(symbols_, dup(f.symbol, [g]))
	return (f_, g_)