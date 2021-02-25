class Symbol:
	def __init__(self, char):
		if not isinstance(char, str):
			raise TypeError("can generate Symbol only with str, not %s", char.__class__.__name__)
		self.char = char

	def __repr__(self):
		return self.char

	def __str__(self):
		return self.char

	def __add__(a, b):
		if isinstance(b, Symbol):
			return a.as_dup() + b.as_dup()
		else:
			return a.as_dup() + b

	def __sub__(a, b):
		if isinstance(b, Symbol):
			return a.as_dup() - b.as_dup()
		else:
			return a.as_dup() - b

	def __mul__(a, b):
		if isinstance(b, Symbol):
			return a.as_dup() * b.as_dup()
		else:
			return a.as_dup() * b

	def __radd__(a, b):
		return a + b

	def __rsub__(a, b):
		return - a + b

	def __rmul__(a, b):
		return a * b

	def __neg__(a):
		return a * -1
 
	def as_dup(self):
		return dup(self, [0, 1])

def symbol(c):
	if not isinstance(c, str):
		raise TypeError("symbol of variable must be str, not %s", c.__class__.__name__)
	return Symbol(c)

def symbols(c):
	if not isinstance(c, str):
		raise TypeError("symbol of variable must be str, not %s", c.__class__.__name__)
	symbols_ = c.replace(" ", "").split(",")
	return tuple([symbol(s) for s in symbols_])

class DUP:
	"""

	Represents a dense univariate polynomial.
	Coefficients need to be integers.

	"""


	def __init__(self, symbol, coeffs=[]):
		if not isinstance(symbol, Symbol):
			raise TypeError("symbol must be Symbol, not %s", symbol.__class__.__name__)
		elif not isinstance(coeffs, list):
			raise TypeError("coefficients must be list, not %s", coeffs.__class__.__name__)
		else:
			self.symbol = symbol
			if coeffs == []:
				self.set([0])
			else:
				self.set(coeffs)

	def set(self, coeffs):
		self.coeffs = coeffs
		self.deg = len(coeffs) - 1
		self.trdeg = self.trailing_deg()
		self.unset = False

	def __repr__(self):
		return self.sparse_rep()

	def __str__(self):
		return self.sparse_rep()

	"""
	Binary operations:
	add, sub, mul, flootdiv, mod, eq
	"""

	def __add__(f, g):
		if isinstance(g, DUP):
			if g.symbol == f.symbol:
				m = min(g.deg, f.deg) + 1
				coeffs_ = [f.coeffs[i] + g.coeffs[i] for i in range(m)]
				if f.deg > g.deg:
					coeffs_ += f.coeffs[m:]
				else:
					coeffs_ += g.coeffs[m:]
				return dup(f.symbol, coeffs_)
			else:
				# return DMP
				pass
		if isinstance(g, int):
			coeffs_ = [f.coeffs[0] + g] + f.coeffs[1:]
			return dup(f.symbol, coeffs_)

	def __radd__(f, g):
		return f + g

	def __sub__(f, g):
		return f + g * (-1)

	def __rsub__(f, g):
	 	return - (f - g)

	def __mul__(f, g):
		if isinstance(g, DUP):
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
				var = (f.symbol, g.symbol)

				pass
		if isinstance(g, int):
			coeffs_ = [f.coeffs[i] * g for i in range(f.deg + 1)]
			return dup(f.symbol, coeffs_)

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
			raise TypeError("argument must be dict, not %s", subsdict.__class__.__name__)
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
		raise TypeError("first argument must be tuple, not %s", var.__class__.__name__)
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
		raise TypeError("first argument must be tuple, not %s", var.__class__.__name__)
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
				raise ValueError("elements of coeffs list must be int or DUP, not %s", r.__class__.__name__)
		return dup(var[0], coeffs_)

class DMP:
	def __init__(self, symbols, rep):
		if not isinstance(rep, DUP):
			raise TypeError("rep must be DUP, not %s", rep.__class__.__name__)
		self.symbols = symbols
		self.rep = rep

	def __repr__(self):
		return self.rep.sparse_rep()

	def __str__(self):
		return self.rep.sparse_rep()

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