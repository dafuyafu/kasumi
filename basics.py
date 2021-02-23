class DUP:
	"""

	Represents a dense univariate polynomial.
	Coefficients need to be integers.

	"""
	def __init__(self, symbol, coeffs=[]):
		self.symbol = symbol
		self.coeffs = []
		self.unset = True
		if not isinstance(coeffs, list):
			raise TypeError("needs to be list")
		else:
			if coeffs == []:
				pass
			else:
				self.set(coeffs)

	def set(self, coeffs):
		if self.unset:
			self.coeffs = coeffs
			self.deg = len(coeffs) - 1
			self.trdeg = self.trailing_deg()
			self.unset = False
		else:
			raise ValueError("Already set.")

	def __repr__(self):
		return self.sparse_rep()

	def __str__(self):
		return self.sparse_rep()

	def trailing_deg(self):
		for d in range(self.deg + 1):
			if self.coeffs[d] != 0:
				return d

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
				rep += self.symbol + "**" + str(deg)
			elif deg == 1:
				if c != "1":
					rep += c + "*"
				rep += self.symbol
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

	"""
	Binary operations:
	add, sub, mul, flootdiv, mod, eq
	"""

	def __add__(f, g):
		if isinstance(g, DUP):
			if g.symbol == f.symbol:
				m = min(g.deg, f.deg) + 1
				coeffs_ = [f.coeffs[i] + f.coeffs[i] for i in range(m)]
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

	def __sub__(f, g):
		return f + g * (-1)

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
				# return DMP
				pass
		if isinstance(g, int):
			coeffs_ = [f.coeffs[i] * g for i in range(f.deg + 1)]
			return dup(f.symbol, coeffs_)

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

class DMP:
	def __init__(self, symbols, rep):
		if not isinstance(rep, DUP):
			raise TypeError("needs to be a DUP object")
		self.symbols = symbols
		self.rep = rep

def dup(s, co):
	return DUP(s, co)

def expr(rep):
	return Expr(rep)

def symbol(c):
	return expr(dup(c, [0,1]))

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
