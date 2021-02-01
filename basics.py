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
		for deg in range(self.deg + 1)[::-1]:
			if self.coeffs[deg] == 0:
				continue		
			if deg > 1:
				if self.coeffs[deg] != 1:
					rep += str(self.coeffs[deg]) + "*"
				rep += self.symbol + "**" + str(deg)
			elif deg == 1:
				if self.coeffs[deg] != 1:
					rep += str(self.coeffs[deg]) + "*"
				rep += self.symbol
			else:
				rep += str(self.coeffs[0])
			if deg != self.trdeg:
				rep += " + "
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
			coeffs_ = [f.coeffs[0] + g] + f.coeffs[:1]
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
			coeffs_ = [f.coeffs[i] * g for i in range(f.deg)]
			return dup(f.symbol, coeffs_)

class DMP:
	def __init__():
		pass

def dup(s, co):
	return DUP(s, co)

def expr(rep):
	return Expr(rep)

def symbol(c):
	return expr(dup(c, [0,1]))