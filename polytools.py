from basics import dp, DP

class Poly:
	"""
	represents an element of polynomial ring

	
	"""

	def __new__(cls, var, rep, mod=0, quo=0, rel=0):
		if isinstance(rep, int):
			self = super().__new__(Integer)
		elif isinstance(rep, DP):
			self = super().__new__(cls)
		return self

	def __init__(self, var, rep, mod=0, quo=0, rel=0):
		self.rep = rep

	def __repr__(self):
		return "%s: %s" % (self.__class__.__name__, str(self.rep))

class Constant(Poly):
	pass

class Integer(Constant):
	pass

def poly(f):
	return Poly(f)