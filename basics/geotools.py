from pys.pytools import validate_type
from basics.basictools import Symbol, dp, DP
from basics.domains import ring, Ring

class Point:
	def __init__(self, value, *var, **options):
		validate_type(value, dict, tuple, list)
		if isinstance(value, dict):
			self.value_dict = value
		else:
			self.value_dict = dict()
			if not var:
				raise ValueError("require variable with value tuple")
			i = 0
			for v in var:
				self.value_dict[v] = valur[0]
				i += 1

		if "dom" in options:
			self.dom = options["dom"]
		else:
			self.dom = ring()

		if "proj" in options and options["proj"]:
			self.affine = False
		else:
			self.affine = True

	def __hash__(self):
		pass

	def __eq__(p, q):
		if self.affine:
			if p.value_dict == q.value_dict:
				return True
			else:
				return False
		else:
			if p.value_dict == q.value_dict:
				return True
			else:
				q_ = r_  = q.value_dict
				for c in self.dom.it_elements():
					if c == 0:
						continue
					for k in q_:
						q_[k] *= c
					if p.value_dict == q_:
						return True
					else:
						q_ = r_
				return False

	"""
	* Representation methods
	"""

	def __repr__(self):
		return "Point(" + str(self) + ", dom = " + str(self.dom) + ")"

	def __str__(self):
		i = 0
		if self.affine:
			str_ = "("
		else:
			str_ = "["
		for p in self:
			str_ += str(self[p])
			if i < len(self) - 1:
				if self.affine:
					str_ += ", "
				else:
					str_ += ": "
			i += 1
		if self.affine:
			return str_ + ")"
		else:
			return str_ + "]"

	def as_dict(self):
		return self.value_dict

	def as_tuple(self):
		return tuple(v for v in self.value_dict.values())

	def __iter__(self):
		return iter(self.value_dict)

	def __len__(self):
		return len(self.value_dict)

	def __getitem__(self, key):
		return self.value_dict[key]

	def __setitem__(self, key, value):
		self.value_dict[key] = value

	def __mul__(p, n):
		validate_type(n, int, DP)
		for v in self:
			self[v] *= n

	def get_domain(self):
		return self.dom

def point(value, *var, **options):
	return Point(value, *var, **options)