from basics import DP, dp

def validator(func):
	def vldtr_(f):
		if isinstance(f, DP):
			pass
		else:
			raise TypeError("argument must be DP object, not '%s'" % f.__class__.__name__)
	return vldtr_

@validator
def LT(f):
	return 
	