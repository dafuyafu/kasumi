from math import sqrt, floor

def is_prime(n):
	if n == 2:
		return True
	if n % 2 == 0:
		return False
	sqrt_ = floor(sqrt(n)) + 1
	for i in range(3, sqrt_, 2):
		if n % i == 0:
			return False
	return True