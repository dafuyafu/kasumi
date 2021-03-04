def tuple_union(t, s):
	union_ = t
	for s_ in s:
		if s_ in union_:
			pass
		else:
			union_ += (s_, )
	return union_

def tuple_minus(t, s):
	minus_ = list(t)
	for s_ in s:
		if s_ in minus_:
			minus_.remove(s_)
		else:
			pass
	return tuple(minus_)

def tuple_intersection(t, s):
	inter_ = list(t)
	for s_ in s:
		if s_ in inter_:
			pass
		else:
			inter_.remove()
	return tuple(inter_)

def zero_tuple(i):
	return tuple(0 for i in range(i))

def validate_type(x, *types):
	if any(isinstance(x, t) for t in types):
		return None
	if len(types) == 1:
		str_ = types[0].__name__ + ", "
	else:
		str_, i = str(), 0
		while True:
			str_ += types[i].__name__
			if i == len(types) - 1:
				str_ += ", "
				break
			elif i == len(types) -2:
				str_ += " or "
				i += 1
			else:
				str_ += ", "
				i += 1
	type_ = x.__class__.__name__
	raise TypeError("argument must be " + str_ + "not " + type_)