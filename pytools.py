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