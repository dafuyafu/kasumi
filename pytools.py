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