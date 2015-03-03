def update(rn):
	return 1, rn + 1

readNum = 0
while(True):
	print readNum
	x, y = update(readNum)
	readNum = y
