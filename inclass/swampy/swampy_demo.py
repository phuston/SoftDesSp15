""" TurtleWorld Demo """

from swampy.TurtleWorld import *

def make_square(init_x, init_y, side):
	world = TurtleWorld()
	pat = Turtle()
	pat.set_delay(.001)

	for i in range(400):
		pat.x = init_x-(i*5/2)
		pat.y = init_y-(i*5/2)
		for j in range(4):
			pat.fd(side+i*5)
			pat.lt(90)
	wait_for_user()

def make_polygon(init_x, init_y, sides, length):
	world = TurtleWorld()
	pat = Turtle()
	pat.set_delay(.001)
	angle_int = 360/sides

	pat.x = init_x
	pat.y = init_y
	for i in range(sides):
		pat.fd(length)
		pat.lt(angle_int)

	wait_for_user()

def make_square_2(init_x, init_y, length):
	make_polygon(init_x, init_y, 4, length)


make_square(0,0,20)
make_polygon(0,0,20,20)
make_square_2(0,0,24)