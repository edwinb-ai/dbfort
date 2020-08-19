CC = gfortran
FLAGS = -O3 -g -llapack -Wextra
objects = utils.o tensor.o positions.o
target = dbrown

dbrown: $(objects)
	$(CC) $(target).f90 -o $(target) $(objects) $(FLAGS)

positions.o: positions.f90
	$(CC) -c positions.f90

utils.o: utils.f90
	$(CC) -c utils.f90

tensor.o: tensor.f90
	$(CC) -g -c tensor.f90

clean:
	rm $(objects) $(target)