CC = gfortran
FLAGS = -Wall -Wextra -Wimplicit-interface -fPIC -fmax-errors=1 -g -fcheck=all -fbacktrace
LIBS = -llapack
objects = positions.o tensor.o utils.o outerprod.o random.o parameters.o
target = dbrown
INCLUDE = include/
SOURCES = src/

dbrown: $(objects)
	$(CC) $(FLAGS) main/$(target).f90 -o $(target) $(objects) $(LIBS) -I$(INCLUDE)

parameters.o: $(SOURCES)parameters.f90
	$(CC) $(FLAGS) -c $(SOURCES)parameters.f90 -J$(INCLUDE) -I$(INCLUDE)

outerprod.o: $(SOURCES)outerprod.f90
	$(CC) $(FLAGS) -c $(SOURCES)outerprod.f90 -J$(INCLUDE) -I$(INCLUDE)

positions.o: $(SOURCES)positions.f90
	$(CC) $(FLAGS) -c $(SOURCES)positions.f90 -J$(INCLUDE) -I$(INCLUDE)

utils.o: $(SOURCES)utils.f90
	$(CC) $(FLAGS) -c $(SOURCES)utils.f90 -J$(INCLUDE) -I$(INCLUDE)

tensor.o: $(SOURCES)tensor.f90
	$(CC) $(FLAGS) -c $(SOURCES)tensor.f90 -J$(INCLUDE) -I$(INCLUDE)

random.o: $(SOURCES)random.f90
	$(CC) $(FLAGS) -c $(SOURCES)random.f90 -J$(INCLUDE) -I$(INCLUDE)


clean:
	rm $(objects) $(target)