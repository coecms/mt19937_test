FC=ifort
FLAGS=-i8

test : test_mt.o mt19937_mod.o
	$(FC) $(FLAGS) -o $@ $^

.SUFFIXES: .o .F90

.F90.o :
	$(FC) $(FLAGS) -c -o $@ $<

clean:
	@rm -rf *.o *.mod test

.PHONY: clean

test_mt.o : mt19937_mod.o
