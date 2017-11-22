FC=ifort

test : test_mt.o mt19937_mod.o
	$(FC) -o $@ $^

.SUFFIXES: .o .F90

.F90.o :
	$(FC) -c -o $@ $<

clean:
	@rm -rf *.o *.mod test

.PHONY: clean

test_mt.o : mt19937_mod.o
