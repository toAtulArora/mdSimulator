#Makefile adapted from (has been modified)
#http://stackoverflow.com/questions/5871780/creating-a-fortran-makefile

#Compiler
FC = gfortran -fbackslash

#source files, use the most dependent in the end
SRCS = $(addsuffix .o, gnuplot_fortran bounce)

PROGRAM = bounce

all: $(PROGRAM)

$(PROGRAM): $(SRCS)
	$(FC) -o $@ $^

%.o : %.f95
	$(FC) -o $@ -c $<

clean:
	rm -f *.o *.mod
