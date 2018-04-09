#####################################################################
#
# Typing the command below => results in the following actions
#    make             =>      Executable a.out
#    make clean       =>      Remove object files, executables
#
#####################################################################

FF = ifort
#CFLAG = -O3
#LIB = -mkl

FILES = h2opot.o main.o

all: pes.out
pes.out: $(FILES)
	$(FF) $(LIB) $(CFLAG) $(FILES) -o $@

clean:
	rm -f *.o pes.out

.SUFFIXES:

h2opot.o: h2opot.f90
	$(FF) -c h2opot.f90

main.o: main.f90
	$(FF) -c main.f90
