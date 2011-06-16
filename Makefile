#--------------------------------------------------------------------
F90 = ftn
LFLAGS = 
FFLAGS = -Mpre -Minline=levels:11 -Mipa -Munroll
FFLAGS_IO = 
FFLAGS_MODEL = -Mpre -Mdse -Minline=levels:9,size:80 -Mvect=sse,noaltcode,uniform -Mnozerotrip -Munroll=n:4 
DEBUG = -g -Mbounds
MPI = -D MPI
#--------------------------------------------------------------------

OBJ = Mersenne.o Params.o Parallel.o Quaternions.o Functions.o Model.o IO.o Cashew.o

default:
	@echo " par    - parallel MPI version with optimal compiler flags"
	@echo " debug  - parallel MPI version, no optimization and with"
	@echo "          debugger handles"
	@echo " serial - serial version, mild compiler optimization"
	@echo " clean  - deletes existing object files"

#----------------------------------------------------------

par	: $(OBJ)
	$(F90) $(LFLAGS) $(OBJ) -o cashew

debug	: 
	$(F90) $(DEBUG) Mersenne.F90 Params.F90 Parallel.F90 Quaternions.F90 Functions.F90 Model.F90 IO.F90 Cashew.F90   -o cashew-debug

serial: 
	$(F90) Mersenne.F90 Params.F90 Parallel.F90 Quaternions.F90 Functions.F90 Model.F90 IO.F90 Cashew.F90 -o cashew-serial

.PHONY: clean

clean:
	\rm *.o *.mod 

#----------------------------------------------------------



Model.o: Model.F90
	$(F90) $(MPI) $(FFLAGS_MODEL)  -c $<
IO.o: IO.F90
	$(F90) $(MPI) $(FFLAGS_IO)  -c $<
%.o: %.F90
	$(F90) $(MPI) $(FFLAGS) -c $<



