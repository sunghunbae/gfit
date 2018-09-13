CC = g++
CCOPT = -O
GSLLIB = -static -L/usr/lib -lgsl -lgslcblas
HEADERS = gfit.h
OBJECTS = func.o gridsearch.o fitcurve.o levmar.o brent.o file.o \
	data.o print.o main.o

gfit: $(OBJECTS)
	$(CC) -o gfit $(OBJECTS) $(GSLLIB) 

func.o: func.c $(HEADERS)
	$(CC) -c $(CCOPT) func.c
fitcurve.o: fitcurve.c $(HEADERS)
	$(CC) -c $(CCOPT) fitcurve.c
levmar.o: levmar.c $(HEADERS)
	$(CC) -c $(CCOPT) levmar.c
brent.o: brent.c $(HEADERS)
	$(CC) -c $(CCOPT) brent.c
file.o: file.c $(HEADERS)
	$(CC) -c $(CCOPT) file.c
gridsearch.o: gridsearch.c $(HEADERS)
	$(CC) -c $(CCOPT) gridsearch.c
data.o: data.c $(HEADERS)
	$(CC) -c $(CCOPT) data.c
print.o: print.c $(HEADERS)
	$(CC) -c $(CCOPT) print.c
main.o: main.c $(HEADERS)
	$(CC) -c $(CCOPT) main.c
