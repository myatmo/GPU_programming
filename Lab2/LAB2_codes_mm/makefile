# use "make main" and "mpirun -np 4 -hostfile phi_cluster main.ex" to run main.c
# use "make mainStats" and "mpirun -np 8 -hostfile phi_cluster mainStats.ex" to run mainStats.c
.SUFFIXES: .c .o .f .F

CC			=  mpicc
CFLAGS			=  -g -Wall

FILES =  mainStats.o aux.o hqsort.o
FILES1 =  main.o aux.o hqsort.o

all: mainStats main

mainStats: $(FILES)
	${CC} ${CFLAGS} -o mainStats.ex -lm $(FILES)

main: $(FILES1)
	${CC} ${CFLAGS} -o main.ex -lm $(FILES1)

.c.o:
	${CC} ${CFLAGS} $< -c -o $@

clean:
	rm *.o *.ex
