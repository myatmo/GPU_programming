.SUFFIXES: .c .o .f .F

CC			=  mpicc
CFLAGS			=  -g -Wall

FILES =  main.o aux.o hqsort.o mainStats.o

main.ex: $(FILES)
	${CC} ${CFLAGS} -o main.ex -lm $(FILES)

.c.o:
	${CC} ${CFLAGS} $< -c -o $@

clean:
	rm *.o *.ex
