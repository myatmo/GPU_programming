.SUFFIXES: .c .o .f .F

CC			=  mpicc  
CFLAGS			=  -g -Wall

FILES =  main.o pipe_ge.o back_solve.o

main.ex: $(FILES) 
	${CC} ${CFLAGS} -o main.ex -lm $(FILES) 

.c.o:
	${CC} ${CFLAGS} $< -c -o $@

clean:
	rm *.o *.ex
