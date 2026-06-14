# Flags de Compilação
CFLAGS = -O3 -march=native -mavx -fopt-info-vec -g
FLAGS_VALGRIND = --leak-check=full

nomePrograma = broyden

all: $(nomePrograma)

$(nomePrograma): main.o broyden.o utils.o
	gcc main.o broyden.o utils.o -o $(nomePrograma) -lm

main.o: main.c
	gcc -c main.c $(CFLAGS)

broyden.o: broyden.h broyden.c
	gcc -c broyden.c $(CFLAGS)

utils.o: utils.h utils.c
	gcc -c utils.c $(CFLAGS)

valgrind: $(nomePrograma)
	valgrind $(FLAGS_VALGRIND) ./$(nomePrograma)

clean:
	rm -f *.o *.gch $(nomePrograma)