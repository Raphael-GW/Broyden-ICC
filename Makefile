# Adicionado -g para incluir símbolos de depuração, essencial para o Valgrind
parametrosCompilacao=-O3 -march=native -mavx -fopt-info-vec -lm;
nomePrograma=broyden

# Flags recomendadas para uma análise completa de memória
FLAGS_VALGRIND=--leak-check=full --show-leak-kinds=all --track-origins=yes

all: $(nomePrograma)

$(nomePrograma): main.o broyden.o utils.o
	gcc -o $(nomePrograma) main.o broyden.o utils.o $(parametrosCompilacao)

main.o: main.c
	gcc -c main.c $(parametrosCompilacao)

broyden.o: broyden.h broyden.c
	gcc -c broyden.c $(parametrosCompilacao)

utils.o: utils.h utils.c
	gcc -c utils.c $(parametrosCompilacao)

# Nova regra: compila (se necessário) e executa com valgrind
valgrind: $(nomePrograma)
	valgrind $(FLAGS_VALGRIND) ./$(nomePrograma)

clean:
	rm -f *.o *.gch $(nomePrograma)