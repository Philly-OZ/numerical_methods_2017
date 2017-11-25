CC = gcc
CFLAGS = -W -Wall
SRC = $(wildcard ./src/*.c)
OBJ = $(SRC:.c=.o)

main : $(OBJ)
	$(CC) -o $@ $^ 

$.o : $.c
	$(CC) -o $@ -c $< $(CFLAGS)

clean :
	rm -rf ./src/*.o

mrproper : clean
	rm main 