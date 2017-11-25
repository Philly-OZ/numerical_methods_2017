CC = gcc
CFLAGS = -W -Wall
SRC = $(wildcard *.c)
OBJ = $(SRC:.c=.o)

main : $(OBJ)
	$(CC) -o $@ $^ 

$.o : $.c
	$(CC) -o $@ -c $< $(CFLAGS)

clean :
	rm -rf *.o

mrproper : clean
	rm main 