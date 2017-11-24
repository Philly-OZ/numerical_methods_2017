clean:
	rm -f *.o
	rm -f main

main: main.c linear_project_generator.o
	cc 