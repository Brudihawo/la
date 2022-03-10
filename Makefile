CFLAGS := -Wall -Wextra -pedantic -std=c11 -g
LIBS := -lm

.DEFAULT_GOAL = test_solv

log:
	$(CC) log.c -c $(CFLAGS) $(LIBS)

la:
	$(CC) la.c -c $(CFLAGS) $(LIBS)

la_lib: la log
	ar rcs la.a log.o la.o

test_solv: la_lib
	$(CC) test/test_solv.c la.a -o test/test_solv $(CFLAGS) $(LIBS) -I..

clean:
	rm -vf *.o
	rm -vf *.a
	rm -vf test_solv
