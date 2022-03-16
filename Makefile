CFLAGS := -Wall -Wextra -pedantic -std=c11 -g
LIBS := -lm

SRCDIR := ./src
TESTDIR := ./test

.DEFAULT_GOAL = test

log:
	$(CC) $(SRCDIR)/log.c -c -o $(SRCDIR)/log.o $(CFLAGS) $(LIBS)

la:
	$(CC) $(SRCDIR)/la.c -c -o $(SRCDIR)/la.o $(CFLAGS) $(LIBS)

la_lib: la log
	ar rcs $(SRCDIR)/la.a $(SRCDIR)/log.o $(SRCDIR)/la.o

test_solv: la_lib
	$(CC) $(TESTDIR)/test_solv.c $(SRCDIR)/la.a -o $(TESTDIR)/test_solv \
		$(CFLAGS) $(LIBS) -I./src

test_sparse:
	$(CC) $(TESTDIR)/test_sparse.c $(SRCDIR)/la.a -o $(TESTDIR)/test_sparse \
		$(CFLAGS) $(LIBS) -I./src

test: test_solv test_sparse

coverage: clean test
	./coverage.sh

clean:
	rm -vf $(SRCDIR)*.o
	rm -vf $(SRCDIR)*.a
	rm -vf $(TESTDIR)/test_solv
	rm -vf $(TESTDIR)/test_sparse
	rm -vfr coverage
