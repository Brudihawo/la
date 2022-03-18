CFLAGS := -Wall -Wextra -pedantic -Wuninitialized -std=c11 -g
TEST_FLAGS := -fsanitize=address -fsanitize=leak -fsanitize=undefined -fsanitize=bounds
LIBS := -lm
CC := clang

SRCDIR := ./src
TESTDIR := ./test

.DEFAULT_GOAL = test

log:
	$(CC) $(SRCDIR)/log.c -c -o $(SRCDIR)/log.o $(CFLAGS)

la:
	$(CC) $(SRCDIR)/la.c -c -o $(SRCDIR)/la.o $(CFLAGS)

matf:
	$(CC) $(SRCDIR)/matf.c -c -o $(SRCDIR)/matf.o $(CFLAGS)

sparse:
	$(CC) $(SRCDIR)/sparse.c -c -o $(SRCDIR)/sparse.o $(CFLAGS) $(TEST_FLAGS)

la_lib: la matf sparse log
	ar rcs $(SRCDIR)/la.a $(SRCDIR)/*.o

test_solv: la_lib
	$(CC) $(TESTDIR)/test_solv.c $(SRCDIR)/la.a -o $(TESTDIR)/test_solv \
		$(CFLAGS) $(TEST_FLAGS) $(LIBS) -I./src

test_sparse:
	$(CC) $(TESTDIR)/test_sparse.c $(SRCDIR)/la.a -o $(TESTDIR)/test_sparse \
		$(CFLAGS) $(TEST_FLAGS) $(LIBS) -I./src

test: test_solv test_sparse

coverage: clean test
	./coverage.sh

clean:
	rm -vf $(SRCDIR)*.o
	rm -vf $(SRCDIR)*.a
	rm -vf $(TESTDIR)/test_solv
	rm -vf $(TESTDIR)/test_sparse
	rm -vfr coverage
