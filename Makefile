SRCDIR := ./src
TESTDIR := ./test
EXDIR := ./examples

objs := $(patsubst %.c, %.o, $(wildcard $(SRCDIR)/*.c))
benchs := $(patsubst %.c, %, $(wildcard $(TESTDIR)/bench_*.c))
tests := $(patsubst %.c,%,$(wildcard $(TESTDIR)/test_*.c))
examples := $(patsubst %.c,%,$(wildcard $(EXDIR)/*.c))

CFLAGS := -Wall -Wextra -pedantic -Wuninitialized -std=c11
SANIFLAGS := -fsanitize=address -fsanitize=leak -fsanitize=undefined -fsanitize=bounds
LIBS := -lm
CC := clang

.PHONY = example test bench
.DEFAULT_GOAL = test

ifndef ECHO
	HIT_TOTAL != ${MAKE} ${MAKECMDGOALS} --dry-run ECHO="HIT_MARK" | grep -c "HIT_MARK"
	HIT_COUNT = $(eval HIT_N != expr ${HIT_N} + 1) ${HIT_N}

	ECHO = echo -e "\033[1;32m[`expr ${HIT_COUNT} '*' 100 / ${HIT_TOTAL}`%]\033[0m"
endif

ifeq ("$(DEBUG)", "true")
  CFLAGS += -g -Og
else
  CFLAGS += -O2
endif

ifeq ("$(SANITIZE)", "true")
  CFLAGS += $(SANIFLAGS)
endif


$(objs): %.o: %.c
	@$(ECHO) $@
	$(CC) $(CFLAGS) $< -c -o $@

$(SRCDIR)/la.a: $(objs)
	@$(ECHO) $@
	ar rcs $(SRCDIR)/la.a $(SRCDIR)/*.o

$(tests): %: %.c $(SRCDIR)/la.a
	@$(ECHO) $@
	$(CC) $< $(SRCDIR)/la.a -o $@ $(CFLAGS) $(LIBS) -I$(SRCDIR)

$(benchs): %: %.c $(SRCDIR)/la.a
	@$(ECHO) $@
	$(CC) $< $(SRCDIR)/la.a -o $@ $(CFLAGS) $(LIBS) -I$(SRCDIR)

$(examples): %: %.c $(SRCDIR)/la.a
	@$(ECHO) $@
	$(CC) $< $(SRCDIR)/la.a -o $@ $(CFLAGS) $(LIBS) -I$(SRCDIR)

test: $(tests)

runtest: test
	./run_tests.sh $(TESTDIR)

example: $(examples)

bench: $(benchs)

all: bench example test

coverage: clean test
	./coverage.sh

clean:
	rm -vf $(SRCDIR)/*.o
	rm -vf $(SRCDIR)/*.a
	rm -vf $(tests)
	rm -vf $(benchs)
	rm -vf $(examples)
	rm -vfr coverage
