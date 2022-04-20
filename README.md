# Linear Algebra Library in C

Provides Matrix, Vector Multiplication etc. This is a work in progress.

# Building

## Quickstart

How to get started and run the tests:
```shell
$ git clone https://github.com/brudihawo/la
$ cd la
$ make
```

## Debug / Sanitizing
To turn on debug mode / activate sanitizing flags:

```shell
$ make DEBUG=true SANITIZE=true
```

# TODO

- [ ] SMatF Initialisation functions
- [x] SMatF Deinitialisation functions
- [x] Verify SMatF product correctness with tests (maybe rely on MatF Multiplication)
- [ ] SMatF Iterative Methods for linear systems of equations
- [x] Benchmark SMatF vs MatF product
- [ ] valgrind for profiling
- [x] SMatF check for nvals < rows * cols
- [x] Build System improvement with debug and release targets
- [ ] Optimize forward and backward substituion for SMatF

## In-File-TODOs:

- [ ] `./src/matf.c:37`: TODO: can i somehow check that rows * cols = sizeof vals?
- [ ] `./src/matf.c:102`: TODO: Handle size checks nicely
- [ ] `./src/matf.c:139`: TODO: Proper error reporting for size checks
- [ ] `./src/matf.c:177`: TODO: parallelisation via openmp
- [ ] `./src/matf.c:301`: TODO: implement stability measures in gauss eliminiation
- [ ] `./src/matf.c:354`: TODO: implement pivoting and stability measures in gauss eliminiation
- [ ] `./src/matf.c:383`: TODO: check for invertibility of A
- [ ] `./src/matf.c:403`: TODO: implement pivoting and stability measures in gauss eliminiation
- [ ] `./src/matf.c:404`: TODO: do i want to zero the right lower half of A?

- [ ] `./src/sparse.c:300`: TODO: check correctness
- [ ] `./src/sparse.c:333`: TODO: I think this can be optimized. I dont think we have to do all these
- [ ] `./src/sparse.c:369`: TODO: do i want a non-panicking setter?
- [ ] `./src/sparse.c:465`: TODO: think about checking querying nonzero in A or B first, but i
- [ ] `./src/sparse.c:482`: TODO: think about checking querying nonzero in A or B first, but i
- [ ] `./src/sparse.c:564`: TODO: Do i need column-wise iteration in SMatF?

- [ ] `./src/sparse.h:5`: TODO: Update Sparse Matrix structure documentation

- [ ] `./test/test_solv.c:122`: TODO: Write proper tests
