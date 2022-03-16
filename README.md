# Linear Algebra Library in C

Provides Matrix, Vector Multiplication etc. This is a work in progress.

# TODO

- [ ] SMatF Initialisation functions
- [ ] SMatF Deinitialisation functions
- [ ] Verify SMatF product correctness with tests (maybe rely on MatF Multiplication)
- [ ] SMatF Iterative Methods for linear systems of equations
- [ ] Benchmark SMatF vs MatF product
- [ ] Sorting algorithms so i dont require some inputs to be sorted

## In-File-TODOs:

- [ ] `./src/la.c:64`: TODO: can i somehow check that rows * cols = sizeof vals?
- [ ] `./src/la.c:129`: TODO: Handle size checks nicely
- [ ] `./src/la.c:166`: TODO: Proper error reporting for size checks
- [ ] `./src/la.c:204`: TODO: parallelisation via openmp
- [ ] `./src/la.c:328`: TODO: implement stability measures in gauss eliminiation
- [ ] `./src/la.c:381`: TODO: implement pivoting and stability measures in gauss eliminiation
- [ ] `./src/la.c:410`: TODO: check for invertibility of A
- [ ] `./src/la.c:430`: TODO: implement pivoting and stability measures in gauss eliminiation
- [ ] `./src/la.c:431`: TODO: do i want to zero the right lower half of A?
- [ ] `./src/la.c:641`: TODO: check correctness
- [ ] `./src/la.c:673`: TODO: do i want a non-panicking setter?
- [ ] `./src/la.c:819`: TODO: Do i need column-wise iteration in SMatF?
- [ ] `./src/la.h:302`: TODO: Update Sparse Matrix structure documentation
- [ ] `./src/la.h:336`: TODO: Implement SM_empty_from_pos
- [ ] `./src/la.h:340`: TODO: Implement SM_empty_from_pos_vals
- [ ] `./test/test_solv.c:121`: TODO: Write proper tests

