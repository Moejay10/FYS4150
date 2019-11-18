CPPflags = c++ -std=c++11

LIB = -larmadillo

# Test of functions

test :	  		test.o jacobi.o
				${CPPflags} test.o jacobi.o ${LIB} -o test.exe

test.o :		test.cpp 
				${CPPflags} -c test.cpp

jacobi.o :		jacobi.cpp 
				${CPPflags} -c jacobi.cpp

# Buckling beam

buckling_beam:		buckling_beam.o jacobi.o
				${CPPflags} buckling_beam.o jacobi.o ${LIB} -o buckling_beam.exe

buckling_beam.o :	buckling_beam.cpp 
		        	${CPPflags} -c buckling_beam.cpp


# Quantum dot with one electron

quantum_one:		quantum_one.o jacobi.o
				${CPPflags} quantum_one.o jacobi.o ${LIB} -o quantum_one.exe

quantum_one.o :		quantum_one.cpp 
		        	${CPPflags} -c quantum_one.cpp

# Quantum dot with two electrons

quantum_two:		quantum_two.o jacobi.o
				${CPPflags} quantum_two.o jacobi.o ${LIB} -o quantum_two.exe

quantum_two.o :		quantum_two.cpp 
		        	${CPPflags} -c quantum_two.cpp
