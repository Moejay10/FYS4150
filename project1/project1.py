import sys
import numpy as np
import matplotlib.pyplot as plt

# Read the highest power of 10 from the command line
if len(sys.argv) <= 1:
    print "Error:", sys.argv[0], "reads max power of 10"
    exit(1)
else:
    N = int(sys.argv[1])

# Define arrays for log10(h) and max relative error
h = np.zeros(N)
max_error = np.zeros(N)

for I in range(1, N+1):
    # Define the number of grid points
    n = 10**I + 2

    # Define an array which stores the data from data-I-.dat
    data = np.zeros((n,4))

    # Open input data file and read in the results
    datafile = "data{}.dat".format(I)
    with open(datafile, 'r') as infile:
        infile.readline()
        for i in range(n):
            data[i] = infile.readline().split()
        infile.readline()
        h[I-1] = infile.readline()
    
    # Find the max value of the relative error
    max_error[I-1] = max(data[1:-1,3])
    
    # Plot the exact and approximate solution
    plt.figure(I)
    plt.plot(data[:,0],data[:,1],label="Approx.")
    plt.plot(data[:,0],data[:,2], label="Exact")
    plt.xlabel("x")
    plt.ylabel("u(x)")
    plt.legend()
    plt.savefig("plot{}.png".format(I))

# Plot epsilon vs log10(h)
plt.figure(N+1)
plt.plot(h,max_error)
plt.xlabel("log10(h)")
plt.ylabel("epsilon")
plt.savefig("plot_error.png")
