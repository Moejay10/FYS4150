import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from matplotlib.pyplot import cm
from tqdm import tqdm
from matplotlib import rc


#Latex font for plots
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
rc('font', family='serif')
plt.rcParams.update({'font.size': 10}) # Setting all font sizes


print("Which Task you want to run 1-dim , 2-dim or Lithosphere?")
print("Write 1 , 2 or 3")

Task = input("Write here: ")

w = 5.78851          # Latex document text width
if Task == "1":

    L = 1.0
    u_list = []
    x_list = []
    t_list = []

    for method in ["Analytic:","FE:", "BE:", "CN:"]:
        for dx in [0.1, 0.01]:
            dt = 0.5*dx*dx
            #Generate t-mesh
            T = int(1.0/dt) #Number of time steps till final time
            t = np.zeros(T)
            for l in range(len(t)):
                t[l] = l*dt
            #Generate x-mesh
            N = int(1.0/dx)   #Number of integration points along x-axis (inner points only)
            x = np.zeros(N+2)
            for k in range(len(x)):
                x[k] = k/(N+1)
            if method == "FE:":
                x_list.append(x)
                t_list.append(t)

            with open (method+str(dx)) as file:
                lines = file.readlines()
                u = np.zeros((len(lines),len(lines[0].split())))
                for i in range(len(lines)):
                    u[i,:] = lines[i].split()

                u_list.append(u)

                fig = plt.figure();
                x,t = np.meshgrid(x,t)


                ax = fig.gca(projection='3d');
                # Plot the surface.
                surf = ax.plot_surface(x, t, u, cmap=cm.coolwarm,
                                   linewidth=0, antialiased=False);
                                   # Customize the z axis.
                #ax.set_zlim(-0.10, 1.40);
                for angle in range(0,230):
                    ax.view_init(40,angle)
                ax.zaxis.set_major_locator(LinearLocator(10));
                ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'));
                plt.xlabel("x")
                plt.ylabel("t")
                name = method+" dx = "+str(dx)
                plt.title(name)
                #fig.savefig("plots/"+name+".png")


    dx = [0.1, 0.01]
    for i in range(2):
        fig = plt.figure();
        fig.set_size_inches(w=w*0.8,h= 4.5)
        plt.title("Computed solutions at time = 0.1 \n dx = %g" % (dx[i]))
        plt.plot(x_list[i], u_list[i][int(len(t_list[i])/10)], ".")
        plt.plot(x_list[i], u_list[2+i][int(len(t_list[i])/10)], ".")
        plt.plot(x_list[i], u_list[4+i][int(len(t_list[i])/10)], ".")
        plt.plot(x_list[i], u_list[6+i][int(len(t_list[i])/10)])
        #plt.plot(x_analytic, u_analytic[int(len(t_analytic)/10)])
        plt.legend(["FE", "BE", "CN", "Analytic"])
        plt.xlabel("x")
        plt.ylabel("u(x,t=0.1)")
        plt.savefig("plots/1dim/Comparison_0.1_"+str(dx[i])+".pgf")


        fig = plt.figure();
        plt.title("Absolute difference between computed and analytical at time = 0.1 \n dx = %g" % (dx[i]))
        plt.plot(x_list[i], abs(u_list[i][int(len(t_list[i])/10)] - u_list[6+i][int(len(t_list[i])/10)]), ".")
        plt.plot(x_list[i], abs(u_list[2+i][int(len(t_list[i])/10)] - u_list[6+i][int(len(t_list[i])/10)]), ".")
        plt.plot(x_list[i], abs(u_list[4+i][int(len(t_list[i])/10)] - u_list[6+i][int(len(t_list[i])/10)]), ".")
        plt.plot(x_list[i], abs(u_list[6+i][int(len(t_list[i])/10)] - u_list[6+i][int(len(t_list[i])/10)]))
        #plt.plot(x_analytic, u_analytic[0])
        plt.legend(["FE", "BE", "CN", "Analytic"])
        plt.xlabel("x")
        plt.ylabel("u(x,t=0.1) - $u_{exact}$(x,t=0.1)")
        #plt.ylim((-10**(-5),10**(-3)))
        plt.savefig("plots/1dim/Differences_0.1_"+str(dx[i])+".pgf")
    plt.show()

    for i in range(2):
        fig = plt.figure();
        plt.title("Computed solutions at time = 0.2 \n dx = %g" % (dx[i]))
        plt.plot(x_list[i], u_list[i][int(len(t_list[i])/5)], ".")
        plt.plot(x_list[i], u_list[2+i][int(len(t_list[i])/5)], ".")
        plt.plot(x_list[i], u_list[4+i][int(len(t_list[i])/5)], ".")
        plt.plot(x_list[i], u_list[6+i][int(len(t_list[i])/5)])
        #plt.plot(x_analytic, u_analytic[int(len(t_analytic)/5)])
        plt.legend(["FE", "BE", "CN", "Analytic"])
        plt.xlabel("x")
        plt.ylabel("u(x,t=0.2)")
        plt.savefig("plots/1dim/Comparison_0.2_"+str(dx[i])+".pgf")

        fig = plt.figure();
        plt.title("Absolute difference between computed and analytical at time = 0.2 \n dx = %g" % (dx[i]))
        plt.plot(x_list[i], abs(u_list[i][int(len(t_list[i])/5)] - u_list[6+i][int(len(t_list[i])/5)]), ".")
        plt.plot(x_list[i], abs(u_list[2+i][int(len(t_list[i])/5)] - u_list[6+i][int(len(t_list[i])/5)]), ".")
        plt.plot(x_list[i], abs(u_list[4+i][int(len(t_list[i])/5)] - u_list[6+i][int(len(t_list[i])/5)]), ".")
        plt.plot(x_list[i], abs(u_list[6+i][int(len(t_list[i])/5)] - u_list[6+i][int(len(t_list[i])/5)]))
        #plt.plot(x_analytic, u_analytic[int(len(t_analytic)/5)])
        plt.legend(["FE", "BE", "CN", "Analytic"])
        plt.xlabel("x")
        plt.ylabel("u(x,t=0.2) - $u_{exact}$(x,t=0.2)")
        #plt.ylim((-10**(-5),10**(-3)))
        plt.savefig("plots/1dim/Differences_0.2_"+str(dx[i])+".pgf")
    plt.show()


elif Task == "2":

    L = 1.0
    u_list = []
    dxlist = [0.1,0.01]

    for method in ["Analytic","Implicit"]:
        for dx in dxlist:
            #2 different dt for analysing stability of scheme
            dtlist = [dx,dx/10,dx/100]
            for dt in dtlist:
                T = int(1.0/dt)
                #Generate t-mesh
                t = np.linspace(0,1,T)
                #Generate x- and y-mesh
                N = int(1.0/dx)
                x = np.linspace(0,1,N+2)
                y = np.linspace(0,1,N+2)
                filename = "2dim_"+method+":dx="+str(dx)+"dt="+str(dt)
                time = int(0.01*T)  #The time we choose to sample the solution at
                with open(filename) as file:

                    lines = file.readlines()
                    u = np.zeros((len(x),len(y)))
                    for i in range(len(y)):
                        data = lines[time*len(x)+i].split()

                        u[i] = data

                    if dt == dx:
                        fig = plt.figure();
                        fig.set_size_inches(w=w*0.8,h= 3.5)
                        x_,y_ = np.meshgrid(x,y)

                        ax = fig.gca(projection='3d',xlim = (0,1.0),ylim = (0,1.0),zlim = (0,1.0));
                        # Plot the surface.
                        surf = ax.plot_surface(x_, y_, u, cmap=cm.coolwarm,
                                           linewidth=0, antialiased=False);
                                           # Customize the z axis.
                        #ax.set_zlim(-0.10, 1.40);
                        for angle in range(0,230):
                            ax.view_init(40,angle)
                        ax.zaxis.set_major_locator(LinearLocator(10));
                        ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'));
                        plt.xlabel("x")
                        plt.ylabel("y")
                        name = "2-dim "+method+": dx = "+str(dx)
                        plt.title(name)
                        fig.savefig("plots/2dim/"+method+"/"+str(dx)+"/"+name.replace(" ","")+".pgf")
                        plt.close()

                        #Printing maximum absolute differences
                    u_list.append(u)
    combinations = len(dxlist)*3
    print("Mean absolute differences between implicit and analytical")
    for i in range(combinations):
        diff = np.mean(abs(u_list[i]-u_list[combinations+i]))

        print(diff)

elif Task == "3":
    w = 5.78851          # Latex document text width
    #fig = plt.figure();
    #fig.set_size_inches(w=w*1.0,h= 4.0)

    # Analytical solution to heat production
    analytic = []
    x = np.linspace(0,120,101)
    x_ = [x[:17],x[17:34],x[34:]]
    z_ = [
        [-0.28,-23.66,8],
        [-0.07,-15.26,92],
        [-0.01,-10.46,188]
        ]
    y2 = []
    for x,zone in zip(x_,z_):
        y = np.polyval(zone,-x)
        y2.append(y)


    analytic.append(np.concatenate((y2[0],y2[1],y2[2])))
    # Analytical solution to heat production + radioactive enrichment
    z_ = [
        [-0.28,-29,8],
        [-0.07,-20.6,92],
        [-0.11,-23.8,28]
        ]
    y3 = []
    for x,zone in zip(x_,z_):
        y = np.polyval(zone,-x)
        y3.append(y)
    analytic.append(np.concatenate((y3[0],y3[1],y3[2])))
    """
    for i in range(len(analytic)):
        # Analytical solution of heat production plot
        plt.plot(np.linspace(0,120,101),analytic[i],"--")
    plt.xlabel("Depth [km]")
    plt.ylabel(r"Temperature $[^\circ C]$")
    plt.grid()
    plt.legend(["Heat","Enriched mantle $\\&$ No Decay"])
    plt.title("Analytical")
    plt.savefig("plots/Lithosphere/Analytical.pgf")
    plt.show()
    """

    #fig = plt.figure();
    #fig.set_size_inches(w=w*1.0,h= 4.0)
    Nx = 126
    Ny = 101
    u = np.zeros((Nx,Ny))
    numerical = []
    for filename in ["Heat","No_Decay", "Decay"]:
        dx = 0.01
        dt = dx
        T = int(1.0/dt)
        #Generate t-mesh
        t = np.linspace(0,1,T)
        #Generate x- and y-mesh
        N = int(1.0/dx)
        with open(filename) as file:
            lines = file.readlines()
            for t in tqdm(range(T)):
                for i in range(Nx):
                    data = lines[t*Nx+i].split()

                    u[i] = data


        temp = u[int(Nx/2)]*1292 + 8
        depth = np.linspace(0,120,Ny)
        numerical.append(temp)
        # Numerical simulation of heat production plot
        """
        plt.plot(depth,temp)
        plt.xlabel("Depth [km]")
        plt.ylabel(r"Temperature $[^\circ C]$")
        """

    """
    plt.legend(["Heat", "Enriched mantle $\\&$ No Decay"])
    plt.title("Numerical")
    plt.grid()
    plt.savefig("plots/Lithosphere/Numerical.pgf")
    plt.show()
    """

    # Numerical vs Analytical for case 1 & 2 plot
    fig = plt.figure();
    fig.set_size_inches(w=w*1.0,h= 4.0)
    x = np.linspace(0,120,101)
    for i in range(len(analytic)):
        a = np.array(analytic[i])
        n = np.array(numerical[i])
        plt.plot(x,a)
        plt.plot(x,n,"--")
    plt.xlabel('Depth [km]')
    plt.ylabel(r"Temperature $[^\circ C]$")
    plt.grid()
    plt.legend(["Analytical: Heat","Numerical: Heat","Analytical: Enriched mantle $\\&$ No Decay","Numerical: Enriched mantle $\\&$ No Decay"])
    plt.savefig("plots/Lithosphere/Comparison.pgf")
    plt.show()


    # Relative error plot
    fig = plt.figure();
    fig.set_size_inches(w=w*1.0,h= 4.0)
    x = np.linspace(0,120,101)
    for i in range(len(analytic)):
        a = np.array(analytic[i])
        n = np.array(numerical[i])
        z = abs((a-n)/(a))
        plt.plot(x[1:-1],z[1:-1])
        plt.yscale('log')
    plt.xlabel('Depth [km]')
    plt.ylabel('Relative error')
    plt.grid()
    plt.legend(["Heat","Enriched mantle $\\&$ No Decay"])
    plt.savefig("plots/Lithosphere/Relative_Error.pgf")
    plt.show()


    # Comparison between Numerical simulation of Decay and No Decay plot
    fig = plt.figure();
    fig.set_size_inches(w=w*1.0,h= 4.0)
    x = np.linspace(0,120,101)
    for i in range(1,len(numerical)):
        n = np.array(numerical[i])
        plt.plot(x,n)
    plt.xlabel('Depth [km]')
    plt.ylabel(r"Temperature $[^\circ C]$")
    plt.grid()
    plt.title("Numerical Simulation")
    plt.legend(["Enriched mantle $\\&$ No Decay", "Enriched mantle $\\&$ Decay"])
    plt.savefig("plots/Lithosphere/Decay_NoDecay.pgf")
    plt.show()

    #Plots the difference at the end of the simulation with and without decay
    nodecay = np.zeros((Nx,Ny))
    decay = np.zeros((Nx,Ny))
    for filename in ["No_Decay", "Decay"]:
        dx = 0.01
        dt = dx
        T = int(1.0/dt)
        #Generate t-mesh
        t = np.linspace(0,1,T)
        #Generate x- and y-mesh
        N = int(1.0/dx)
        with open(filename) as file:
            lines = file.readlines()
            for t in tqdm(range(T)):
                for i in range(Nx):
                    data = lines[t*Nx+i].split()
                    if filename == "No_Decay":
                        nodecay[i] = data

                    elif filename == "Decay":
                        decay[i] = data

    decay = decay*1292 + 8 # Scaling [0,1] --> [8,1300]
    nodecay = nodecay*1292 + 8 # Scaling [0,1] --> [8,1300]
    diff = (nodecay - decay)
    x = np.linspace(0,150,Nx)
    y = np.linspace(0,120,Ny)
    levels = np.linspace(0,35,100)
    plt.contourf(x,y,diff.T,levels=levels)
    cbar = plt.colorbar(ticks=[0,7,14,21,28,35])
    cbar.set_label(r"$T_{diff} [^\circ C]$")
    plt.xlabel("Width [km]")
    plt.ylabel("Depth [km]")
    plt.savefig("plots/Lithosphere/Tdiff.pgf")
    plt.show()

    maxtemp = np.max(diff.T)
    depth = np.argmax(diff.T[:,int(Nx/2)],axis=0)
    width = np.argmax(diff.T[int(Ny/2),:],axis=0)
    print("Max Temperature Difference =  %g, happens at Depth = %g & Width = %g" % (maxtemp, depth*1.2, width*1.2))


#, "Enriched mantle $\\&$ Decay"
else:
    print("Please write either 1, 2 or 3")
