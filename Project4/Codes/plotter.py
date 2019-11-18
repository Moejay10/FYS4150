import matplotlib.pyplot as plt
import numpy as np
import os
print("Which Project Task do you want to run")
print("Task C - Equilibrium State: Write c")
print("Task D - Probability Histogram: Write d")
print("Task E & F - Phase Transitions & Critical Temperature: Write e")
Task = input("Write here: ")

tablesdir = os.path.join(os.path.dirname(__file__), '..', 'Tables')
plotsdir = os.path.join(os.path.dirname(__file__), '..', 'Results/Plots')
"""
-------------
Equilibrium
-------------
"""

if Task == "c":
    filenames = ["Ordered1","Ordered"]
    filenames2 = ["Unordered1", "Unordered"]

    MCcycles = []
    energyO = []
    energyO2 = []
    energyU = []
    energyU2 = []
    magO = []
    magO2 = []
    magU = []
    magU2 = []
    NconfigsU = []
    NconfigsU2 = []
    for i in (filenames):
        with open(os.path.join(tablesdir,i)) as file:
            lines = file.readlines()
            #Skip the first two lines
            for j in range(2,len(lines)):
                line = lines[j]
                pieces = line.split()
                if i == "Ordered1":
                    MCcycles.append(float(pieces[0]))
                    energyO.append(float(pieces[1]))
                    magO.append(float(pieces[2]))
                else:
                    energyO2.append(float(pieces[1]))
                    magO2.append(float(pieces[2]))

    for i in (filenames2):
        with open(os.path.join(tablesdir,i)) as file:
            lines = file.readlines()
            #Skip the first two lines
            for j in range(2,len(lines)):
                line = lines[j]
                pieces = line.split()
                if i == "Unordered1":
                    energyU.append(float(pieces[1]))
                    magU.append(float(pieces[2]))
                    NconfigsU.append(float(pieces[3]))

                else:
                    energyU2.append(float(pieces[1]))
                    magU2.append(float(pieces[2]))
                    NconfigsU2.append(float(pieces[3]))
    plt.figure()
    plt.title("Ordered")
    plt.plot(MCcycles, energyO)
    plt.plot(MCcycles,energyO2)
    plt.legend(["T = 1.0","T = 2.4"])
    plt.xlabel("# of Monte Carlo cycles")
    plt.ylabel("Energy expectation value $\langle$E$\\rangle$ [J]")
    plt.savefig(os.path.join(plotsdir,"Energy_exp_ordered.png"))

    plt.figure()
    plt.title("Unordered")
    plt.plot(MCcycles, energyU)
    plt.plot(MCcycles,energyU2)
    plt.legend(["T = 1.0","T = 2.4"])
    plt.xlabel("# of Monte Carlo cycles")
    plt.ylabel("Energy expectation value $\langle$E$\\rangle$ [J]")
    plt.savefig(os.path.join(plotsdir,"Energy_exp_unordered.png"))

    plt.figure()
    plt.title("Ordered")
    plt.plot(MCcycles, magO, "")
    plt.plot(MCcycles, magO2, "")
    plt.legend(["T = 1.0","T = 2.4"])
    plt.xlabel("# of Monte Carlo cycles")
    plt.ylabel("Magnetization expectation value $\langle$|M|$\\rangle$ [1]")
    plt.savefig(os.path.join(plotsdir,"Magn_exp_ordered.png"))


    plt.figure()
    plt.title("Unordered")
    plt.plot(MCcycles, magU, "")
    plt.plot(MCcycles, magU2, "")
    plt.legend(["T = 1.0","T = 2.4"])
    plt.xlabel("# of Monte Carlo cycles")
    plt.ylabel("Magnetization expectation value $\langle$|M|$\\rangle$ [1]")
    plt.savefig(os.path.join(plotsdir,"Magn_exp_unordered.png"))

    plt.figure()
    plt.title("Unordered")
    plt.plot(MCcycles, NconfigsU, "")
    plt.plot(MCcycles, NconfigsU2, "")
    plt.legend(["T = 1.0","T = 2.4"])
    plt.xlabel("# of Monte Carlo cycles")
    plt.ylabel("Accepted configurations (normalized)")
    plt.savefig(os.path.join(plotsdir,"Accepted_configs_unordered.png"))

    Temp = []
    configs = []
    with open(os.path.join(tablesdir,"Nconfig_vs_Temp")) as file:
        lines = file.readlines()
        for i in range(2,len(lines)):

            pieces = lines[i].split()
            Temp.append(float(pieces[0]))
            configs.append(float(pieces[1]))
    plt.figure()
    plt.plot(Temp,configs)
    plt.xlabel("Temperature [kT/J]")
    plt.ylabel("Accepted number of configurations (normalized)")
    plt.title("Accepted number of configurations (normalized) as a function of T")
    plt.savefig(os.path.join(plotsdir,"Accepted_configs_temperature.png"))
    plt.show()

"""
-------------
Probabilities
-------------
"""

if Task == "d":
    filenames = ["Probability_1","Probability_24"]


    for i in filenames:
        with open(os.path.join(tablesdir,i)) as file:
            lines = file.readlines()
        Energies = []
        counts = []
        max_count = 0
        most_probable_energy = 0
        for j in range(1,len(lines)):
            line = lines[j]
            pieces = line.split()
            energy = float(pieces[0])
            count = float(pieces[1])
            Energies.append((energy))
            counts.append((count))
            if count > max_count:
                max_count = count
                most_probable_energy = energy
        plt.bar(Energies,counts,width = 4 if i == "Probability_1" else 3)
        plt.xlim(-805,-770) if i == "Probability_1" else plt.xlim(-705,-305)
        plt.xlabel("Energy [J]")
        plt.ylabel("Energy counts")
        plt.tight_layout()
        plt.subplots_adjust(top=0.88)

        if i == "Probability_1":
            plt.title("T = 1.0")
        else:
            plt.title("T = 2.4")
        props = dict(boxstyle='round', facecolor='wheat', alpha=1)
        plt.text(0.05*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0] ,plt.ylim()[1]*0.85, "Most probable energy:\n" + str(most_probable_energy), bbox = props)
        plt.savefig(os.path.join(plotsdir,i+".png"))
        plt.show()

if Task == "e":
    with open(os.path.join(tablesdir,"Temperature_100")) as file:
        lines = file.readlines()
    temps = []
    energylist = []
    maglist = []
    Cvlist = []
    Suscplist = []
    indeks = 0
    for i in range(1, len(lines)):
        pieces = lines[i].split()
        temps.append(float(pieces[0]))
        energylist.append(float(pieces[1]))
        maglist.append(float(pieces[2]))
        Cvlist.append(float(pieces[3]))
        Suscplist.append(float(pieces[4]))

    firstTemp = temps[0]
    for i in range(1,len(temps)):
        if temps[i] == firstTemp:
            temps = temps[0:i]
            break

    TCCv = []
    TCX = []
    for i in range(int(len(energylist)/len(temps))):
        max_temp = 0
        sublistCv = Cvlist[i*len(temps):len(temps)*(i+1)]
        sublistSuscp = Suscplist[i*len(temps):len(temps)*(i+1)]
        maxCv = max(sublistCv)
        maxSuscp = max(sublistSuscp)
        TCCv.append(temps[sublistCv.index(maxCv)])
        TCX.append(temps[sublistSuscp.index(maxSuscp)])
        print("Tc for Cv =",temps[sublistCv.index(maxCv)])
        print("Tc for X =",temps[sublistSuscp.index(maxSuscp)])

    plt.figure()
    plt.title("Mean Energy")
    plt.xlabel("T [kT/J]")
    plt.ylabel("Energy expectation value $\langle$E$\\rangle$ [J]")
    for i in range(int(len(energylist)/len(temps))):
        plt.plot(temps,energylist[i*len(temps):len(temps)*(i+1)],"")
    plt.legend(["L = 40","L = 60","L = 80","L = 100"])
    plt.savefig(os.path.join(plotsdir,"Phase_trans_energy.png"))

    plt.figure()
    plt.title("Absolute mean Magnetization")
    plt.xlabel("T [kT/J]")
    plt.ylabel("Magnetization expectation value $\langle$|M|$\\rangle$ [1]")
    for i in range(int(len(energylist)/len(temps))):
        plt.plot(temps,maglist[i*len(temps):len(temps)*(i+1)],"")
    plt.legend(["L = 40","L = 60","L = 80","L = 100"])
    plt.savefig(os.path.join(plotsdir,"Phase_trans_mag.png"))


    plt.figure()
    plt.title("Specific heat")
    plt.xlabel("T [kT/J]")
    plt.ylabel("Specific heat $\langle$$C_v$$\\rangle$ [$J^2/kT^2$]")
    for i in range(int(len(energylist)/len(temps))):
        plt.plot(temps,Cvlist[i*len(temps):len(temps)*(i+1)],"")
    plt.legend(["L = 40","L = 60","L = 80","L = 100"])
    plt.savefig(os.path.join(plotsdir,"Phase_trans_Cv.png"))


    plt.figure()
    plt.title("Susceptibility")
    plt.xlabel("T [kT/J]")
    plt.ylabel("Susceptibility $\langle$$\chi$$\\rangle$ [1/kT]")
    for i in range(int(len(energylist)/len(temps))):
        plt.plot(temps,Suscplist[i*len(temps):len(temps)*(i+1)],"")
    plt.legend(["L = 40","L = 60","L = 80","L = 100"])
    plt.savefig(os.path.join(plotsdir,"Phase_trans_suscp.png"))

    plt.show()

    """
    Task f)
    """
    #Performing a linear regression to find critical temp in thermodyn. limit
    TCCv = np.array(TCCv)
    TCX = np.array(TCX)
    Llist = np.array([40,60,80,100])
    Llist = 1.0/Llist

    linreg1 = np.polyfit(Llist,TCCv,1)
    linreg2 = np.polyfit(Llist,TCX,1)

    plt.figure()
    plt.title("Specific heat $C_V$")
    plt.xlabel("$\\frac{1}{L}$")
    plt.ylabel("$T_C$ [kT/J]")
    plt.plot(Llist,TCCv,"o")
    plt.plot(Llist,np.polyval(linreg1,Llist))
    plt.legend(["$T_C$(L) from simulations","$T_C(L)$ = a$\\cdot$ $\\frac{1}{L}$ + $T_C(L = \infty)$ $\\to$ %g$\\cdot$x + %g" % (linreg1[0],linreg1[1])])
    plt.savefig(os.path.join(plotsdir,"linregCv.png"))

    plt.figure()
    plt.title("Susceptibility $\chi$")
    plt.xlabel("$\\frac{1}{L}$")
    plt.ylabel("$T_C$ [kT/J]")
    plt.plot(Llist,TCX,"o")
    plt.plot(Llist,np.polyval(linreg2,Llist))
    plt.legend(["$T_C$(L) from simulations","$T_C(L)$ = a$\\cdot$ $\\frac{1}{L}$ + $T_C(L = \infty)$ $\\to$ %g$\\cdot$x + %g" % (linreg2[0],linreg2[1])])
    plt.savefig(os.path.join(plotsdir,"linregX.png"))

    print("\n")
    print("The estimated Critical Temperature from our simulations is Tc = %g " % (0.5*(linreg1[1]+linreg2[1])))

    plt.show()
