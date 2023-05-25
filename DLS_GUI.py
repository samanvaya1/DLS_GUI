#Import the tkinter package to create the GUI and Initialize the GUIfrom tkinter import  *
import tkinter
from tkinter import *
from tkinter import filedialog, ttk
from tkinter import messagebox
from tkinter import ttk

#Import Package to allow for the finding of files
import os

#Initialize the GUI
root = Tk()
root.title('DLS GUI')
#root.configure(bg='white')
root.geometry("1050x725")#(x,y)


#Import Maplotlib to plot data
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
matplotlib.use("TkAgg")
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,NavigationToolbar2Tk)

#Import Statistic packages for the Cumulant Methods
import statsmodels.sandbox.distributions.extras as extras

#Import numpy and math to help with the math
import numpy as np
import math

#Imports for Packages for the Analysis Algorthuims (This is done in the main file to help in the freezing process to generate the executible)
from scipy.optimize import curve_fit
from scipy import optimize

#For timing of events
import time

#Import the Analysis Algorthuims
import DLS_linear_fit
import DLS_quad_fit
import DLS_cubic_fit
import DLS_quartic_fit
import DLS_contin
import DLS_REPES
import DLS_DYNALS
import DLS_NNLS

#Initialize the tau column and g2 column of the data
global tau
global g2
global TimeChange
tau = []
g2 = []
TimeChange = False
#Function that checks if it is normalized and if it is not it normalizes the data

#Function that checks if it is normalized and if it is not it normalizes the data
#If normalized = 1 If not normalized = 0
def Normalize_Check(tauLC,g2LC,B,beta):
    #Considered normalize if the baseline is less than 0.1
    if (select_dtype.get() == "g1\u00b2"):

        #Double Check to Ensure System is normalized with baseline at zero
        g2Norm = np.subtract(np.divide(g2LC,B),1)
        beta = beta/B-1
        g1_squared = np.divide(g2Norm, beta)
        Normalized = 0
        print("Works")

    #If the data is gw use: C(tau) = (G2(tau)-1-B)/B to normalize the data and account for (g2-1)/beta
    else:
        g2Norm = np.subtract(np.divide(np.subtract(g2LC,1),B),1)
        beta = beta/B-1
        g1_squared = np.divide(g2Norm, beta)
        Normalized = 0
    return g2Norm , Normalized,beta,g1_squared

#Generate a modified boltzman distibution based on the coefficiients provided based of a Gram-Charlier expansion on a Gaussian distribution found on github
def Distrubution(mean,sigma,skew,kurt,min,max):
    #This is done to ensure that the average is always included break the dist into two halfs meeting at the mean
    points = int(Entry_points.get())
    half = int(points/2)

    #If Log Scale is desired
    if(select_dist_type.get() == 'Logarithmic'):
        x0L = np.logspace(np.log10(max), np.log10(mean), half)
        x0R = np.logspace(np.log10(mean), np.log10(min), half)
        Plot_prob.set_xscale('log')
    #If Linear Scale is desired
    else:
        x0L = np.linspace(max,mean,half)
        x0R = np.linspace(mean,min,half)
        Plot_prob.set_xscale('linear')
    x0 = np.concatenate((x0L,x0R))

    #Get the distribution data and normalize it to sum up to 1
    pdf = extras.pdf_mvsk([mean,sigma,skew,kurt]) 
    y0 = np.abs(pdf(x0))
    y0 = np.divide(y0,np.sum((y0)))
    return x0 , y0

#Load data and get the plots G2 plot up and running
def LoadData():
    #Initialize global varibles
    global filename,Figure_g2, count_update,tau,g2,beta,Plot_g2,TimeChange
    #Case of a time change
    #Consider the case already a plot need to fix the axes of the plot
    if(count_update > 0):
        #Remove all the current information on the G2 plot
        count_update = 0
        Clear()
        Entry_low.delete(0,'end')
        Entry_beta.delete(0,'end')
        Entry_high.delete(0,'end')
        len_g2 = len(Plot_g2.lines)
        for i in range(len_g2):
            del Plot_g2.lines[0]
        tau  = []
        g2 = []

        #Updated figure displayed on the GUI to show the cleared plot
        canvas_g2 = FigureCanvasTkAgg(Figure_g2, root)
        canvas_g2.get_tk_widget().place(height=400, width=500, x=510, y=50)
        Plot_g2.legend().remove()
        Plot_g2.legend()
        Plot_g2.plot(tau, g2, color='red')

    #Check If Switching Time
    if(TimeChange):
        TimeChange = False

    else:#Ask user for the file
        root.filename = filedialog.askopenfilename(initialdir= "/" , title = "Select A File",filetypes=(("txt files","*.txt"),("all files","*.*")))
        filename = root.filename
    #Check if a file is selected
    if(filename != "" ):
        File = 1
        #label_filename = Label(root,text = str(filename))
        #label_filename.place(height=30, width=100, x=110, y=0)
        filenameshort = os.path.basename(filename)
        button_data.configure(text=filenameshort)
        file = open(filename,"r")
        file_data = np.loadtxt(filename,usecols=(0,1))

        if(select_time_type.get() == "s"):
            TimeFactor = 1
        else:
            TimeFactor = 10**(-6)

        #Import the data into the global vectors preallocated
        tau = file_data[:,0]*TimeFactor
        g2 = file_data[:,1]
        file.close()
        #Find some default taus for the baseline calculation ,betas
        length = len(tau)
        Entry_low.insert(0,str(float('%.3g' %tau[length-22])))
        Entry_high.insert(0,str(float('%.3g' %tau[length-2])))
        #get the beta value
        B = 0
        index = 0
        # Go until lower limit is reached and record only the given data
        count = len(tau)-1
        tau_low = tau[length-22]
        tau_high = tau[length-2]
        while (tau[count] > tau_low):
            if (tau[count] < tau_high):
                B = g2[count] + B
                index = index + 1
            count = count - 1
        B = B / index
        beta = 0.5*(g2[0]+g2[1])
        Entry_beta.insert(0,str("{:.3e}".format(beta)))
        #With the imported data input the data into G2 graph using red open circles
        #delete the place holder graph
        del Plot_g2.lines[0]
        label = 'Data'
        Plot_G2(tau,g2,'ro',label)
        Update()

#Multi-Angle code Requires a folder with only the textfiles from a multiangle experiment. This requires that the experiment be conducted in the same conditions as listed in teh main GUI
#Need multiple data sets to ensure that Multiangle code is actually functional and somehow is broken?
def MultiAngle():
    #Load Multiple Data Sets
    #Creates a second window to input the Multiangle Data
    AngleWindow =  Toplevel(root)
    AngleWindow.title('DLS GUI Multi Angle')
    AngleWindow.geometry("600x200")
    #FUnction to handle the multiple data files that are needed
    def LoadMulti():
        global G2_Multi
        global Tau_Multi
        global x1
        global y1
        global NumFiles
        global filesMulti
        global File_Names_Long
        G2_Multi = dict()
        Tau_Multi = dict()
        count = 0
        Num = str(count)
        Time = Drop_time_type.get()
        if(Time == "s"):
            TimeFactor = 1
        else:
            TimeFactor = 10**(-6)

        #Parameters for placement
        h=20; w =200; x1 = 140 ; y1 = 0;
        root.directory = filedialog.askdirectory()
        path = root.directory
        filesMulti = os.listdir(path)
        if os.path.isfile(path+'/.DS_Store'):
            filesMulti.remove('.DS_Store')
        NumFiles = len(filesMulti)
        File_Names_Long = []
        #Loop through the files to gather the G2 and Tau data
        for index in range(NumFiles):
            button_data.configure(text=filesMulti[index])
            file = open(path+"/"+filesMulti[index], "r")
            File_Names_Long.append(path+"/"+filesMulti[index])
            file_data = np.loadtxt(file, usecols=(0, 1))
            file.close()
            Tau_Multi[Num] = file_data[:, 0] * TimeFactor
            G2_Multi[Num] = file_data[:, 1]
            Entry_FileNames = Entry(AngleWindow)
            Entry_FileNames.place(height=h, width=w, x=x1, y=y1)
            Entry_FileNames.insert(0, filesMulti[index])
            Entry_FileNames.config(state='readonly')

            #Entry for multiple angles and update the main GUI with plots of the data
            Entry_MultiAngle.append(Entry(AngleWindow))
            Entry_MultiAngle[index].place(height=h,width=.25*w,x=x1+1.5*w,y=y1)
            Label_MultiAngle = Label(AngleWindow,text = "Angle ("+u"\N{DEGREE SIGN}"+'):')
            Label_MultiAngle.place(height = h , width = w*.5,x=x1+w,y=y1)
            label = 'Data'

            Plot_G2(Tau_Multi[Num], G2_Multi[Num], 'ro', label)
            y1 = y1 + 20
            #Keep track of the number of files used
            count = count + 1
            Num = str(count)
    #This function is used to get the size distribution from the data. It is an
    def RunMulti():
        global tau, error, r_dist, x, Fit_data, gamma_prob, beta, B
        global filename,filesMulti,File_Names_Long
        #Get the dist type
        if (select_dist_type.get() == "Logarithmic"):
            dist = 0
        else:
            dist = 1
        if (Entry_regularization.get() == ""):
            Entry_regularization.insert(0, '0.5')
        global wavelength
        wavelength = int(Entry_wavelength.get()) * 10 ** (-9)
        # refractive
        global refractive
        refractive = float(Entry_refractive.get())
        # Viscosity
        global viscosity
        viscosity = float(Entry_viscosity.get())
        # Temperant
        global temperature
        temperature = float(Entry_temperant.get()) + 273
        # Boltzman Constant
        global boltzman_con
        boltzman_con = 1.38064852 * 10 ** (-23)
        Fit_Multi = dict()
        #Find the Initial Guess depending on Method
        if (select_methods.get() == "Quad Fit"):
            # Initial guess
            initialGuess = [1, 1]
        elif (select_methods.get() == 'Cubic Fit'):
            # Initial guess
            initialGuess = [1, 1, 1]
        elif (select_methods.get() == 'Quartic Fit'):
            # Initial guess
            initialGuess = [1, 1, 1, 1]
        #need the gamma dist set up first need to redo figure out later
        #Get the intial guess and rmin and rmax ranges
        elif (select_methods.get() == "DYNALS" or select_methods.get() == "Contin" or select_methods.get() == "NNLS" or select_methods.get()=="Choose..."):
            # Initial Conditions
            g2_holder = G2_Multi['0']
            tau_holder = Tau_Multi['0']
            Angle  = int(Entry_MultiAngle[0].get())
            #Baseline set up
            B = 0
            length = len(g2_holder)
            for index in range(20):
                B = B + g2_holder[length-2-index]
            B = B/20
            beta = (g2_holder[0] + g2_holder[1])/2
            g2Norm, normal, beta, g1_squared = Normalize_Check(tau_holder, g2_holder, B, beta)
            #g1_squared = g2Norm#np.divide(g2Norm, beta)
            n_pts = int(Entry_points.get())
            initialGuess  = [1,1,1,1]
            Params1, Fit = DLS_quartic_fit.DLS_quartic_fit(tau_holder, g2Norm, B, beta,initialGuess)
            q = (4 * np.pi * refractive / wavelength) * np.sin(np.deg2rad( int(Angle) / 2))

            Diffusive_Coef = Params1[0] / (np.square(q))
            r_mean = (boltzman_con * temperature) / (6 * np.pi * viscosity * Diffusive_Coef) * 10 ** 9

            if (select_dist_type.get() == 'Logarithmic'):
                r_min = r_mean / 100
                r_max = r_mean * 100
                r_dist = np.logspace(np.log10(r_min), np.log10(r_max), n_pts)
                Plot_prob.set_xscale('log')
                dist =0
            else:
                r_min = r_mean / 10
                r_max = r_mean * 10
                r_dist = np.linspace(r_min, r_max, n_pts)
                dist =1
            gamma = np.multiply(np.divide(1, r_dist), (16 * boltzman_con * temperature * np.pi * (refractive ** 2)*((np.sin(np.deg2rad(int(Angle) / 2))) ** 2)/(6 * viscosity * wavelength * wavelength)) * 10 ** 9)
            W = np.diag(np.ones(len(g1_squared)))
            # Initial Guess
            x0_in = np.ones(len(gamma)) + (10 ** (-3)) * np.random.rand(1, len(gamma))
            x0 = np.divide(x0_in, np.sum(x0_in))

        elif (select_methods.get() == "REPES"):
            #Initial Conditions
            # Initial Conditions
            g2_holder = G2_Multi['1']
            tau_holder = Tau_Multi['1']
            Angle = int(Entry_MultiAngle[0].get())
            # Baseline set up
            B = 0
            length = len(g2_holder)
            for index in range(20):
                B = B + g2_holder[length - 2 - index]
            B = B / 20
            beta = (g2_holder[0] + g2_holder[1]) / 2
            g2Norm, normal, beta,g1_squared = Normalize_Check(tau_holder, g2_holder, B, beta)
            n_pts = int(Entry_points.get())
            initialGuess  = [1,1,1,1]
            Params1, Fit = DLS_quartic_fit.DLS_quartic_fit(tau_holder, g2Norm, B, beta,initialGuess)
            q = (4 * np.pi * refractive / wavelength) * np.sin(np.deg2rad( int(Angle) / 2))

            Diffusive_Coef = Params1[0] / (np.square(q))
            r_mean = (boltzman_con * temperature) / (6 * np.pi * viscosity * Diffusive_Coef) * 10 ** 9
            r_min = r_mean/100
            r_max = r_mean*100
            n_pts = int(Entry_points.get())
            a0 = np.ones(n_pts)
            a0 = np.divide(a0, np.sum(a0))
        #Cycle through and get the data
        for index in range(NumFiles):
            #Get all the fit data points and plot them on graphs accordingly
            #Get new q
            Angle = int(Entry_MultiAngle[index].get())
            q = (4 * np.pi * refractive / wavelength) * np.sin(np.deg2rad( Angle / 2))
            #Get new baseline and beta
            g2_holder = G2_Multi[str(index)]
            tau_holder = Tau_Multi[str(index)]
            beta_holder = (g2_holder[0]+g2_holder[1])/2
            B_holder =0
            length_holder = len(g2_holder)
            for i in range(10):
                B_holder = B_holder+g2_holder[length_holder-2-i]
            B_holder = B_holder/10

            # Check if the data is linearized
            g2Norm_holder, normal, beta_holder,g1_squared = Normalize_Check(tau_holder, g2_holder, B_holder, beta_holder)
            # Use the selected data method
            if (select_methods.get() == "Linear Fit"):
                g2Norm_holder = np.multiply(g2Norm_holder, beta_holder)

                Fit_data, Fit_Param = DLS_linear_fit.DLS_linear_fit(tau_holder, g2Norm_holder, B_holder, beta_holder)
                Diffusive_Coef = Fit_Param[0] / (np.square(q))
                color = '#00' + str(index) + '000'
                label = 'Data Set ' + str(index)
            elif (select_methods.get() == "Quad Fit"):
                g2Norm_holder = np.multiply(g2Norm_holder, beta_holder)
                Fit_Param, Fit_data = DLS_quad_fit.DLS_quad_fit(tau_holder, g2Norm_holder, B_holder, beta_holder, initialGuess)
                initialGuess = Fit_Param
            elif (select_methods.get() == 'Cubic Fit'):
                g2Norm_holder = np.multiply(g2Norm_holder, beta_holder)
                Fit_Param, Fit_data = DLS_cubic_fit.DLS_cubic_fit(tau_holder, g2Norm_holder, B_holder, beta_holder, initialGuess)
                initialGuess = Fit_Param
            elif (select_methods.get() == 'Quartic Fit'):
                g2Norm_holder = np.multiply(g2Norm_holder, beta_holder)
                Fit_Param, Fit_data = DLS_quartic_fit.DLS_quartic_fit(tau_holder, g2Norm_holder, B_holder, beta_holder, initialGuess)
                initialGuess = Fit_Param
            elif (select_methods.get() == 'Contin'or select_methods.get()=="Choose..."):
                g1_squared_holder = g2Norm_holder
                reg_param = float(Entry_regularization.get())
                n_pts = int(Entry_points.get())
                if (select_dist_type.get() == 'Logarithmic'):
                    r_dist = np.logspace(np.log10(r_min), np.log10(r_max), n_pts)
                    Plot_prob.set_xscale('log')
                else:
                    r_dist = np.linspace(r_min, r_max, n_pts)

                gamma = np.multiply(np.divide(1, r_dist), (16 * boltzman_con * temperature * np.pi * (refractive ** 2) * ((np.sin(np.deg2rad(Angle / 2))) ** 2) / (6 * viscosity * wavelength * wavelength)) * 10 ** 9)
                W = np.diag(np.ones(len(g1_squared_holder)))
                Fit_data, gamma, x = DLS_contin.DLS_contin(tau_holder, g1_squared_holder, gamma, reg_param, W, x0)
                Fit_data = Fit_data*beta_holder
                x0 = x
                x = x / np.sum(x)
                color = '#00' + str(index) + '000'
                label = 'Data Set ' + str(index)
            elif (select_methods.get() == "NNLS"):
                g1_squared_holder = g2Norm_holder
                reg_param = 0
                n_pts = int(Entry_points.get())
                if (select_dist_type.get() == 'Logarithmic'):
                    r_dist = np.logspace(np.log10(r_min), np.log10(r_max), n_pts)
                    Plot_prob.set_xscale('log')
                    dist = 0
                else:
                    r_dist = np.linspace(r_min, r_max, n_pts)
                    dist = 1
                gamma = np.multiply(np.divide(1, r_dist), (
                            16 * boltzman_con * temperature * np.pi * (refractive ** 2) * (
                                (np.sin(np.deg2rad(Angle / 2))) ** 2) / (
                                        6 * viscosity * wavelength * wavelength)) * 10 ** 9)
                W = np.diag(np.ones(len(g1_squared)))
                Fit_data, gamma, x = DLS_contin.DLS_contin(tau_holder, g1_squared_holder, gamma, reg_param, W, x0)
                Fit_data = Fit_data*beta_holder
                x0 = x
                x = x / np.sum(x)
                color = '#00' + str(index) + '000'
                label = 'Data Set ' + str(index)
                Plot_Prob(r_dist, x, color, label)
                gamma_prob = x

            elif (select_methods.get() == "REPES"):
                g1_squared_holder = g2Norm_holder
                g1_squared_holder = [0 if i < 0 else i for i in g1_squared_holder]
                g1 = g1_squared_holder#np.sqrt(g1_squared_holder)
                reg_param = float(Entry_regularization.get())
                n_pts = int(Entry_points.get())
                if (select_dist_type.get() == 'Logarithmic'):
                    r_dist = np.logspace(np.log10(r_min), np.log10(r_max), n_pts)
                    Plot_prob.set_xscale('log')
                    dist =0
                else:
                    r_dist = np.linspace(r_min, r_max, n_pts)
                    dist =1
                gamma = np.multiply(np.divide(1, r_dist), (
                            16 * boltzman_con * temperature * np.pi * (refractive ** 2) * (
                                (np.sin(np.deg2rad(Angle / 2))) ** 2) / (
                                        6 * viscosity * wavelength * wavelength)) * 10 ** 9)
                Tau = np.divide(1, gamma)
                Fit_data, gamma, x = DLS_REPES.DLS_REPES(tau_holder, g1, Tau, reg_param, a0)
                Fit_data = np.multiply(np.square(Fit_data), beta_holder)
                color = '#00' + str(index) + '000'
                label = 'Data Set ' + str(index)
                Plot_Prob(r_dist, x, color, label)
                gamma_prob = x

            elif (select_methods.get() == "DYNALS"):
                g1_squared = g2Norm_holder
                g1_squared = [0 if i < 0 else i for i in g1_squared]
                g1 = np.sqrt(g1_squared)
                n_pts = int(Entry_points.get())
                if (select_dist_type.get() == 'Logarithmic'):
                    r_dist = np.logspace(np.log10(r_min), np.log10(r_max), n_pts)
                    Plot_prob.set_xscale('log')
                    dist = 0
                else:
                    r_dist = np.linspace(r_min, r_max, n_pts)
                    dist = 1
                gamma = np.multiply(np.divide(1, r_dist), (
                            16 * boltzman_con * temperature * np.pi * (refractive ** 2) * (
                                (np.sin(np.deg2rad(Angle / 2))) ** 2) / (
                                        6 * viscosity * wavelength * wavelength)) * 10 ** 9)
                #Fix the DYNALS initial guess bug
                Fit_data, gamma, x = DLS_DYNALS.DLS_DYNALS(tau_holder, g1, gamma, x0)
                Fit_data = beta_holder*(Fit_data**2)
                x0 = x
                color = '#00' + str(index) + '000'
                label = 'Data Set ' + str(index)
                Plot_Prob(r_dist,x,color,label)
                gamma_prob = x

            #plotting the data
            if(select_methods.get() == "Linear Fit" or select_methods.get() == "Quad Fit" or select_methods.get() == 'Cubic Fit' or select_methods.get() == 'Quartic Fit'):
                # Get the diffusive coeff and average radius
                Diffusive_Coef = Fit_Param[0] / (np.square(q))
                r_mean = (boltzman_con * temperature) / (6 * np.pi * viscosity * Diffusive_Coef) * 10 ** 9
                r_min = r_mean/100
                r_max = r_mean*100
                if (select_methods.get() == "Linear Fit"):
                    color = '#00'+ str(index) +'000'
                    label = 'Data Set ' + str(index)
                elif (select_methods.get() == "Quad Fit" or select_methods.get() == "Cubic Fit" or select_methods.get() == 'Quartic Fit'):
                    # Get upper and lower bounds
                    D_min = boltzman_con * temperature / (6 * np.pi * viscosity * r_min)
                    Gamma_min = D_min * (q * q) * (10 ** 9)
                    D_max = boltzman_con * temperature / (6 * np.pi * viscosity * r_max)
                    Gamma_max = D_max * (q * q) * (10 ** 9)
                    if (select_methods.get() == "Quad Fit"):
                        # Initial guess
                        color = '#00' + str(index) + '000'
                        label = 'Data Set ' + str(index)
                        # Get the sigma of constant and find sigma of radius
                        sigma = np.sqrt(Fit_Param[1])
                        gamma_dist, gamma_prob = Distrubution(Fit_Param[0], sigma, 0, 0, Gamma_min, Gamma_max)

                    elif (select_methods.get() == 'Cubic Fit'):
                        color = '#00'+ str(index) +'000'
                        label = 'Data Set ' + str(index)
                        # Get the sigma of constant and find sigma of radius
                        sigma = np.sqrt(Fit_Param[1])
                        skew = Fit_Param[2]**(1/3)
                        gamma_dist, gamma_prob = Distrubution(Fit_Param[0], sigma, skew, 0, Gamma_min, Gamma_max)

                    elif (select_methods.get() == 'Quartic Fit'):
                     # Initial guess
                        color = '#00'+ str(index) +'000'
                        label = 'Data Set ' + str(index)
                        # Get the sigma of constant and find sigma of radius
                        sigma = np.sqrt(Fit_Param[1])
                        skew = Fit_Param[2] ** (1 / 3)
                        kert = Fit_Param[3] **(1/4)
                        gamma_dist, gamma_prob = Distrubution(Fit_Param[0], sigma, skew, kert, Gamma_min,Gamma_max)
                D_dist = np.divide(gamma_dist, (q ** 2)) * (10 ** -9)
                r_dist = np.divide(boltzman_con * temperature, np.multiply(6 * np.pi * viscosity, D_dist))
                # If there is no value in min or max size to autofill a default value in. min is .01 times smaller and max is 100 times larger
                Plot_Prob(r_dist, gamma_prob, color, label)
            elif(select_methods.get()=='Contin' or select_methods.get()=="Choose..."):
                Plot_Prob(r_dist,x,color,label)
                gamma_prob = x
            # Plot the Normal Data Based on the model

            if (normal == 0):
                Fit_data = np.add(np.multiply(Fit_data, B_holder), B_holder)
            Plot_G2(tau_holder, Fit_data, color, label)
            # Plot the normal error
            error = np.subtract(g2_holder, Fit_data)
            Plot_Error(tau_holder, np.divide(error, g2_holder), color, label)

            #Exporting the fit data
            filename  = File_Names_Long[index]
            tau = tau_holder
            g2 = g2_holder
            error = np.divide(error, g2_holder)
            beta = beta_holder
            B = B_holder
            MultiExport()

    #Load Data Button
    Button_MultiLoad =Button(AngleWindow,text = "Load Data Sets",command=LoadMulti)
    Button_MultiLoad.place(height = 30,width =120,x=10,y=0)

    #Button Run Multi
    Button_RunMulti = Button(AngleWindow,text = ' Run Multi', command = RunMulti)
    Button_RunMulti.place(height= 30, width = 80, x = 500, y=0)

    Entry_MultiAngle = []
    #Export Button


#Function to close the window for L-curve
def on_closing():
    #if messagebox.askokcancel("Quit", "Do you want to quit?"):
   root.quit()
root.protocol("WM_DELETE_WINDOW", on_closing)

#Set up the function for L-curve
def LCurve():
    # Function to close the window for L-curve
    #def on_closing_Top():
     #   if messagebox.askokcancel("Quit", "Do you want to quit?"):
      #      Top.quit()
    #Warn user of long computation time
    result = messagebox.askquestion("L-Curve", "The L-Curve will take time, Do you still want to proceed")
    if (result == 'yes'):
        global wavelength
        wavelength = int(Entry_wavelength.get()) * 10 ** (-9)
        # refractive
        global refractive
        refractive = float(Entry_refractive.get())
        # Viscosity
        global viscosity
        viscosity = float(Entry_viscosity.get())
        # Temperant
        global temperature
        temperature = float(Entry_temperant.get()) + 273
        # Scattering
        global scattering
        scattering = float(Entry_scattering.get())
        # Boltzman Constant
        global boltzman_con
        boltzman_con = 1.38064852 * 10 ** (-23)
        global q
        q = (4 * np.pi * refractive / wavelength) * np.sin(np.deg2rad(scattering / 2))
        global g2, beta
        tau_high = float(Entry_high.get())
        tau_low = float(Entry_low.get())
        beta = float(Entry_beta.get())
        # Get baseline from low tau and high tau
        count = len(tau) - 1
        B = 0
        index = 0
        # Go until lower limit is reached and record only the given data
        while (tau[count] > tau_low):
            if (tau[count] < tau_high):
                B = g2[count] + B
                index = index + 1
            count = count - 1
        B = B / index

        #Get an average from the Cumulant expansions incase the user did not already specify a size range
        g2Norm, normal, beta,g1_squared = Normalize_Check(tau, g2, B, beta)
        g2Norm = np.multiply(g2Norm,beta)
        Params1, Fit = DLS_linear_fit.DLS_linear_fit(tau, g2, B, beta)
        initialGuess = [Params1[0], 1]
        OldParams, OldData = DLS_quad_fit.DLS_quad_fit(tau, g2, B, beta, initialGuess)
        initialGuess = [OldParams[0], OldParams[1], 1]
        OldParams, OldData = DLS_cubic_fit.DLS_cubic_fit(tau, g2, B, beta, initialGuess)
        initialGuess = [OldParams[0], OldParams[1], OldParams[2], 1]
        Fit_Param, Fit_data = DLS_quartic_fit.DLS_quartic_fit(tau, g2Norm, B, beta, initialGuess)

        # Get the diffusive coeff and average radius
        Diffusive_Coef = Fit_Param[0] / (np.square(q))
        r_mean = (boltzman_con * temperature) / (6 * np.pi * viscosity * Diffusive_Coef) * 10 ** 9

        # If there is no value in min or max size to autofill a default value in. min is .01 times smaller and max is 100 times larger
        # Get the sigma of constant and find sigma of radius
        if (select_dist_type.get() == "Linear"):
            PlotRange = 10
        else:
            PlotRange = 100
        if (Entry_min.get() == '' or Entry_min.get() == 'N/A'):
            Entry_min.delete(0, END)
            r_min = r_mean / PlotRange
            Entry_min.insert(0, str(float('%.3g' % r_min)))
        if (Entry_max.get() == '' or Entry_max.get() == 'N/A'):
            Entry_max.delete(0, END)
            r_max = r_mean * PlotRange
            Entry_max.insert(0, str(float('%.3g' % r_max)))

        #Initialize the size range for the regularization methods
        g1_squared = g2Norm
        n_pts = int(Entry_points.get())
        r_min = float(Entry_min.get())
        r_max = float(Entry_max.get())

        #Get the initial guess for the first calculation in the L-curve
        if (select_dist_type.get() == 'Logarithmic'):
            r_dist = np.logspace(np.log10(r_min), np.log10(r_max), n_pts)
        else:
            r_dist = np.linspace(r_min, r_max, n_pts)
        gamma = np.multiply(np.divide(1, r_dist), (16 * boltzman_con * temperature * np.pi * (refractive ** 2) * (
                    (np.sin(np.deg2rad(scattering / 2))) ** 2) / (6 * viscosity * wavelength * wavelength)) * 10 ** 9)
        W = np.diag(np.ones(len(g1_squared)))

        # set up alpha distrtibution and place holders for dist norm and residual norm
        Npoints = 20
        alpha_dist = np.logspace(np.log10(10 ** (-6)), np.log10((10 ** (2))), Npoints)
        x_square = np.ones(Npoints)
        res_square = np.ones(Npoints)

        #Set up the initial guesses for the algorthuims
        #Initial Guess for the CONTIN algorithuim
        x0_in = np.ones(len(gamma)) + (10 ** (-3)) * np.random.rand(1, len(gamma))
        x0 = np.divide(x0_in, np.sum(x0_in))

        #Initial Guess for the REPES algorthuim
        if (select_methods.get() == "REPES"):
            Tau = np.divide(1, gamma)
            a0 = np.ones(len(Tau))
            a0 = np.divide(a0, np.sum(a0))
            g1_squared = [0 if i < 0 else i for i in g1_squared]
            g1 = np.sqrt(g1_squared)

        #For loop to run through the alpha parameters and get the values of the dist and norm sizes
        for index in range(len(x_square)):
            #L-Curve code if REPES method is selected
            if (select_methods.get() == "REPES"):
                Fit_data, gamma, x = DLS_REPES.DLS_REPES(tau, g1, Tau, alpha_dist[index], a0)
                MethodName = " for REPES"
                Fit_data = np.square(Fit_data)
            #L-Curve code if CONTIN is selected
            else:
                Fit_data, gamma, x = DLS_contin.DLS_contin(tau, g1_squared, gamma, alpha_dist[index], W, x0)
                MethodName = " for CONTIN"
           #Ensure that the solution for the previos problem becomes the starting point for the next itteration. x is always the distribution solution
            a0 = x
            x0 = x
            x = x / np.sum(x)

            Fit_data = np.multiply(Fit_data, beta)
            Error = np.subtract(Fit_data, g2Norm)
            Error = np.square(Error)
            res_square[index] = np.sqrt(np.sum(Error))
            x_square[index] = np.sqrt((np.sum(np.square(x))))

        # Solving for two different Optimums. The CONTIN algorthuim peferes wide distributions and as such that is the trade-off with increasing alpha (i.e. Wide Dist. Smalller Error or Smaller Error and wider Dist).
        # The REPES algorithum wants to eliminate error all together and as such requires the opposite trend (i.e. when to stop getting rid of the noise of the data)
        if (select_methods.get() == "REPES"):
            x_int = min(res_square)
            y_int = min(x_square)
        # If the algorthuim is CONTIN
        else:
            x_int = min(res_square)
            y_int = min(x_square)


        # Find the Interval of Interest. Must Ensure that it is within the point range that has a differing point (i.e not the extremes)
        #Need to rescale the erors using logs and the [-10,10] scale
        #Change the scale to log scale to a [-10 , 10] scale
        x_square_log       = np.log10(np.copy(x_square))
        res_square_log     = np.log10(np.copy(res_square))
        x_square_log_new   = np.copy(x_square)
        res_square_log_new = np.copy(res_square)
        for index in range(len(x_square)):
            res_square_log_new[index] = ((10+10)/(max(res_square_log)-min(res_square_log)))*(res_square_log[index] - (10*min(res_square_log) - (-10)*max(res_square_log))/((10+10)))
            x_square_log_new[index]   = ((10+10)/(max(x_square_log)-min(x_square_log)))*(x_square_log[index] - (10*min(x_square_log) - (-10)*max(x_square_log))/((10+10)))
        lowest_dist = 10**8
        #Set up the intercept if its REPES
        if(select_methods.get() == "REPES"):
            x_int_log_new   = 0#min(x_square_log_new)
            res_int_log_new = 0#max(res_square_log_new)
        #Set yp the interscept if its CONTIN
        else:
            x_int_log_new   = 0#min(x_square_log_new)
            res_int_log_new = 0#min(res_square_log_new)

        #Find the closest point to the intercept
        for index in range(len(x_square)):
            dist = np.sqrt(np.square((x_square_log_new[index] - x_int_log_new) / np.max(x_square_log_new)) + np.square(
                (res_square_log_new[index] - res_int_log_new) / np.max(res_square_log_new)))
            if (dist < lowest_dist):
                lowest_dist = dist
                index_opt = index
        #Create a new window that will pop up with L-curve information
        Top = Toplevel()
        #Top.protocol("WM_DELETE_WINDOW", on_closing_Top)
        Top.title("DLS GUI L-Curve" + MethodName)
        Top.geometry("1000x450")#(x,y)
        # Display optimal value
        Label_opt = Label(Top, text='Optimal Value is: ' + str(float('%.2g' % alpha_dist[index_opt])))
        Label_opt.place(x=400,y=0,width=210,height=45)#place(height=20,width=200, x=200, y=600)

        #Plot the L-curve
        Figure_L, ax_L = plt.subplots(figsize=(5, 4))
        ax_L.ticklabel_format(style='sci', axis='y',scilimits = (0,0))
        ax_L.ticklabel_format(style='sci', axis='x',scilimits = (0,0))
        ax_L.plot(res_square, x_square, 'bo', fillstyle='none')
        ax_L.plot(res_square[index_opt], x_square[index_opt], 'r*')
        ax_L.axhline(y_int)
        ax_L.axvline(x_int)
        ax_L.plot(x_int, y_int, 'k*')
        ax_L.set_title("L-Curve" + MethodName)
        ax_L.set_xlabel("Residual Norm")
        ax_L.set_ylabel("Dist. Norm")
        Figure_L.tight_layout()
        canvas_L = FigureCanvasTkAgg(Figure_L, Top)
        canvas_L.get_tk_widget().place(x=0,y=46,width=500,height=400)#place(height=600, width=500, x=10, y=0)

        # Plot the residuals/dist vs  alpha
        Figure_Res, ax1 = plt.subplots(figsize=(5, 4))
        #Plot the Residual Norm
        ax1.ticklabel_format(style='sci', axis='y',scilimits = (0,0))
        ax1.set_xscale('log')
        ax1.plot(alpha_dist, x_square, 'bo', fillstyle='none')
        ax1.plot(alpha_dist[index_opt], x_square[index_opt], 'b*')
        ax1.set_title('Residuals')
        ax1.set_xlabel('Reguralization Parameter')
        ax1.set_ylabel('Dist. Residual', color='blue')
        ax2 = ax1.twinx()

        #Plot the distribution norm
        ax2.ticklabel_format(style='sci', axis='y',scilimits = (0,0))
        ax2.plot(alpha_dist, res_square, 'ro', fillstyle='none')
        ax2.plot(alpha_dist[index_opt], res_square[index_opt], 'r*')
        ax2.set_ylabel('Norm Residual', color='red')
        Figure_Res.tight_layout()
        canvas_Res = FigureCanvasTkAgg(Figure_Res, Top)
        canvas_Res.get_tk_widget().place(x=500,y=46,width=500,height = 400)

        # Export the LCurve data to a txt file
        # Set up to export the L-Curve values
        Export_L_Curve(alpha_dist, res_square, x_square)

#When asked to run the system will save all the input data as a global variable

#When asked to run the system will save all the input data as a global variable
def Run():
    #Bring in the global variable
    global g2, tau, tau_high , tau_low,beta
    # Globalize the data for extraction
    global error
    global r_dist
    global x,gamma_prob
    global Fit_data, Fit_Param
    global B

    tau_high = float(Entry_high.get())
    tau_low = float(Entry_low.get())
    beta = float(Entry_beta.get())
    #Get baseline from low tau and high tau
    count = len(tau)-1
    B = 0
    index =0
    #Go until lower limit is reached and record only the given data
    while(tau[count]>tau_low):
        if(tau[count]<tau_high):
            B = g2[count]+B
            index = index +1
        count =count-1
    B = B/index

    #wavelength
    global wavelength
    wavelength = int(Entry_wavelength.get()) * 10**(-9)
    #refractive
    global refractive
    refractive = float(Entry_refractive.get())

    #Viscosity
    global viscosity
    viscosity = float(Entry_viscosity.get())

    #Temperant
    global temperature
    temperature = float(Entry_temperant.get()) +273

    #Scattering
    global scattering
    scattering = float(Entry_scattering.get())

    #Boltzman Constant
    global boltzman_con
    boltzman_con = 1.38064852*10**(-23)

    #Get the scattering vectore q
    global q
    q = (4*np.pi*refractive/wavelength)*np.sin(np.deg2rad(scattering/2))

    #Check if the data is linearized
    g2Norm,normal,beta,g1_squared = Normalize_Check(tau,g2,B,beta)


    #Use the selected data method
    if(select_methods.get()=="Linear Fit"):
        Fit_data,Fit_Param = DLS_linear_fit.DLS_linear_fit(tau,g2Norm,B,beta)
        #Input Information to the GUI based on Model
        Entry_min.delete(0,'end')
        Entry_min.insert(0,'N/A')
        Entry_max.delete(0,'end')
        Entry_max.insert(0,'N/A')
    #Get the Size range (Min and Max)
    else:
        if(select_methods.get()=="Quad Fit"or select_methods.get()=="DYNALS" or select_methods.get()=="Contin" or select_methods.get()=="NNLS" or select_methods.get()=="REPES"):
            Params1, Fit = DLS_linear_fit.DLS_linear_fit(tau, g2, B, beta)
            initialGuess = [Params1[0], 1]
            Fit_Param, Fit_data = DLS_quad_fit.DLS_quad_fit(tau, g2Norm, B, beta,initialGuess)
        elif(select_methods.get()=='Cubic Fit'):
            # Initial guess
            Params1, Fit = DLS_linear_fit.DLS_linear_fit(tau, g2, B, beta)
            initialGuess = [Params1[0], 1]
            OldParams, OldData = DLS_quad_fit.DLS_quad_fit(tau, g2, B, beta,initialGuess)
            initialGuess = [OldParams[0], OldParams[1], 1]
            Fit_Param,Fit_data = DLS_cubic_fit.DLS_cubic_fit(tau,g2Norm,B,beta,initialGuess)
        elif(select_methods.get() == 'Quartic Fit'):
            # Initial guess
            Params1, Fit = DLS_linear_fit.DLS_linear_fit(tau, g2, B, beta)
            initialGuess = [Params1[0], 1]
            OldParams, OldData = DLS_quad_fit.DLS_quad_fit(tau, g2, B, beta, initialGuess)
            initialGuess = [OldParams[0], OldParams[1], 1]
            OldParams, OldData = DLS_cubic_fit.DLS_cubic_fit(tau, g2, B, beta,initialGuess)
            initialGuess = [OldParams[0], OldParams[1], OldParams[2], 1]
            Fit_Param , Fit_data = DLS_quartic_fit.DLS_quartic_fit(tau,g2Norm,B,beta,initialGuess)
        #Get the diffusive coeff and average radius
        Diffusive_Coef = Fit_Param[0]/(np.square(q))
        r_mean = (boltzman_con*temperature)/(6*np.pi*viscosity*Diffusive_Coef)*10**9

       #Get the min and max range for the R-dist. Linear dist needs a smaller range because of the way the plot looks
        if(select_dist_type.get()=="Linear"):
            PlotRange = 10
        else:
            PlotRange = 100
        #If there is no value in min or max size to autofill a default value in. min is .01 times smaller and max is 100 times larger
        if(Entry_min.get()== '' or Entry_min.get()=='N/A'):
            #Consider if the dist type is linear or logarthimic
            Entry_min.delete(0,END)
            r_min = r_mean/PlotRange
            Entry_min.insert(0, str(float('%.3g' % r_min)))
        if(Entry_max.get()== ''or Entry_max.get()=='N/A'):
            Entry_max.delete(0,END)
            r_max = r_mean*PlotRange
            Entry_max.insert(0, str(float('%.3g' % r_max)))
        if(Entry_regularization.get()==''):
            Entry_regularization.insert(0,str('0.5'))

    #Get the distributions for the Method Selected
    if(select_methods.get()=="Linear Fit"):
        color = 'b'
        label = 'linear'
        #Input Information to the GUI based on Model
        Entry_min.delete(0,'end')
        Entry_min.insert(0,'N/A')
        Entry_max.delete(0,'end')
        Entry_max.insert(0,'N/A')
    #Use the distribution function to get the Cumulant Size Distributions
    elif(select_methods.get()=="Quad Fit" or select_methods.get()=="Cubic Fit" or select_methods.get()=='Quartic Fit'):
        r_min = float(Entry_min.get())
        r_max = float(Entry_max.get())
        # Get upper and lower bounds
        D_min = boltzman_con * temperature / (6 * np.pi * viscosity * r_min)
        Gamma_min = D_min * (q * q) * (10 ** 9)
        D_max = boltzman_con * temperature / (6 * np.pi * viscosity * r_max)
        Gamma_max = D_max * (q * q) * (10 ** 9)
        if(select_methods.get()=="Quad Fit"):
            color = 'g'
            label = 'Quad'
            # Get the sigma of constant and find sigma of radius
            sigma = np.sqrt(Fit_Param[1])
            gamma_dist,gamma_prob = Distrubution(Fit_Param[0],sigma,0,0, Gamma_min,Gamma_max)

        elif(select_methods.get()=='Cubic Fit'):
            color = 'c'
            label = 'Cubic'
            # Get the sigma of constant and find sigma of radius
            sigma = np.sqrt(Fit_Param[1])
            gamma_dist, gamma_prob = Distrubution(Fit_Param[0], sigma, Fit_Param[2], 0, Gamma_min,Gamma_max)

        elif(select_methods.get() == 'Quartic Fit'):
            color = 'm'
            label = 'Quartic'
            # Get the diffusive coeff and average radius
            Diffusive_Coef = Fit_Param[0] / (np.square(q))
            # Get the sigma of constant and find sigma of radius
            sigma = np.sqrt(Fit_Param[1])

            gamma_dist, gamma_prob = Distrubution(Fit_Param[0], sigma, Fit_Param[2], Fit_Param[3], Gamma_min,Gamma_max)
        D_dist = np.divide(gamma_dist,(q**2))*(10**-9)
        r_dist = np.divide(boltzman_con*temperature,np.multiply(6*np.pi*viscosity,D_dist))
        #If there is no value in min or max size to autofill a default value in. min is .01 times smaller and max is 100 times larger
        Plot_Prob(r_dist, gamma_prob, color,label)
    #If the Algorithuim desired is CONTIN or NNLS (alpha = 0)
    elif(select_methods.get()=="Contin" or select_methods.get()=="NNLS"):
        #g1_squared = np.divide(g2Norm,beta)
        reg_param = float(Entry_regularization.get())
        n_pts = int(Entry_points.get())
        r_min = float(Entry_min.get())
        r_max = float(Entry_max.get())
        if(select_dist_type.get() == 'Logarithmic'):
            r_dist = np.logspace(np.log10(r_min),np.log10(r_max),n_pts)
            Plot_prob.set_xscale('log')
        else:
            r_dist = np.linspace(r_min,r_max,n_pts)
        gamma = np.multiply(np.divide(1,r_dist),(16*boltzman_con*temperature*np.pi*(refractive**2)*((np.sin(np.deg2rad(scattering/2)))**2)/(6*viscosity*wavelength*wavelength))*(10**9))
        W =np.diag(np.ones(len(g1_squared)))
        #Initial Conditions
        x0_in = np.ones(len(gamma)) + (10 ** (-3)) * np.random.rand(1, len(gamma))
        x0 = np.divide(x0_in, np.sum(x0_in))
        #If NNLS
        if(select_methods.get()=="NNLS"):
            reg_param = 0
            color = 'k'
            label = 'NNLS'
            Fit_data, gamma, x = DLS_NNLS.DLS_NNLS(tau, g1_squared, gamma, reg_param, W, x0)
        else:
            color = 'y'
            label = 'Contin'
            Fit_data, gamma, x = DLS_contin.DLS_contin(tau, g1_squared, gamma, reg_param, W, x0)

        x = x/np.sum(x)
        Fit_data = np.multiply(Fit_data,beta)
        Plot_Prob(r_dist, x, color,label)
    #If REPES is sleceted
    elif(select_methods.get()=="REPES"):
        #g1_squared = np.divide(g2Norm,beta)
        g1_squared = [0 if i<0 else i for i in g1_squared]
        g1 = g1_squared #np.sqrt(g1_squared)
        reg_param = float(Entry_regularization.get())
        n_pts = int(Entry_points.get())
        r_min = float(Entry_min.get())
        r_max = float(Entry_max.get())
        if(select_dist_type.get() == 'Logarithmic'):
            r_dist = np.logspace(np.log10(r_min),np.log10(r_max),n_pts)
            Plot_prob.set_xscale('log')
        else:
            r_dist = np.linspace(r_min,r_max,n_pts)
        gamma = np.multiply(np.divide(1,r_dist),(16*boltzman_con*temperature*np.pi*(refractive**2)*((np.sin(np.deg2rad(scattering/2)))**2)/(6*viscosity*wavelength*wavelength))*10**9)
        Tau = np.divide(1,gamma)
        a0 = np.ones(len(Tau))
        a0 = np.divide(a0, np.sum(a0))
        Fit_data,gamma,x = DLS_REPES.DLS_REPES(tau,g1,Tau,reg_param,a0)
        Fit_data = np.multiply((Fit_data),beta)
        color = 'orangered'
        label = 'REPES'
        Plot_Prob(r_dist, x, color,label)
    #If DYNALS is sleceted
    elif(select_methods.get()=="DYNALS"):
        #g1_squared = np.divide(g2Norm,beta)
        g1_squared = [0 if i<0 else i for i in g1_squared]
        g1 = np.sqrt(g1_squared)
        n_pts = int(Entry_points.get())
        r_min = float(Entry_min.get())
        r_max = float(Entry_max.get())
        if(select_dist_type.get() == 'Logarithmic'):
            r_dist = np.logspace(np.log10(r_min),np.log10(r_max),n_pts)
            Plot_prob.set_xscale('log')
            dist = 0
        else:
            r_dist = np.linspace(r_min,r_max,n_pts)
            dist = 1
        gamma = np.multiply(np.divide(1,r_dist),(16*boltzman_con*temperature*np.pi*(refractive**2)*((np.sin(np.deg2rad(scattering/2)))**2)/(6*viscosity*wavelength*wavelength))*10**9)
        x0_in = np.ones(len(gamma)) + (10 ** (-3)) * np.random.rand(1, len(gamma))
        x0 = np.divide(x0_in, np.sum(x0_in))
        Fit_data,gamma,x = DLS_DYNALS.DLS_DYNALS (tau,g1,gamma,x0)
        Fit_data = np.multiply(np.square(Fit_data),beta)
        color = 'deepskyblue'
        label = 'DYNALS'
        Plot_Prob(r_dist, x, color,label)

    #Plot the Normal Data Based on the model.
    #Unnormalize the data if it is normalized
    if(normal==0):
        Fit_data = np.add(np.multiply(Fit_data,B),B)
    #Plot the fit onto the data
    Plot_G2(tau,Fit_data,color,label)

   #Plot the normal error
    error = np.subtract(g2,Fit_data)
    error = np.divide(error, g2)
    Plot_Error(tau,error,color,label)

#This function fills in the Min and Max size Entries if a regularization method is selected from the drop down menu
def Drop_Fun(Input):
    if (select_methods.get() == "DYNALS" or select_methods.get() == "Contin" or select_methods.get() == "NNLS" or select_methods.get() == "REPES"):

        # Bring in the global variable
        global g2, tau, tau_high, tau_low, beta

        tau_high = float(Entry_high.get())
        tau_low = float(Entry_low.get())
        beta = float(Entry_beta.get())
        # Get baseline from low tau and high tau
        count = len(tau) - 1
        B = 0
        index = 0
        # Go until lower limit is reached and record only the given data
        while (tau[count] > tau_low):
            if (tau[count] < tau_high):
                B = g2[count] + B
                index = index + 1
            count = count - 1
        B = B / index

        # wavelength
        global wavelength
        wavelength = int(Entry_wavelength.get()) * 10 ** (-9)

        # refractive
        global refractive
        refractive = float(Entry_refractive.get())

        # Viscosity
        global viscosity
        viscosity = float(Entry_viscosity.get())

        # Temperant
        global temperature
        temperature = float(Entry_temperant.get()) + 273

        # Scattering
        global scattering
        scattering = float(Entry_scattering.get())

        # Boltzman Constant
        global boltzman_con
        boltzman_con = 1.38064852 * 10 ** (-23)

        global q
        q = (4 * np.pi * refractive / wavelength) * np.sin(np.deg2rad(scattering / 2))

        # Check if the data is linearized
        g2Norm, normal, beta,g1_squared = Normalize_Check(tau, g2, B, beta)
        # Initial guess
        Params1, Fit = DLS_linear_fit.DLS_linear_fit(tau, g2, B, beta)
        initialGuess = [Params1[0], 1]
        OldParams, OldData = DLS_quad_fit.DLS_quad_fit(tau, g2, B, beta, initialGuess)
        initialGuess = [OldParams[0], OldParams[1], 1]
        OldParams, OldData = DLS_cubic_fit.DLS_cubic_fit(tau, g2, B, beta,initialGuess)
        initialGuess = [OldParams[0], OldParams[1], OldParams[2], 1]
        Fit_Param, Fit_data = DLS_quartic_fit.DLS_quartic_fit(tau, g2Norm, B, beta,initialGuess)
         # Get the diffusive coeff and average radius
        Diffusive_Coef = Fit_Param[0] / (np.square(q))
        r_mean = (boltzman_con * temperature) / (6 * np.pi * viscosity * Diffusive_Coef) * 10 ** 9

        # Get the sigma of constant and find sigma of radius
        # If there is no value in min or max size to autofill a default value in. min is .01 times smaller and max is 100 times larger
        if (Entry_min.get() == '' or Entry_min.get() == 'N/A'):
            Entry_min.delete(0, END)
            r_min = r_mean / 100
            Entry_min.insert(0, str(float('%.3g' % r_min)))
        if (Entry_max.get() == '' or Entry_max.get() == 'N/A'):
            Entry_max.delete(0, END)
            r_max = r_mean * 100
            Entry_max.insert(0, str(float('%.3g' % r_max)))
        if (Entry_regularization.get() == ''):
            Entry_regularization.insert(0, str('0.5'))

#Clear only the graphs and eventually allow for the clearing of select graphs
def Clear():
    #Clear the Min and Max Peaks etc.
    Entry_min.delete(0,'end')
    Entry_max.delete(0,'end')

    #Clear the error, prob dist, g2 graph
    len_err = len(Plot_error.lines)
    len_prob = len(Plot_prob.lines)
    len_g2 = len(Plot_g2.lines)
    #Clear the error
    for i in range(len_err-1):
        del Plot_error.lines[1]
        Plot_error.legend().remove()
        canvas_Error = FigureCanvasTkAgg(Figure_Error, root)
        canvas_Error.get_tk_widget().place(height=400, width=500, x=510, y=455)
    Plot_error.legend().remove()

    #Clear the prob dist
    for i in range(len_prob-1):
        del Plot_prob.lines[1]
        canvas_prob = FigureCanvasTkAgg(Figure_prob, root)
        canvas_prob.get_tk_widget().place(height=415, width=510, x=0, y=285)
    Plot_prob.legend().remove()
    Plot_prob.relim()

    #Clear the G2 graph
    for i in range(len_g2-4):
        del Plot_g2.lines[4]
        canvas_g2 = FigureCanvasTkAgg(Figure_g2, root)
        canvas_g2.get_tk_widget().place(height=400, width=500, x=510, y=50)
    #Removes any writing in the legends but if there was a legend it will leave a blank space that is not very noticable
    Plot_g2.legend().remove()
    Plot_g2.legend()

#For the rounding up the linear plot area
def round_up(n, decimals=0):
    multiplier = 10 ** decimals
    return math.ceil(n * multiplier) / multiplier

#For roundind down down the linear plot area
def round_down(n, decimals=0):
    multiplier = 10 ** decimals
    return math.floor(n * multiplier) / multiplier

#Plot the R-dsit
def Plot_Prob(xdata,ydata,color,labels):

    if(select_dist_type.get() == "Logarithmic"):
        dist = 0
        Plot_prob.set_xscale('log')
    else:
        dist = 1
        Plot_prob.set_xscale('linear')

    Plot_prob.plot(xdata, ydata, color =color,marker='o',label = labels,markersize = 3)
    Plot_prob.legend(loc = 'upper right',frameon=True)
    #Find the range of points that satisfy criteria
    #FInd index of max and min range
    #Figure out how to adapt to log plot
    minX =0
    maxX =0
    tol = 10**-4
    for index in range(len(ydata)):
        if(tol<ydata[index] and minX==0):
            minX = xdata[index]
            maxX = xdata[index]
        elif(tol < ydata[index] and minX > 0):
            maxX = xdata[index]
    #Get the rounding based on dist type
    #old___ will be used to find the place holder (i.e tens ones hundreds etc)
    oldmax = maxX
    oldmin  = minX
    #Get the lowest Exponent Value
    c = 1
    pos_nums = []
    while minX !=0 :
            z = minX % 10
            pos_nums.append(z*c)
            minX = minX //10
            c = c*10
    minX = 10**(len(pos_nums)-1)
    #Check if Max is too close to a power of ten add decrease by a power of ten

    if(oldmin/minX < 2.5):
        minX = 10**(len(pos_nums)-2)

    #Get the highest exponent value
    c = 1
    pos_nums = []
    while maxX != 0 :
            z = maxX % 10
            pos_nums.append(z*c)
            maxX = maxX // 10
            c = c*10
    maxX = 10**(len(pos_nums))

    #Check if Max is too close to a power of ten add another power of ten
    if(oldmax/maxX > 0.75):
        maxX = 10**(len(pos_nums)+1)

    #Check if the power falls above the set point ranges
    if(minX < min(xdata)):
        minX = min(xdata)
    if(maxX > max(xdata)):
        maxX = max(xdata)
    #If the distribution is log type
    if(dist==0):
        x_labels = np.logspace(np.log10(minX),np.log10(maxX),5)
    #If the dist type is linear
    else:
        #Give the graph a 20 percent over and under shoot of the non-zero components
        #Get the minimum and maximum for the linear dist
        # Check if the power falls above the set point ranges
        oldmin = round_down(oldmin,-len(str(int(oldmin))))
        oldmax = round_up(oldmax,1-len(str(int(oldmax))))
        #oldmin = oldmin - 0.50 * oldmin
        #oldmax = oldmax+0.5*oldmax
        if (oldmin-0.50*oldmin < min(xdata)):
            oldmin = min(xdata)
        if (oldmax+0.50*oldmax > max(xdata)):
            oldmax = max(xdata)
        x_labels = np.linspace(oldmin , oldmax,5)

    x_labels[0] = np.round(x_labels[0])
    x_labels[1] = np.round(x_labels[1])
    x_labels[2] = np.round(x_labels[2])
    x_labels[3] = np.round(x_labels[3])
    x_labels[4] = np.round(x_labels[4])



    #Insert Code to Remove previos axis (okay if none before the algorthuim will push forward
    Plot_prob.xaxis.set_minor_formatter(NullFormatter())
    Plot_prob.yaxis.set_minor_formatter(NullFormatter())
    Plot_prob.set_xlim([x_labels[0],x_labels[4]])
    #Prevent Minor Tick Labeling From Interfering with Plots
    Plot_prob.set_xticklabels(x_labels)
    Plot_prob.set_xticks(x_labels)
    canvas_prob = FigureCanvasTkAgg(Figure_prob, root)
    canvas_prob.get_tk_widget().place(height=415, width=510, x=0, y=285)

#Plot the data or data firs
def Plot_G2(xdata,ydata,color,labels):
    global g2

    Plot_g2.plot(xdata, ydata, color,label = labels,fillstyle = 'none')
    Plot_g2.legend(loc = 'lower left',frameon=True)
    Plot_g2.autoscale()

    #Set up the rescaling of the y-axis
    minX  = min(tau)*0.95
    maxX  = max(tau)*1.05
    Plot_g2.set_xlim([minX , maxX])

    #Set up the rescaling of the y-axis
    if(len(g2) != 0):
        minY  = min(g2)*0.996
        maxY  = max(g2)*1.004
        Plot_g2.set_ylim([minY , maxY])

    canvas_g2 = FigureCanvasTkAgg(Figure_g2, root)
    canvas_g2.get_tk_widget().place(height=400, width=500, x=510, y=50)

#Plot the error of the data fits. Limit code could be simplified but still works well
def Plot_Error(xdata,ydata,color,labels):
    miny = round((min(ydata)*1.25),0 - int(np.ceil(np.log10(abs((min(ydata)*1.25)))-1)))
    maxy = round((max(ydata)*1.25),0 - int(np.ceil(np.log10(abs((max(ydata)*1.25)))-1)))
    Plot_error.plot(xdata, ydata,color,label = labels,linestyle = 'none',marker = 'o' ,fillstyle = 'none')
    Plot_error.legend(loc = 'upper right',frameon=True)
    Plot_error.set_ylim([1.25*miny,1.25*maxy])
    Plot_error.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    Plot_error.tick_params(axis='y', labelsize=9)
    canvas_Error = FigureCanvasTkAgg(Figure_Error, root)
    canvas_Error.get_tk_widget().place(height=400, width=500, x=510, y=455)

#Button to update the beta and the  baseline range lines when the user changes the values
def Update():
    global tau_low,tau_high,beta,count_update
    count_update= count_update+1
    if(count_update>1):
        del Plot_g2.lines[3] ,Plot_g2.lines[2] ,Plot_g2.lines[1]
    tau_low = float(Entry_low.get())
    tau_high= float(Entry_high.get())
    beta = float(Entry_beta.get())
    #tau low
    Plot_g2.lines.insert(1,Plot_g2.axvline(x=tau_low, color='b', linestyle='--'))
    del Plot_g2.lines[len(Plot_g2.lines)-1]
    canvas_g2 = FigureCanvasTkAgg(Figure_g2, root)
    canvas_g2.get_tk_widget().place(height=400, width=500, x=510, y=50)

    #tau high
    Plot_g2.lines.insert(1,Plot_g2.axvline(x=tau_high, color='r', linestyle='--'))
    del Plot_g2.lines[len(Plot_g2.lines)-1]
    canvas_g2 = FigureCanvasTkAgg(Figure_g2, root)
    canvas_g2.get_tk_widget().place(height=400, width=500, x=510, y=50)

    #beta
    Plot_g2.lines.insert(1,Plot_g2.axhline(y=beta, color='g', linestyle='--'))
    del Plot_g2.lines[len(Plot_g2.lines)-1]
    canvas_g2 = FigureCanvasTkAgg(Figure_g2, root)
    canvas_g2.get_tk_widget().place(height=400, width=500, x=510, y=50)

#Export rhe Datqa Generated From the L-curve. This is done automatically
def Export_L_Curve(RegPar,Resid,DistNorm):
    #Set up to export the L-Curve values
    global filename
    path = os.path.dirname(filename)
    file = filename
    #File Name
    file = file.replace(".txt","_L_Curve.txt")

    #Print the values gathered from the L-curve algorithuim
    with open(os.path.join(path, file), 'w') as fp:
        fp.write("|Reg. Param |Residual |Dist. Norm|\n\n")
        for index in range(len(RegPar)):
            # Account for when tau is larger than R
            fp.write("|" + str(float('%.4g' % RegPar[index])) + " " * (
                10 - len(str(float('%.4g' % RegPar[index])))) + "|" + str(
            float('%.4g' % Resid[index])) + " " * (
                10 - len(str(float('%.4g' % Resid[index])))) + "|" + str(
            float('%.4g' % DistNorm[index])) + " " * (10 - len(str(float('%.4g' % DistNorm[index])))) + "|\n")
    fp.close

    return

#Export the fit data,error,R-dist, Prob and Tau from the multiangle algorithuim
def MultiExport():
    #Identify the path to download it, need to figure out to figure out where to put the file
    global filename
    global tau, error , r_dist,x , Fit_data,gamma_prob,beta,B

    if (select_methods.get() == "DYNALS" or select_methods.get() == "Contin" or select_methods.get() == "NNLS" or select_methods.get() == "REPES"):
        x =x
    else:
        x = gamma_prob

    path = os.path.dirname(filename)
    file = filename
    #File Name
    file = file.replace(".txt","_FitData.txt")
    R_len = len(r_dist)
    Tau_len = len(tau)
    largest = Tau_len
    if(Tau_len<R_len):
        largest = R_len

    #Three cases equal length, tau larger and R-size larger
    with open(os.path.join(path, file), 'w') as fp:
        fp.write("|Method: "+ select_methods.get()+"|Alpha: "+ str(Entry_regularization.get())+"|"+"|Beta: "+ str("{:.5e}".format(B*beta+B))+"|"+"\n\n")
        fp.write("|Tau      |G2 Fit      |Error     |R-Size   |Prob \n\n")

        for index in range(largest):
            #Account for when tau is larger than R
            if(index > R_len-1):
                fp.write("|"+str(float('%.5g' % tau[index]))+" "*(9-len(str(float('%.5g' % tau[index]))))+"|"+str(float('%.5g' % Fit_data[index]))+" "*(9-len(str(float('%.5g' % Fit_data[index]))))+"|"+str(float('%.4g' % error[index])) + " " * (10 - len(str(float('%.4g' % error[index])))) + "|\n")

            elif(index > Tau_len-1):
                fp.write("|"+"         "+"|"+"         "+"|"+ "         "
                         + "|" + str(float('%.5g' % r_dist[index])) + " " * (9 - len(str(float('%.5g' % r_dist[index]))))
                         + "|" + str(float('%.5g' % x[index])) + " " * (12 - len(str(float('%.5g' % x[index]))))
                         + "|\n")
            else:
                fp.write("|"+str(float('%.5g' % tau[index]))+" "*(9-len(str(float('%.5g' % tau[index]))))+"|"+str(float('%.5g' % Fit_data[index]))+" "*(9-len(str(float('%.5g' % Fit_data[index]))))+"|"+str(float('%.4g' % error[index])) + " " * (10 - len(str(float('%.4g' % error[index]))))
                         + "|" + str(float('%.5g' % r_dist[index])) + " " * (9 - len(str(float('%.5g' % r_dist[index]))))
                         + "|" + str(float('%.5g' % x[index])) + " " * (12 - len(str(float('%.5g' % x[index]))))
                         + "|\n")
    fp.close
    return

#Export the fit data, error, R-dist, Prob AND TAU
def Export():
    #Identify the path to download it, need to figure out to figure out where to put the file
    global filename
    global tau, error , r_dist,x , Fit_data,gamma_prob,beta,B

    if (select_methods.get() == "DYNALS" or select_methods.get() == "Contin" or select_methods.get() == "NNLS" or select_methods.get() == "REPES"):
        x =x
    else:
        x = gamma_prob

    path = os.path.dirname(filename)
    file = filename
    #File Name
    file = file.replace(".txt","_FitData.txt")
    R_len = len(r_dist)
    Tau_len = len(tau)
    largest = Tau_len
    if(Tau_len<R_len):
        largest = R_len

    #Three cases equal length, tau larger and R-size larger
    with open(os.path.join(path, file), 'w') as fp:
        fp.write("|Method: "+ select_methods.get()+"|Alpha: "+ str(Entry_regularization.get())+"|"+"|Beta: "+ str("{:.5e}".format(B*beta+B))+"|"+"\n\n")
        fp.write("|Tau      |G2 Fit      |Error     |R-Size   |Prob \n\n")

        for index in range(largest):
            #Account for when tau is larger than R
            if(index > R_len-1):
                fp.write("|"+str(float('%.5g' % tau[index]))+" "*(9-len(str(float('%.5g' % tau[index]))))+"|"+str(float('%.5g' % Fit_data[index]))+" "*(9-len(str(float('%.5g' % Fit_data[index]))))+"|"+str(float('%.4g' % error[index])) + " " * (10 - len(str(float('%.4g' % error[index])))) + "|\n")

            elif(index > Tau_len-1):
                fp.write("|"+"         "+"|"+"         "+"|"+ "         "
                         + "|" + str(float('%.5g' % r_dist[index])) + " " * (9 - len(str(float('%.5g' % r_dist[index]))))
                         + "|" + str(float('%.5g' % x[index])) + " " * (12 - len(str(float('%.5g' % x[index]))))
                         + "|\n")
            else:
                fp.write("|"+str(float('%.5g' % tau[index]))+" "*(9-len(str(float('%.5g' % tau[index]))))+"|"+str(float('%.5g' % Fit_data[index]))+" "*(9-len(str(float('%.5g' % Fit_data[index]))))+"|"+str(float('%.4g' % error[index])) + " " * (10 - len(str(float('%.4g' % error[index]))))
                         + "|" + str(float('%.5g' % r_dist[index])) + " " * (9 - len(str(float('%.5g' % r_dist[index]))))
                         + "|" + str(float('%.5g' % x[index])) + " " * (12 - len(str(float('%.5g' % x[index]))))
                         + "|\n")
    fp.close
    return

#Export the parameters for the cumulant expansion, will spit out the parameters used for the generation of min and max size if a regularization method is selceted
def ExportParams():
    #Identify the path to download it, need to figure out to figure out where to put the file
    global filename
    global Fit_Param

    path = os.path.dirname(filename)
    file = filename
    #File Name
    file = file.replace(".txt","_FitParam.txt")
    length = len(Fit_Param)
    ParamPop = ''

    for index in range(4-length):
        Fit_Param = np.append(Fit_Param,0)


    #Three cases equal length, tau larger and R-size larger
    with open(os.path.join(path, file), 'w') as fp:
        fp.write('Fit Parameters for Gamma Distribution \n')
        fp.write("|Mean        |Standard Dev|Skew        |Kurt    \n\n")
        fp.write("|" + str(float('%.3g' % Fit_Param[0])) + " " * (12 - len(str(float('%.3g' % Fit_Param[0]))))+"|" + str(float('%.3g' % np.sqrt(Fit_Param[1]))) + " " * (12 - len(str(float('%.3g' % np.sqrt(Fit_Param[1])))))+"|" + str(float('%.3g' % Fit_Param[2]**(1/3))) + " " * (12 - len(str(float('%.3g' % Fit_Param[2]**(1/3)))))+"|" + str(float('%.3g' % Fit_Param[3]**(1/3))) + " " * (12 - len(str(float('%.3g' % Fit_Param[3]**(1/3))))))

    fp.close

    #ParamPop = ParamPop + "Fit Parameters for Gamma Distribution \n"
    ParamPop = 'Fit Parameters for Gamma Distribution \n'
    ParamPop = ParamPop + " |Mean        |Stand. Dev|Skew        |Kurt"
    ParamPop = ParamPop + "\n |" + str(float('%.3g' % Fit_Param[0])) + " " * (12 - len(str(float('%.3g' % Fit_Param[0]))))+"|" + str(float('%.3g' % np.sqrt(Fit_Param[1]))) + " " * (15 - len(str(float('%.3g' % np.sqrt(Fit_Param[1])))))+"|" + str(float('%.3g' % Fit_Param[2]**(1/3))) + " " * (15 - len(str(float('%.3g' % Fit_Param[2]**(1/3)))))+"|" + str(float('%.3g' % Fit_Param[3]**(1/3))) + " " * (14 - len(str(float('%.3g' % Fit_Param[3]**(1/3)))))

    #Create a Window With Parameters
    # Create a new window that will pop up with L-curve information
    ParamTop = Toplevel()
    # Top.protocol("WM_DELETE_WINDOW", on_closing_Top)
    ParamTop.title("Fit Parameters")
    ParamTop.geometry("500x60")  # (x,y)
    # Display optimal value
    ParamPop = str(ParamPop)
    Label_opt = Label(ParamTop, text=ParamPop , justify = 'left')
    Label_opt.place(x=0, y=0, width=280, height=60)  # place(height=20,width=200, x=200, y=600)

    return

#Fnction to control the wait time of the function that clears outputs if thr inputs are changed
def  handle_wait(event):
    #cancel the old job
    if root._after_id is not None:
        root.after_cancel(root._after_id)
    #create a new job
    #Responds in 650 ms
    root._after_id = root.after(650,Clear)

#Update the G2 beta and tau lines
def handle_wait_Beta(event):
    #cancel the old job
    if root._after_id is not None:
        root.after_cancel(root._after_id)
    #create a new job
    #Respons in 1,000 ms
    root._after_id = root.after(1000,Update)

#Function to clear the inputs when the plot type is changed
def DistTypeChanges(Input):
    Clear()
    return
#Function to clear the inputs when the plot type is changed
def TimeTypeChanges(Input):
    global TimeChange
    TimeChange = True
    LoadData()
    TimeChange = False
    return


#have counter for update
global count_update
count_update=0

#Set Up for the Inputs of the system
#Set up for the load data button
button_data = Button(root,text="No File")
button_data.place(height=30,width=210, x=130, y=0)
button_file = Button(root,text="Load Data", command= LoadData)
button_file.place(height=30,width=100, x=10, y=0)

#Set the time scale
time_type = [
    's',
    'µs'
]
#Hold the value that the user selects
select_time_type = StringVar()
select_time_type.set('µs')

#Create the drop down menu and label
Drop_time_type = OptionMenu(root,select_time_type,*time_type,command= TimeTypeChanges)
Drop_time_type.config(width=20)
Drop_time_type.place(height=25,width=45, x=415, y=32)
Drop_time_type = Label(root, text = "Time:")
Drop_time_type.place(height=30,width=35 , x=375, y=30)


#Set up for the user to select the data type if it s g1 or g2
dTypes = [
    "g1\u00b2",
    "g2"
]
#Hold the value that the user selects
select_dtype = StringVar()
select_dtype.set("g2")

#Create the drop down menu and label
Drop_dtype = OptionMenu(root,select_dtype,*dTypes)
Drop_dtype.config(width=20)
Drop_dtype.place(height=25,width=50, x= 425 , y=3)
Dtype_label = Label(root,text = "Data Type:")
Dtype_label.place(height=25,width=65, x= 355, y=0)

#Set up the wavelength entry
Entry_wavelength = Entry(root)
Entry_wavelength.place(height=25,width=110, x=10, y=60)
Entry_wavelength.insert(0,'637')
label_wavelength = Label(root,text="Wavelength (nm)")
label_wavelength.place(height=25,width=110, x=10, y=35)

#Bind the entries to the live updating algorithuim
#This is just used for the intializing of the entries
root._after_id = None
Entry_wavelength.bind('<Key>',handle_wait)


#Set up for Refrative Input
Entry_refractive = Entry(root)
Entry_refractive.place(height=25,width=110, x=175, y=60)
Entry_refractive.insert(0,'1.333')
Entry_refractive.bind('<Key>',handle_wait)
label_refractive = Label(root,text="Refractive Index")
label_refractive.place(height=25,width=110, x=175, y=35)

#Set up for Viscosity Input
Entry_viscosity = Entry(root)
Entry_viscosity.place(height=25,width=110, x=340, y=110)
Entry_viscosity.insert(0,'0.001')
Entry_viscosity.bind('<Key>',handle_wait)
label_viscosity = Label(root,text="Viscosity (Pa*s)")
label_viscosity.place(height=25,width=110, x=340, y=85)

#Set up for scattering Input
Entry_scattering = Entry(root)
Entry_scattering.place(height=25,width=110, x=10, y=110)
Entry_scattering.insert(0,'90')
Entry_scattering.bind('<Key>',handle_wait)
label_scattering = Label(root,text="Scattering Angle ("+u"\N{DEGREE SIGN}"+')')
label_scattering.place(height=25,width=130, x=10, y=85)

#Set up for Temperate Input
Entry_temperant = Entry(root)
Entry_temperant.place(height=25,width=110, x=175, y=110)
Entry_temperant.insert(0,'25')
Entry_temperant.bind('<Key>',handle_wait)
label_temperant = Label(root,text = "Temperature ("+u"\N{DEGREE SIGN}"+'C)')
label_temperant.place(height=25,width=110, x=175, y=85)

#Button to export data fit data
Button_export = Button(root,text = 'Export Fit',command = Export)
Button_export.place(height=25, width=110, x = 350,  y =260)

#Button to export Params for Fits
Button_exportParams = Button(root,text = 'Export Fit Parameters', command = ExportParams)
Button_exportParams.place(height=25, width=140, x = 170,  y =260)

#Set up for choosing fitting method

#List of Fitting Methods to choose from
methods = [
    "Linear Fit",
    "Quad Fit",
    "Cubic Fit",
    "Quartic Fit",
    "NNLS",
    "Contin",
    "REPES",
    "DYNALS"
]
#Hold the value that the user selects
select_methods = StringVar()
select_methods.set("Choose...")

#Create the drop down menu and label
Drop_Methods = OptionMenu(root,select_methods,*methods,command=Drop_Fun)
Drop_Methods.config(width=20)
Drop_Methods.place(height=25,width=120, x=110, y=150)
label_methods = Label(root, text = "Choose Fitting\nAlgorithum      ")
label_methods.place(height=40,width=110, x=1, y=140)
label_methods_colon = Label(root, text = ":")
label_methods_colon.place(height=40,width=10, x=100, y=140)

#Create the Run Button
button_enter = Button(root, text="Run",padx=20,command= Run)
button_enter.place(height=25,width=60, x=340, y=150)

#Create the L-Curve button
button_Lcurve = Button(root,text = "L-Curve",command = LCurve)
button_Lcurve.place(height=25, width =80, x = 250, y = 150)

#Create the button for Multiple Angle
button_Multiple = Button(root,text = 'Multi-Angle',command = MultiAngle)
button_Multiple.place(height=25,width =100, x =30, y = 260)

#Create the Clear Button
button_clear = Button(root,text="Clear",padx=20, command = Clear)
button_clear.place(height=25,width=60, x=410, y=150)

#Set up the Regularization value holder
Entry_regularization = Entry(root)
Entry_regularization.place(height=25,width=60, x=110, y=225)
label_regularization = Label(root,text="Regularization\nParameter       ")
label_regularization.place(height=40,width=100, x=10, y=215)
label_regularization_colon = Label(root,text = ":")
label_regularization_colon.place(height=40,width=10, x=100, y=215)

#Set Up for the Number of Points holder
Entry_points = Entry(root)
Entry_points.place(height=25,width=110, x=350, y=185)
Entry_points.insert(0,'200')
label_points = Label(root , text="# of Points:")
label_points.place(height=25,width=70, x=280, y=185)

#Set up for a logorithmic size dist or normal
dist_type = [
    'Linear',
    'Logarithmic'
]

#Hold the value that the user selects
select_dist_type = StringVar()
select_dist_type.set('Logarithmic')

#Create the drop down menu and label
Drop_dist_type = OptionMenu(root,select_dist_type,*dist_type,command= DistTypeChanges)
Drop_dist_type.config(width=20)
Drop_dist_type.place(height=25,width=110, x=120, y=185)
label_dist_type = Label(root, text = "Choose Dist Type:")
label_dist_type.place(height=25,width=115, x=5, y=185)

#Set up for Min size
Entry_min = Entry(root)
Entry_min.place(height=25,width=55, x=270, y=225)
label_min = Label(root,text="Min Size(nm):")
label_min.place(height=25,width=90, x=180, y=225)

#Set Up for Max size
Entry_max = Entry(root)
Entry_max.place(height=25,width=55, x=420, y=225)
label_max = Label(root,text="Max Size(nm):")
label_max.place(height=25,width=90, x=330, y=225)

#Choose Beta
Entry_beta = Entry(root)
Entry_beta.place(height=25,width=150, x=550, y=25)
label_beta = Label(root,text="Choose Beta:")
label_beta.place(height=25,width=150, x=550, y=0)
Entry_beta.bind('<Key>',handle_wait_Beta)

#Enter the Low text entery
Entry_low = Entry(root)
Entry_low.place(height=25,width=150, x=700, y=25)
label_low = Label(root,text="Low tau:")
label_low.place(height=25,width=150, x=700, y=0)
Entry_low.bind('<Key>',handle_wait_Beta)

#Enter the High tau text entery
Entry_high = Entry(root)
Entry_high.place(height=25,width=150, x=850, y=25)
label_high = Label(root,text="High tau:")
label_high.place(height=25,width=150, x=850, y=0)
Entry_high.bind('<Key>',handle_wait_Beta)

#Enter button to allow for update of plots
#Update_button = Button(root,text='Update',command = Update)
#Update_button.place(height=25,width=100, x=850, y=25)

#Graph for the Size distribution
global Figure_g2, Plot_g2, canvas_g2

#Initialize the G2 plot
Figure_g2 = Figure(figsize=(5, 4),dpi=85)
#Width is x and y is height
Plot_g2 = Figure_g2.add_subplot(1,15 , (1,14))
Plot_g2.set_xscale('log')
Plot_g2.plot(tau, g2, color='red')
Plot_g2.set_title("Data")
Plot_g2.set_xlabel("Times (s)")
canvas_g2 = FigureCanvasTkAgg(Figure_g2, root)
canvas_g2.get_tk_widget().place(height=400,width=500, x=510, y=50)

#Initialize the Error Plot
global Figure_Error, Plot_error, canvas_Error
Figure_Error = Figure(dpi=85)
Plot_error = Figure_Error.add_subplot(2,15 , (1,14))
Plot_error.set_xscale('log')
#Plot_error.set_yscale('symlog')
Plot_error.tick_params(axis="y" , labelsize=9)
Plot_error.plot(tau, g2, 'ro', fillstyle='none')
Plot_error.set_title("Normal Error")
Plot_error.set_xlabel(" Times (s)")
Plot_error.yaxis.tick_left()
canvas_Error = FigureCanvasTkAgg(Figure_Error, root)
canvas_Error.get_tk_widget().place(height=400,width=500, x=510, y=450)

global Figure_prob,Plot_prob,canvas_prob

#Initialize the Probility Plot
Figure_prob = Figure(figsize=(5, 4),dpi=85)
Plot_prob = Figure_prob.add_subplot(1, 1, 1)
Plot_prob.plot(tau, g2, color='red')
Plot_prob.set_title("Probability (Intensity)")
Plot_prob.set_xlabel("Particle Size (nm)")
canvas_prob = FigureCanvasTkAgg(Figure_prob, root)
canvas_prob.get_tk_widget().place(height=415,width=510, x=0, y=285)

#Display the Copyright
label_high = Label(root,text="(C) 2023 M. Salazar, H. Srivastav, A. Srivastava, S. Srivastava")
label_high.place(height=25,width=510, x=0, y=700)


#Keep the GUI looping to account for changes in entries
root.mainloop()
