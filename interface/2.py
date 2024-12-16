#second step animation
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk as NavigationToolbar2TkAgg
from matplotlib.figure import Figure
from tkinter import *
import numpy as np
import scipy.interpolate
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time

class interface:
    window = Tk()
    frame = Frame()

    def __init__(self):
        self.window.title("Navier Stokes by FVM")
        self.window.geometry("1536x864")
    
    def change_frame(self):
        self.entry_page.forget()
        self.main_page.pack()

    def read_data(self):
        f = open("out3M500K.dat", "r")
        self.X_data = []
        self.Y_data = []
        self.P_data = []
        self.U_data = []
        self.V_data = []
        self.C_data = []

        for x in f:
            list = [float(i) for i in x.split()]
            self.X_data.append(list[0])
            self.Y_data.append(list[1])
            self.P_data.append(list[2])
            self.U_data.append(list[3])
            self.V_data.append(list[4])
            self.C_data.append(list[5])

        f.close()
    
    def animate(self):
        f = open("out3M500K.dat", "r").read()
        self.P_data = []
        self.U_data = []
        self.V_data = []
        self.C_data = []

        for x in f:
            list = [float(i) for i in x.split()]
            self.P_data.append(list[2])
            self.U_data.append(list[3])
            self.V_data.append(list[4])
            self.C_data.append(list[5])
        #graph 
        #interpolation
            #Pressure
        rbf = scipy.interpolate.Rbf(self.X_data, self.Y_data, self.P_data, function = 'linear')
        pi = rbf(self.xi, self.yi)
        self.p.imshow(pi, vmin=min(self.P_data), vmax=max(self.P_data), origin = 'lower', extent=[min(self.X_data), max(self.X_data), min(self.Y_data), max(self.Y_data)])
        self.plot_p = self.p.scatter(self.X_data, self.Y_data, c=self.P_data, vmin=min(self.P_data), vmax=max(self.P_data))
            #U-speed
        rbf = scipy.interpolate.Rbf(self.X_data, self.Y_data, self.U_data, function = 'linear')
        ui = rbf(self.xi, self.yi)
        self.u.imshow(ui, vmin=min(self.U_data), vmax=max(self.U_data), origin = 'lower', extent=[min(self.X_data), max(self.X_data), min(self.Y_data), max(self.Y_data)])
        self.plot_u = self.u.scatter(self.X_data, self.Y_data, c=self.U_data, vmin=min(self.U_data), vmax=max(self.U_data))
            #V-speed
        rbf = scipy.interpolate.Rbf(self.X_data, self.Y_data, self.V_data, function = 'linear')
        vi = rbf(self.xi, self.yi)
        self.v.imshow(vi, vmin=min(self.V_data), vmax=max(self.V_data), origin = 'lower', extent=[min(self.X_data), max(self.X_data), min(self.Y_data), max(self.Y_data)])
        self.plot_v = self.v.scatter(self.X_data, self.Y_data, c=self.V_data, vmin=min(self.V_data), vmax=max(self.V_data))
            #C-speed
        rbf = scipy.interpolate.Rbf(self.X_data, self.Y_data, self.C_data, function = 'linear')
        ci = rbf(self.xi, self.yi)
        self.c.imshow(ci, vmin=min(self.C_data), vmax=max(self.C_data), origin = 'lower', extent=[min(self.X_data), max(self.X_data), min(self.Y_data), max(self.Y_data)])
        self.plot_c = self.c.scatter(self.X_data, self.Y_data, c=self.C_data, vmin=min(self.C_data), vmax=max(self.C_data))
        if(self.current_parameter == "P"):
            self.cb = plt.colorbar(self.plot_p, ax=self.p)
            self.canvas = FigureCanvasTkAgg(self.P, self.main_page)
            self.canvas.draw()
            self.p.quiver(self.X_data, self.Y_data, self.U_data, self.V_data)
            self.canvas.get_tk_widget().grid(row=3,column=0)
        elif(self.current_parameter == "U"):
            self.cb = plt.colorbar(self.plot_u, ax=self.u)
            self.canvas = FigureCanvasTkAgg(self.U, self.main_page)
            self.canvas.draw()
            self.p.quiver(self.X_data, self.Y_data, self.U_data, self.V_data)
            self.canvas.get_tk_widget().grid(row=3,column=0)
        elif(self.current_parameter == "V"):
            self.cb = plt.colorbar(self.plot_v, ax=self.v)
            self.canvas = FigureCanvasTkAgg(self.V, self.main_page)
            self.canvas.draw()
            self.p.quiver(self.X_data, self.Y_data, self.U_data, self.V_data)
            self.canvas.get_tk_widget().grid(row=3,column=0)
        elif(self.current_parameter == "C"):
            self.cb = plt.colorbar(self.plot_c, ax=self.c)
            self.canvas = FigureCanvasTkAgg(self.C, self.main_page)
            self.canvas.draw()
            self.p.quiver(self.X_data, self.Y_data, self.U_data, self.V_data)
            self.canvas.get_tk_widget().grid(row=3,column=0)
       
    def change_parameter(self):
        for item in self.canvas.get_tk_widget().find_all():
            self.canvas.get_tk_widget().delete(item)
        var = self.parameter.get()
        self.current_parameter = var
        self.cb.remove()
        if(var == "U"):
            self.canvas = FigureCanvasTkAgg(self.U, self.main_page)
            self.cb = plt.colorbar(self.plot_u, ax=self.u)
            self.canvas.draw()
            self.u.quiver(self.X_data, self.Y_data, self.U_data, self.V_data)
        elif(var == "V"):
            self.canvas = FigureCanvasTkAgg(self.V, self.main_page)
            self.cb = plt.colorbar(self.plot_v, ax=self.v)
            self.canvas.draw()
            self.v.quiver(self.X_data, self.Y_data, self.U_data, self.V_data)
        elif(var == "C"):
            self.canvas = FigureCanvasTkAgg(self.C, self.main_page)
            self.cb = plt.colorbar(self.plot_c, ax=self.c)
            self.canvas.draw()
            self.c.quiver(self.X_data, self.Y_data, self.U_data, self.V_data)
        elif(var == "P"):
            self.canvas = FigureCanvasTkAgg(self.P, self.main_page)
            self.cb = plt.colorbar(self.plot_p, ax=self.p)
            self.canvas.draw()
            self.p.quiver(self.X_data, self.Y_data, self.U_data, self.V_data)
        self.canvas.get_tk_widget().grid(row=3,column=0)

    def main_frame(self):
        self.animate(1)
        self.entry_page = Frame(self.window)
        self.main_page = Frame(self.window)
        self.entry_page.pack()
        #labels
        entry_title = Label(self.entry_page, text = "Solving Navier Stokes by FVM", font = ("Arial Bold" , 20))
        main_title = Label(self.main_page, text = "Change parameter", font=("Arial Bond", 12))
        Label(self.entry_page).grid(row=0, column=0)
        Label(self.main_page).grid(row=0, column=0)
        entry_title.grid(row=1, column=0)
        main_title.grid(row=1, column=0)
        #buttons 
        self.to_start = Button(self.entry_page, text="Start", font = ("Arial Bold", 12), command = self.change_frame).grid(row=2, column = 0)
        #spinbox
        parameters = ["P", "U", "V", "C"]
        spinbox_var = StringVar()
        self.parameter = Spinbox(self.main_page, textvariable = spinbox_var, values = parameters, command=self.change_parameter)
        self.parameter.grid(row=2, column = 0)
        #setting up a regular grid of interpolation points
        self.xi, self.yi = np.linspace(max(self.X_data), min(self.X_data), 101), np.linspace(max(self.Y_data), min(self.Y_data), 101)
        self.xi, self.yi = np.meshgrid(self.xi, self.yi)
        #Pressure subplot
        self.P = Figure(figsize = (6,6), dpi = 100)
        self.p = self.P.add_subplot(111)
        #U-speed
        self.U = Figure(figsize = (6,6), dpi = 100)
        self.u = self.U.add_subplot(111)
        #V-speed
        self.V = Figure(figsize = (6,6), dpi = 100)
        self.v = self.V.add_subplot(111)
        #Concentration
        self.C = Figure(figsize = (6,6), dpi = 100)
        self.c = self.C.add_subplot(111)
        self.current_parameter = "P"
        # app = SeaofBTCapp()
        # ani = animation.FuncAnimation()
        while True: 
            animate()
            print("ok")
            time.sleep(1)

        self.window.mainloop()

if __name__ == "__main__":
    start = interface()
    start.main_frame()