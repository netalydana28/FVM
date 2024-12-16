#second step animation separate graphs one 
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
#Pressure subplot
P = Figure(figsize = (4,4), dpi = 100)
p = P.add_subplot(111)
#U-speed
U = Figure(figsize = (4,4), dpi = 100)
u = U.add_subplot(111)
#V-speed
V = Figure(figsize = (4,4), dpi = 100)
v = V.add_subplot(111)
#Concentration
C = Figure(figsize = (4,4), dpi = 100)
c = C.add_subplot(111)


def animate_P(self, i):
    p.clear()
    f = open("out3M500K.dat", "r")
    X_data = []
    Y_data = []
    P_data = []
    U_data = []
    V_data = []
    C_data = []

    for x in f:
        list = [float(i) for i in x.split()]
        X_data.append(list[0])
        Y_data.append(list[1])
        P_data.append(list[2])
        U_data.append(list[3])
        V_data.append(list[4])
        C_data.append(list[5])

    f.close()
    print(f"Pressure {P_data[0]}\n")
    xi, yi = np.linspace(max(X_data), min(X_data), 101), np.linspace(max(Y_data), min(Y_data), 101)
    xi, yi = np.meshgrid(xi, yi)
    rbf = scipy.interpolate.Rbf(X_data, Y_data, P_data, function = 'linear')
    pi = rbf(xi, yi)
    self.p.imshow(pi, vmin=min(self.P_data), vmax=max(self.P_data), origin = 'lower', extent=[min(self.X_data), max(self.X_data), min(self.Y_data), max(self.Y_data)])
    self.plot_p = self.p.scatter(self.X_data, self.Y_data, c=self.P_data, vmin=min(self.P_data), vmax=max(self.P_data))
    # try: 
    #     self.cb_P.remove()
    # except: 
    #     pass
    # else:
    #     self.cb_P.remove()
    # self.cb_P = plt.colorbar(self.plot_p, ax=self.p)
    # self.canvas_P.delete("all")
    
    self.p.set_title("Pressure")
    

def animate_U(self, i):
    u.clear()
    f = open("out3M500K.dat", "r")
    X_data = []
    Y_data = []
    P_data = []
    U_data = []
    V_data = []
    C_data = []

    for x in f:
        list = [float(i) for i in x.split()]
        X_data.append(list[0])
        Y_data.append(list[1])
        P_data.append(list[2])
        U_data.append(list[3])
        V_data.append(list[4])
        C_data.append(list[5])

    f.close()
    rbf = scipy.interpolate.Rbf(self.X_data, self.Y_data, self.U_data, function = 'linear')
    ui = rbf(self.xi, self.yi)
    self.u.imshow(ui, vmin=min(self.U_data), vmax=max(self.U_data), origin = 'lower', extent=[min(self.X_data), max(self.X_data), min(self.Y_data), max(self.Y_data)])
    self.plot_u = self.u.scatter(self.X_data, self.Y_data, c=self.U_data, vmin=min(self.U_data), vmax=max(self.U_data))
    # try: 
    #     self.cb_U.remove()
    # except: 
    #     pass
    # else:
    #     self.cb_U.remove()
    # self.cb = plt.colorbar(self.plot_u, ax=self.u)
    # self.canvas_U.delete("all")
    

def animate_V(self, i):
    v.clear()
    f = open("out3M500K.dat", "r")
    X_data = []
    Y_data = []
    P_data = []
    U_data = []
    V_data = []
    C_data = []

    for x in f:
        list = [float(i) for i in x.split()]
        X_data.append(list[0])
        Y_data.append(list[1])
        P_data.append(list[2])
        U_data.append(list[3])
        V_data.append(list[4])
        C_data.append(list[5])

    f.close()
    rbf = scipy.interpolate.Rbf(self.X_data, self.Y_data, self.V_data, function = 'linear')
    vi = rbf(self.xi, self.yi)
    self.v.imshow(vi, vmin=min(self.V_data), vmax=max(self.V_data), origin = 'lower', extent=[min(self.X_data), max(self.X_data), min(self.Y_data), max(self.Y_data)])
    self.plot_v = self.v.scatter(self.X_data, self.Y_data, c=self.V_data, vmin=min(self.V_data), vmax=max(self.V_data))
    # try: 
    #     self.cb_V.remove()
    # except: 
    #     pass
    # else:
    #     self.cb_V.remove()
    # self.cb = plt.colorbar(self.plot_v, ax=self.v)
    # self.canvas_V.delete("all")

def animate_C(self, i):
    c.clear()
    f = open("out3M500K.dat", "r")
    X_data = []
    Y_data = []
    P_data = []
    U_data = []
    V_data = []
    C_data = []

    for x in f:
        list = [float(i) for i in x.split()]
        X_data.append(list[0])
        Y_data.append(list[1])
        P_data.append(list[2])
        U_data.append(list[3])
        V_data.append(list[4])
        C_data.append(list[5])

    f.close()
    rbf = scipy.interpolate.Rbf(self.X_data, self.Y_data, self.C_data, function = 'linear')
    ci = rbf(self.xi, self.yi)
    self.c.imshow(ci, vmin=min(self.C_data), vmax=max(self.C_data), origin = 'lower', extent=[min(self.X_data), max(self.X_data), min(self.Y_data), max(self.Y_data)])
    self.plot_c = self.c.scatter(self.X_data, self.Y_data, c=self.C_data, vmin=min(self.C_data), vmax=max(self.C_data))
    # try: 
    #     self.cb_C.remove()
    # except: 
    #     pass
    # else:
    #     self.cb_C.remove()
    # self.cb = plt.colorbar(self.plot_c, ax=self.c)
    # self.canvas_C.delete("all")

class interface():

    def __init__(self):
        self.window = Tk()
        self.frame = Frame()
        self.window.title("Navier Stokes by FVM")
        self.window.geometry("1536x864")
        self.entry_page = Frame(self.window)
        self.main_page = Frame(self.window)
        self.entry_page.pack()
        self.read_data()
        #labels
        entry_title = Label(self.entry_page, text = "Solving Navier Stokes by FVM", font = ("Arial Bold" , 20))
        main_title = Label(self.main_page, text = "Change parameter", font=("Arial Bond", 12))
        Label(self.entry_page).grid(row=0, column=0)
        Label(self.main_page).grid(row=0, column=0)
        entry_title.grid(row=1, column=0)
        main_title.grid(row=1, column=0)
        #buttons 
        self.to_start = Button(self.entry_page, text="Start", font = ("Arial Bold", 12), command = self.change_frame).grid(row=2, column = 0)
        #setting up a regular grid of interpolation points
        # app_P = SeaofBTCapp()
        # ani_P = animation.FuncAnimation(self.P, self.animate_P, interval=1000)
        # # app_P.mainloop()
        # # app_U = SeaofBTCapp()
        # ani_U = animation.FuncAnimation(self.U, self.animate_U, interval=1000)
        # # app_U.mainloop()
        # # app_V = SeaofBTCapp()
        # ani_V = animation.FuncAnimation(self.V, self.animate_V, interval=1000)
        # # app_V.mainloop()
        # ani_C = animation.FuncAnimation(self.C, self.animate_C, interval=1000)
        # app_C.mainloop()
        self.animate_P(1)
        self.animate_U(1)
        self.animate_V(1)
        self.animate_C(1)
        # app_P = SeaofBTCapp()
        ani_P = animation.FuncAnimation(self.P, self.animate_P, interval=100)
        # ani_P.save('pressure.mp4')
        self.canvas_P = FigureCanvasTkAgg(self.P, self.main_page)
        self.canvas_P.draw()
        # self.p.quiver(self.X_data, self.Y_data, self.U_data, self.V_data)
        self.canvas_P.get_tk_widget().grid(row=3,column=0)
        # # app_P.mainloop()
        # # app_U = SeaofBTCapp()
        ani_U = animation.FuncAnimation(self.U, self.animate_U, interval=1000)
        # ani_U.save('uspeed.mp4')
        self.canvas_U = FigureCanvasTkAgg(self.U, self.main_page)
        self.u.set_title("U-speed")
        self.canvas_U.draw()
        # self.u.quiver(self.X_data, self.Y_data, self.U_data, self.V_data)
        self.canvas_U.get_tk_widget().grid(row=3,column=1)
        # # app_U.mainloop()
        # # app_V = SeaofBTCapp()
        ani_V = animation.FuncAnimation(self.V, self.animate_V, interval=1000)
        # ani_V.save('vspeed.mp4')
        # # app_V.mainloop()
        self.canvas_V = FigureCanvasTkAgg(self.V, self.main_page)
        self.v.set_title("V-speed")
        self.canvas_V.draw()
        # self.v.quiver(self.X_data, self.Y_data, self.U_data, self.V_data)
        self.canvas_V.get_tk_widget().grid(row=4,column=0)
        ani_C = animation.FuncAnimation(self.C, self.animate_C, interval=1000)
        self.canvas_C = FigureCanvasTkAgg(self.C, self.main_page)
        self.c.set_title("Concentration")
        self.canvas_C.draw()
        # self.c.quiver(self.X_data, self.Y_data, self.U_data, self.V_data)
        self.canvas_C.get_tk_widget().grid(row=4,column=1)
        # ani_C.save('concentrarion.mp4')
        # app_C.mainloop()
        # plt.show()


        # self.window.mainloop()
    
    def change_frame(self):
        self.entry_page.forget()
        self.main_page.pack()

    

if __name__ == "__main__":
    start = interface()
    start.mainloop()