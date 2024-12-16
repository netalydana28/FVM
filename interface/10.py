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
        f = open("C:\\Users\\Eldana\\OneDrive - АО Казахстанско-Британский Технический Университет\\Documents\\Project\\fvm\\optimization1\\out.dat", "r")
        self.X_data.clear()
        self.Y_data.clear()
        self.U_data.clear()
        self.P_data.clear()
        self.V_data.clear()
        self.C_data.clear()

        for x in f:
            list = [float(i) for i in x.split()]
            self.X_data.append(list[0])
            self.Y_data.append(list[1])
            self.P_data.append(list[2])
            self.U_data.append(list[3])
            self.V_data.append(list[4])
            self.C_data.append(list[5])

        f.close()
    
    

    
    def animate(self, i):
        self.p.clear()
        self.read_data()
        # print(f"Pressure {self.P_data[0]}\n")
        rbf = scipy.interpolate.Rbf(self.X_data, self.Y_data, self.P_data, function = 'linear')
        pi = rbf(self.xi, self.yi)
        self.img_p.set_data(pi)
        self.plot_p = self.p.scatter(self.X_data, self.Y_data, c=self.P_data,   vmin=min(self.P_data), vmax=max(self.P_data))
        self.vectors_p.set_UVC(self.U_data, self.V_data)
        self.vectors_p = self.p.quiver(self.X_data, self.Y_data, self.U_data, self.V_data)
        self.cb_p.remove()
        self.cb_p = plt.colorbar(self.plot_p, ax=self.p)
        self.u.clear()
        rbf = scipy.interpolate.Rbf(self.X_data, self.Y_data, self.U_data, function = 'linear')
        ui = rbf(self.xi, self.yi)
        self.img_u.set_data(ui)
        self.plot_u = self.u.scatter(self.X_data, self.Y_data, c=self.U_data,   vmin=min(self.U_data), vmax=max(self.U_data))
        self.vectors_u.set_UVC(self.U_data, self.V_data)
        self.vectors_u = self.u.quiver(self.X_data, self.Y_data, self.U_data, self.V_data)
        self.cb_u.remove()
        self.cb_u = plt.colorbar(self.plot_u, ax=self.u)
        self.v.clear()
        rbf = scipy.interpolate.Rbf(self.X_data, self.Y_data, self.V_data, function = 'linear')
        vi = rbf(self.xi, self.yi)
        self.img_v.set_data(vi)
        self.plot_v = self.v.scatter(self.X_data, self.Y_data, c=self.V_data,   vmin=min(self.V_data), vmax=max(self.V_data))
        self.vectors_v.set_UVC(self.U_data, self.V_data)
        self.vectors_v = self.v.quiver(self.X_data, self.Y_data, self.U_data, self.V_data)
        self.cb_v.remove()
        self.cb_v = plt.colorbar(self.plot_v, ax=self.v)
        self.c.clear()
        rbf = scipy.interpolate.Rbf(self.X_data, self.Y_data, self.C_data, function = 'linear')
        ci = rbf(self.xi, self.yi)
        self.img_c.set_data(ci)
        self.plot_c = self.c.scatter(self.X_data, self.Y_data, c=self.C_data,   vmin=min(self.C_data), vmax=max(self.C_data))
        self.vectors_c.set_UVC(self.U_data, self.V_data)
        self.vectors_c = self.c.quiver(self.X_data, self.Y_data, self.U_data, self.V_data)
        self.cb_c.remove()
        self.cb_c = plt.colorbar(self.plot_c, ax=self.c)
        self.p.set_title("Pressure")
        self.u.set_title("U-speed")
        self.v.set_title("V-speed")
        self.c.set_title("Concentration")

    def main_frame(self):
        self.entry_page = Frame(self.window)
        self.main_page = Frame(self.window)
        self.entry_page.pack()
        self.X_data = []
        self.Y_data = []
        self.P_data = []
        self.U_data = []
        self.V_data = []
        self.C_data = []
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
        self.xi, self.yi = np.linspace(max(self.X_data), min(self.X_data), 101), np.linspace(max(self.Y_data), min(self.Y_data), 101)
        self.xi, self.yi = np.meshgrid(self.xi, self.yi)
        #Pressure subplot
        self.Graph = Figure(figsize = (8,8), dpi = 100)
        self.p = self.Graph.add_subplot(221)
        rbf = scipy.interpolate.Rbf(self.X_data, self.Y_data, self.P_data, function = 'linear')
        pi = rbf(self.xi, self.yi)
        self.img_p = self.p.imshow(pi, vmin=min(self.P_data), vmax=max(self.P_data), origin = 'lower', extent=[min(self.X_data), max(self.X_data), min(self.Y_data), max(self.Y_data)])
        self.plot_p = self.p.scatter(self.X_data, self.Y_data, c=self.P_data,   vmin=min(self.P_data), vmax=max(self.P_data))
        self.cb_p = plt.colorbar(self.plot_p, ax=self.p)
        self.vectors_p = self.p.quiver(self.X_data, self.Y_data, self.U_data, self.V_data)
        #U-speed
        # self.U = Figure(figsize = (4,4), dpi = 100)
        self.u = self.Graph.add_subplot(222)
        rbf = scipy.interpolate.Rbf(self.X_data, self.Y_data, self.U_data, function = 'linear')
        ui = rbf(self.xi, self.yi)
        self.img_u = self.u.imshow(ui, vmin=min(self.U_data), vmax=max(self.U_data), origin = 'lower', extent=[min(self.X_data), max(self.X_data), min(self.Y_data), max(self.Y_data)])
        self.plot_u = self.u.scatter(self.X_data, self.Y_data, c=self.U_data,   vmin=min(self.U_data), vmax=max(self.U_data))
        self.cb_u = plt.colorbar(self.plot_u, ax=self.u)
        # self.canvas_U.draw()
        self.vectors_u = self.u.quiver(self.X_data, self.Y_data, self.U_data, self.V_data)
        #V-speed
        self.v = self.Graph.add_subplot(223)
        rbf = scipy.interpolate.Rbf(self.X_data, self.Y_data, self.V_data, function = 'linear')
        vi = rbf(self.xi, self.yi)
        self.img_v = self.v.imshow(vi, vmin=min(self.V_data), vmax=max(self.V_data), origin = 'lower', extent=[min(self.X_data), max(self.X_data), min(self.Y_data), max(self.Y_data)])
        self.plot_v = self.v.scatter(self.X_data, self.Y_data, c=self.V_data,   vmin=min(self.V_data), vmax=max(self.V_data))
        self.cb_v = plt.colorbar(self.plot_v, ax=self.v)
        self.vectors_v = self.v.quiver(self.X_data, self.Y_data, self.U_data, self.V_data)
        #Concentration
        self.c = self.Graph.add_subplot(224)
        rbf = scipy.interpolate.Rbf(self.X_data, self.Y_data, self.C_data, function = 'linear')
        ci = rbf(self.xi, self.yi)
        self.img_c = self.c.imshow(ci, vmin=min(self.C_data), vmax=max(self.C_data), origin = 'lower', extent=[min(self.X_data), max(self.X_data), min(self.Y_data), max(self.Y_data)])
        self.plot_c = self.c.scatter(self.X_data, self.Y_data, c=self.C_data,   vmin=min(self.C_data), vmax=max(self.C_data))
        self.cb_c = plt.colorbar(self.plot_c, ax=self.c)
        self.vectors_c = self.c.quiver(self.X_data, self.Y_data, self.U_data, self.V_data)
        self.canvas = FigureCanvasTkAgg(self.Graph, self.main_page)
        self.p.set_title("Pressure")
        self.u.set_title("U-speed")
        self.v.set_title("V-speed")
        self.c.set_title("Concentration")
        self.canvas.draw()
        self.canvas.get_tk_widget().grid(row=3,column=0)

        ani = animation.FuncAnimation(self.Graph, self.animate, interval=500)
        plt.show()
        # Writer = animation.writers['ffmpeg']
        # writer = Writer(fps=20, metadata=dict(artist='Me'), bitrate=1800)
        # ani.save('animation.mp4', writer = writer)
        self.window.mainloop()
        

if __name__ == "__main__":
    start = interface()
    start.main_frame()