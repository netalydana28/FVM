import numpy as np
from tkinter import *
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import matplotlib.animation as animation 
from matplotlib.widgets import LassoSelector
from matplotlib.path import Path
import scipy.interpolate


class interface:
    window = Tk()
    frame = Frame()

    def __init__(self):
        self.window.title("Navier Stokes by FVM")
        self.window.geometry("1536x864")
    
    def change_frame1(self):
        self.entry_page.forget()
        self.set_C_page.pack()

    def change_frame2(self):
        f = open("C:\\Users\\Eldana\\OneDrive - АО Казахстанско-Британский Технический Университет\\Documents\\Project\\fvm\\optimization1\\concentration.txt", "w")
        txt = ""
        for i in range(len(self.ind)):
            txt+= str(self.ind[i]) + " "
        f.write(txt)
        f.close()
        self.set_C_page.forget()
        self.set_U_page.pack()
    
    def change_frame4(self):
        txt = ""
        for i in range(len(self.vis_list)):
            txt += " ".join(self.vis_list[i])+"\n"
        f = open("C:\\Users\\Eldana\\OneDrive - АО Казахстанско-Британский Технический Университет\\Documents\\Project\\fvm\\optimization1\\uspeed.txt", "w")
        f.write(txt)
        f.close()
        self.set_U_page.forget()
        self.choice_page.pack()

    def change_frame3(self):
        self.parameters = self.e_choice.get().split()
        par_dict = {'P': 0, 'U': 1, 'V':2, 'C':3}
        for i in range(len(self.parameters)):
            self.parameters[i] = par_dict[self.parameters[i]]
        self.choice_page.forget()
            #setting up a regular grid of interpolation points
        self.Graph = Figure(figsize = (8,8), dpi = 100)
        self.xi, self.yi = np.linspace(max(self.X_data), min(self.X_data), 101), np.linspace(max(self.Y_data), min(self.Y_data), 101)
        self.xi, self.yi = np.meshgrid(self.xi, self.yi)
        functions = [[self.draw_P, self.animate_P], [self.draw_U, self.animate_U], [self.draw_V, self.animate_V], [self.draw_C, self.animate_C]]
        if(len(self.parameters)>2):
            subplots = 220
        else:
            subplots = 100 + len(self.parameters)*10
        for i in range(len(self.parameters)):
            if (self.parameters[i] == 0):
                self.p = self.Graph.add_subplot(subplots+i+1)
            elif (self.parameters[i] == 1):
                self.u = self.Graph.add_subplot(subplots+i+1)
            elif (self.parameters[i] == 2):
                self.v = self.Graph.add_subplot(subplots+i+1)
            elif (self.parameters[i] == 3):
                self.c = self.Graph.add_subplot(subplots+i+1)
        self.functions = np.array(functions)
        for i in range(len(self.parameters)):
            self.functions[self.parameters[i]][0]()
        self.Graph.tight_layout()
        #Setting canvas
        self.canvas = FigureCanvasTkAgg(self.Graph, self.vis_page)
        self.canvas.draw()
        self.canvas.get_tk_widget().grid(row=3,column=0)
        self.vis_page.pack()
        ani = animation.FuncAnimation(self.Graph, self.animate, blit = False, interval=400)
        plt.show()

    def check(self):
        self.degrees = float(self.e_degrees.get())
        self.points_X = np.array([0, 1, 0, 1])
        self.points_Y = np.array([0, 0, 1, 1])
        self.initial_u = np.tile(np.array([np.cos((self.degrees)* np.pi/180)])*5, 4)
        self.initial_v = np.tile(np.array([np.sin((self.degrees)* np.pi/180)])*5, 4)
        self.mesh_U = Figure(figsize = (6,6), dpi = 100)
        self.Mesh_U = self.mesh_U.add_subplot(111) 
        self.Mesh_U.scatter(self.X_data, self.Y_data, s = 20)
        self.Mesh_U.quiver(self.points_X, self.points_Y, self.initial_u, self.initial_v)
        self.mesh_U.tight_layout()
        self.canvas_mesh_U = FigureCanvasTkAgg(self.mesh_U, self.set_U_page)
        self.canvas_mesh_U.draw()
        self.canvas_mesh_U.get_tk_widget().grid(row=3,column=0, rowspan = 8, columnspan = 8)

    def read_data(self):
        f = open("C:\\Users\\Eldana\\OneDrive - АО Казахстанско-Британский Технический Университет\\Documents\\Project\\fvm\\optimization1\\out5.dat", "r")
        P_data = []
        U_data = []
        V_data = []
        C_data = []

        for x in f:
            list = [float(i) for i in x.split()]
            P_data.append(list[2])
            U_data.append(list[3])
            V_data.append(list[4])
            C_data.append(list[5])

        f.close()

        self.P_data = np.array(P_data)
        self.U_data = np.array(U_data)
        self.V_data = np.array(V_data)
        self.C_data = np.array(C_data) 

    def draw_P(self):
        rbf = scipy.interpolate.Rbf(self.X_data, self.Y_data, self.P_data, function = 'linear')
        pi = rbf(self.xi, self.yi)
        self.img_p = self.p.imshow(pi, vmin=min(self.P_data), vmax=max(self.P_data), origin = 'lower', extent=[min(self.X_data), max(self.X_data), min(self.Y_data), max(self.Y_data)])
        self.plot_p = self.p.scatter(self.X_data, self.Y_data, c=self.P_data, cmap = self.cmap, vmin=min(self.P_data), vmax=max(self.P_data))
        self.cb_p = plt.colorbar(self.plot_p, ax=self.p, fraction=0.046, pad=0.04)
        self.vectors_p = self.p.quiver(self.X_data, self.Y_data, self.U_data, self.V_data)
        self.p.set_title("Pressure")

    def draw_U(self):
        rbf = scipy.interpolate.Rbf(self.X_data, self.Y_data, self.U_data, function = 'linear')
        ui = rbf(self.xi, self.yi)
        self.img_u = self.u.imshow(ui, vmin=min(self.U_data), vmax=max(self.U_data), origin = 'lower', extent=[min(self.X_data), max(self.X_data), min(self.Y_data), max(self.Y_data)])
        self.plot_u = self.u.scatter(self.X_data, self.Y_data, c=self.U_data, cmap = self.cmap, vmin=min(self.U_data), vmax=max(self.U_data))
        self.cb_u = plt.colorbar(self.plot_u, ax=self.u, fraction=0.046, pad=0.04)
        self.vectors_u = self.u.quiver(self.X_data, self.Y_data, self.U_data, self.V_data)
        self.u.set_title("U-speed")

    def draw_V(self):
        rbf = scipy.interpolate.Rbf(self.X_data, self.Y_data, self.V_data, function = 'linear')
        vi = rbf(self.xi, self.yi)
        self.img_v = self.v.imshow(vi, vmin=min(self.V_data), vmax=max(self.V_data), origin = 'lower', extent=[min(self.X_data), max(self.X_data), min(self.Y_data), max(self.Y_data)])
        self.plot_v = self.v.scatter(self.X_data, self.Y_data, c=self.V_data,cmap = self.cmap, vmin=min(self.V_data), vmax=max(self.V_data))
        self.cb_v = plt.colorbar(self.plot_v, ax=self.v, fraction=0.046, pad=0.04)
        self.vectors_v = self.v.quiver(self.X_data, self.Y_data, self.U_data, self.V_data)
        self.v.set_title("V-speed")

    def draw_C(self):   
        rbf = scipy.interpolate.Rbf(self.X_data, self.Y_data, self.C_data, function = 'linear')
        ci = rbf(self.xi, self.yi)
        self.img_c = self.c.imshow(ci, vmin=min(self.C_data), vmax=max(self.C_data), origin = 'lower', extent=[min(self.X_data), max(self.X_data), min(self.Y_data), max(self.Y_data)])
        self.plot_c = self.c.scatter(self.X_data, self.Y_data, c=self.C_data, cmap = self.cmap, vmin=min(self.C_data), vmax=max(self.C_data))
        self.cb_c = plt.colorbar(self.plot_c, ax=self.c, fraction=0.046, pad=0.04)
        self.vectors_c = self.c.quiver(self.X_data, self.Y_data, self.U_data, self.V_data)
        self.c.set_title("Concentration")

    def animate_P(self):
        rbf = scipy.interpolate.Rbf(self.X_data, self.Y_data, self.P_data, function = 'linear')
        pi = rbf(self.xi, self.yi)
        self.img_p.set_data(pi)
        self.plot_p.set_array(self.P_data)
        self.vectors_p.set_UVC(self.U_data, self.V_data)
        self.cb_p.update_normal(self.plot_p)
        self.p.set_title("Pressure")
        return (self.img_p, ), (self.plot_p, ), (self.vectors_p, ), (self.cb_p, )

    def animate_U(self):
        rbf = scipy.interpolate.Rbf(self.X_data, self.Y_data, self.U_data, function = 'linear')
        ui = rbf(self.xi, self.yi)
        self.img_u.set_data(ui)
        self.plot_u.set_offsets(self.coor)
        self.plot_u.set_array(self.U_data)
        self.vectors_u.set_UVC(self.U_data, self.V_data)
        self.cb_u.update_normal(self.plot_u)
        self.u.set_title("U-speed")
        return (self.img_u, ), (self.plot_u, ), (self.vectors_u, ), (self.cb_u, )

    def animate_V(self):
        rbf = scipy.interpolate.Rbf(self.X_data, self.Y_data, self.V_data, function = 'linear')
        vi = rbf(self.xi, self.yi)
        self.img_v.set_data(vi)
        self.plot_v.set_offsets(self.coor)
        self.plot_v.set_array(self.V_data)
        self.vectors_v.set_UVC(self.U_data, self.V_data)
        self.cb_v.update_normal(self.plot_v)
        self.v.set_title("V-speed")
        return (self.img_v, ), (self.plot_v, ), (self.vectors_v, ), (self.cb_v, )

    def animate_C(self):
        rbf = scipy.interpolate.Rbf(self.X_data, self.Y_data, self.C_data, function = 'linear')
        ci = rbf(self.xi, self.yi)
        self.img_c.set_data(ci)
        self.plot_c.set_offsets(self.coor)
        self.plot_c.set_array(self.C_data)
        self.vectors_c.set_UVC(self.U_data, self.V_data)
        self.cb_c.update_normal(self.plot_c)
        self.c.set_title("Concentration")
        return (self.img_c, ), (self.plot_c, ), (self.vectors_c, ), (self.cb_c, )

    def animate(self, i):
        try:
            self.read_data()
            print("animate")
            for i in range(len(self.parameters)):
                self.functions[self.parameters[i]][1]()
        except: 
            pass
        else:
            pass

    def onselect(self, verts):
        path = Path(verts)
        self.ind = np.nonzero(path.contains_points(self.coor))[0]
        self.fc[:, -1] = 0.3
        self.fc[self.ind, -1] = 1
        self.pts.set_facecolors(self.fc)
        self.canvas_mesh.draw_idle()

    def append_data(self):
        self.degrees = str(self.e_degrees.get())
        self.start = str(self.e_start.get())
        self.end = str(self.e_end.get())
        self.text_widget.insert(END, f"{self.degrees}°\t[{self.start}, {self.end}]\n")
        self.vis_list.append([self.degrees, self.start, self.end])

    def main_frame(self):
        self.cmap = plt.get_cmap("jet")
        self.entry_page = Frame(self.window)
        self.set_C_page = Frame(self.window)
        self.set_U_page = Frame(self.window)
        self.choice_page = Frame(self.window)
        self.vis_page = Frame(self.window)
        self.entry_page.pack()
        #labels
        entry_title = Label(self.entry_page, text = "Solving Navier Stokes by FVM", font = ("Arial Bold" , 20))
        set_C_title = Label(self.set_C_page, text = "Set concentration", font = ("Arial Bold", 20))
        set_U_title = Label(self.set_U_page, text = "Set U-speed", font = ("Arial Bold", 20))
        set_U_degrees = Label(self.set_U_page, text = "Enter degrees", font = ("Arial Bold", 12))
        set_start = Label(self.set_U_page, text = "Enter start time ", font = ("Arial Bold", 12))
        set_end = Label(self.set_U_page, text = "end time", font = ("Arial Bold", 12))
        choice_title = Label(self.choice_page, text = "Enter visualization choice \n (Enter with ' ' separator)", font = ("Arial Bold", 20))
        Label(self.entry_page).grid(row=0, column=0)
        Label(self.set_C_page).grid(row=0, column=0)
        Label(self.set_U_page).grid(row=0, column=0)
        Label(self.choice_page).grid(row=0, column=0)
        entry_title.grid(row=1, column=0)
        choice_title.grid(row=1, column=0)
        set_C_title.grid(row=1, column=0)
        set_U_title.grid(row=1, column=0)
        set_U_degrees.grid(row = 5, column=9)
        set_start.grid(row=6, column=9)
        set_end.grid(row=6, column=11)
        #entries 
        self.e_choice= Entry(self.choice_page)
        self.e_choice.grid(row=2, column = 0)
        self.e_degrees = Entry(self.set_U_page)
        self.e_degrees.grid(row=5, column = 10)
        self.e_start = Entry(self.set_U_page)
        self.e_start.grid(row=6, column=10)
        self.e_end = Entry(self.set_U_page)
        self.e_end.grid(row=6, column=12)
        #scroll bar 
        self.scroll_bar = Scrollbar(self.set_U_page, orient = "vertical") # Create a scroll bar
        self.text_widget = Text(self.set_U_page, yscrollcommand=self.scroll_bar.set) # Create a text_widget 
        self.scroll_bar.grid(row = 8, column = 12, rowspan = 3) # Pack the scroll bar
        self.scroll_bar.config(command=self.text_widget.yview) # Attach the scrollbar with the text widget
        self.text_widget.grid(row = 8, column = 9, rowspan = 3, columnspan = 3) # Pack the text widget
        self.text_widget.insert(END, "Degree\tTime interval\n")# Insert text into the text widget
        #buttons 
        self.to_start = Button(self.entry_page, text="Start", font = ("Arial Bold", 12), command = self.change_frame1).grid(row=2, column = 0)
        self.to_set_C = Button(self.set_C_page, text="Set", font = ("Arial Bold", 12), command = self.change_frame2).grid(row=2, column = 0)
        self.to_set_U = Button(self.set_U_page, text="Set", font = ("Arial Bold", 12), command = self.change_frame4).grid(row=2, column = 0)
        self.to_choose = Button(self.choice_page, text="Make choice", font = ("Arial Bold", 12), command = self.change_frame3).grid(row=7, column = 1)
        self.to_check = Button(self.set_U_page, text="Check direction", font = ("Arial", 12), command = self.check).grid(row = 7, column = 9)
        self.to_append = Button(self.set_U_page, text="Append data", font = ("Arial", 12), command = self.append_data).grid(row = 7, column = 11)
        #mesh drawing
            #reading files 
            #coordinates
        f = open("C:\\Users\\Eldana\\OneDrive - АО Казахстанско-Британский Технический Университет\\Documents\\Project\\fvm\\optimization1\\map_coor.txt", "r")
        X_data = []
        Y_data = []
        for x in f:
            list = [float(i) for i in x.split()]
            X_data.append(list[1])
            Y_data.append(list[2])

        f.close()
        self.X_data = np.array(X_data)
        self.Y_data = np.array(Y_data)
        self.coor = np.stack([self.X_data, self.Y_data]).T
            #links
        f = open("C:\\Users\\Eldana\\OneDrive - АО Казахстанско-Британский Технический Университет\\Documents\\Project\\fvm\\optimization1\\map_links.txt", "r")
        triangles = []
        for x in f: 
            list = x.split()
            list.pop(0)
            list.pop(0)
            list.pop(0)
            list = [int(i)-1 for i in list]
            triangles.append(list)

        f.close()
        self.mesh = Figure(figsize = (7,7), dpi = 100)
        self.Mesh = self.mesh.add_subplot(111)
        self.pts = self.Mesh.scatter(self.X_data, self.Y_data, s = 20)
        self.fc = self.pts.get_facecolors()
        self.fc = np.tile(self.fc, (len(self.X_data), 1))
        self.mesh.tight_layout()
        self.canvas_mesh = FigureCanvasTkAgg(self.mesh, self.set_C_page)
        self.canvas_mesh.draw()
        self.canvas_mesh.get_tk_widget().grid(row=3,column=0)
        self.lasso = LassoSelector(self.Mesh, onselect = self.onselect)
        self.vis_list = []
        self.read_data()
        self.window.mainloop()
        
if __name__ == "__main__":
    start = interface()
    start.main_frame()