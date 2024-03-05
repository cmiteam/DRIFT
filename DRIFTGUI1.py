import os
import csv
import tkinter as tk
from tkinter import ttk
from tkinter import font as tkFont
from DRIFT1 import run_population_model 

data_directory = "data"

class PopulationModelingApp:

    def __init__(self, master):
        self.master = master
        
        master.title("DRIFT 1.0")
        frame_info = {
            'main_parameter_frame': {'title': 'Main Parameters and Settings', 'group': 'main', 'placement': 'top_left'},
            'DNA_parameter_frame': {'title': 'DNA Parameters and Settings', 'group': 'DNA', 'placement': 'top_right'},
            'mutation_parameter_frame': {'title': 'Mutation Parameters and Settings', 'group': 'mutation', 'placement': 'bottom_right'}
        }
        self.initial_parameter_states = self.load_defaults('parameterdefaults.csv')
        self.initial_plot_states = self.load_defaults('plotdefaults.csv')
        self.parameter_states = {}
        self.plot_states = {}
        self.check_var_dict = {}

        for frame_name, info in frame_info.items():
            frame = self.create_parameter_frame(master, info['title'], info['placement'])
            setattr(self, frame_name, frame)
            self.populate_frame(frame, self.initial_parameter_states, info['group'])

        self.toggle_frames()
        self.update_plot_variables()
        tk.Label(self.main_parameter_frame, text="All Set?", font=tkFont.Font(size=12, weight="bold")).grid(row=32, column=0, columnspan=4)
        self.run_button = tk.Button(self.main_parameter_frame, text="Run Model", command=lambda: self.run_model(self.parameter_states, self.plot_states, self.initial_parameter_states, self.check_var_dict))
        self.run_button.grid(row=33, column=0, columnspan=4, pady=(10, 10))
        self.run_button.config(state="normal", bg="green", fg="white")

    def create_parameter_frame(self, master, title, placement):
        frame = tk.Frame(master, bd=2, relief="ridge")
        if placement == 'top_left':
            frame.place(x=25, y=25, anchor="nw")
        elif placement == 'top_right':
            frame.place(x=450, y=25, anchor="nw")
        elif placement == 'bottom_right':
            frame.place(x=450, y=380, anchor="nw")
        tk.Label(frame, text=title, font=tkFont.Font(size=14, weight="bold", underline=True)).grid(row=0, column=0, columnspan=4)
        tk.Label(frame, text="Model Variables", font=tkFont.Font(size=12, weight="bold")).grid(row=1, column=0, columnspan=2)
        tk.Label(frame, text="Plot Variables", font=tkFont.Font(size=12, weight="bold")).grid(row=1, column=2, columnspan=2)
        return frame

    def populate_frame(self, frame, states, group):
        widget_row = 2
        for key, value in states.items():
            if value.get("group") == group:
                if value['type'] == 'Text':
                    tk.Label(frame, text=value["label"] + ":").grid(row=widget_row, column=0)
                    entry = tk.Entry(frame, text="", width=20)
                    entry.grid(row=widget_row, column=1)
                    entry.insert(0, value["default"])
                    self.parameter_states[key] = entry
                    widget_row += 1
                elif value['type'] == 'Check':
                    tk.Label(frame, text=value["label"] + ":").grid(row=widget_row, column=0)
                    var = tk.IntVar()
                    var.set(value['default'])
                    check = tk.Checkbutton(frame, text='', variable=var, command=self.toggle_frames, width=20)
                    check.grid(row=widget_row, column=1, pady=0)
                    self.parameter_states[key] = check
                    self.check_var_dict[key] = var
                    widget_row += 1
                elif value['type'] == 'Dropdown':
                    tk.Label(frame, text=value["label"] + ":").grid(row=widget_row, column=0)
                    options = value['default'].split(',')
                    dropdown = ttk.Combobox(frame, values=options, state='readonly', width=17)
                    dropdown.set(options[0])
                    dropdown.grid(row=widget_row, column=1)
                    self.parameter_states[key] = dropdown
                    widget_row += 1
        widget_row = 2
        for key, value in self.initial_plot_states.items():
            if value.get("group") == group:
                tk.Label(frame, text=value["label"] + ":").grid(row=widget_row, column=2)
                var = tk.BooleanVar()
                var.set(value['default'])
                check = tk.Checkbutton(frame, variable=var, text="", command=self.update_plot_variables)
                check.grid(row=widget_row, column=3)
                if self.initial_plot_states[key]['default'] == True:
                    check.select()
                self.check_var_dict[key] = var
                widget_row += 1
        tk.Label(frame, text='').grid(row=25, column=0, columnspan=4)

    def toggle_frames(self):
        if self.check_var_dict['track_DNA'].get():
            for widget in self.DNA_parameter_frame.winfo_children():
                widget.configure(state="normal")
        else:
            for widget in self.DNA_parameter_frame.winfo_children():
                widget.configure(state="disabled")
        if self.check_var_dict['track_mutations'].get():
            for widget in self.mutation_parameter_frame.winfo_children():
                widget.configure(state="normal")
        else:
            for widget in self.mutation_parameter_frame.winfo_children():
                widget.configure(state="disabled")

    def update_plot_variables(self):
        for key in self.check_var_dict:
            if key in self.initial_plot_states:
                self.plot_states[key] = self.check_var_dict[key].get()

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

    def load_defaults(self, csv_file):
        defaults = {}
        filename = os.path.join(data_directory, csv_file)
        with open(filename, newline='') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                parameter = row['parameter']
                defaults[parameter] = {
                    'label': row['label'],
                    'type': row['type'],
                    'value_format': row['value_format'],
                    'default': row['default'],
                    'alt': row['alt'] if 'alt' in row else None,
                    'group': row['group']
                }
        return defaults

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

    def run_model(self, parameter_states, plot_states, initial_parameter_states, check_var_dict):

        self.run_button.config(state="disabled", bg="red")
        model_parameters = {}
        model_parameters["plot"] = []

        for key, widget in parameter_states.items():
            widget_type = self.initial_parameter_states[key]['type']
            if widget_type == 'Text':
                format = self.initial_parameter_states[key]['value_format']
                value = ''
                if format == 'int':
                    value = int(widget.get())
                elif format == 'float':
                    value = float(widget.get())
                else:
                    value = widget.get()
                model_parameters[key] = value
            elif widget_type == 'Check':
                value  = self.check_var_dict[key].get()
                model_parameters[key] = value
            elif widget_type == 'Dropdown':
                value = widget.get()
                model_parameters[key] = value

        for key, state in self.plot_states.items():
            if state == 1:
                model_parameters["plot"].append(key)

        result = run_population_model(model_parameters)
        self.run_button.config(state="normal", bg="green")

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
root = tk.Tk()
root.geometry("900x750") # width x height
app = PopulationModelingApp(root)
root.mainloop()