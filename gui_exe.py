
# This file was generated by the Tkinter Designer by Parth Jadhav
# https://github.com/ParthJadhav/Tkinter-Designer


from pathlib import Path

import tkinter
from sklearn.utils import _typedefs
from sklearn.neighbors import _partition_nodes
from tkinter import *
from tkinter import ttk
from pathlib import Path
from tkinter import filedialog
from tkinter import messagebox
# from tkinter import *
# Explicit imports to satisfy Flake8
from tkinter import Tk, Canvas, Entry, Text, Button, PhotoImage
import analysis_pipeline

import os
import sys

OUTPUT_PATH = Path(__file__).parent
ASSETS_PATH = OUTPUT_PATH / Path("./assets")


def getFolderPath():
    '''
    get selected folder path
    :return:
    '''
    folder_selected = filedialog.askdirectory()
    folderPath.set(folder_selected)

def doStuff():
    '''
    'Update' function, main functions
    '''
    folder = folderPath.get()
    if folder:
        try:
            folder_check = analysis_pipeline.folder_check(folder)
            #print(folder_check)
            if folder_check is None: # pass the check
                replace_status = checkbox_val()
                try:
                    #analysis_pipeline.report_analysis(folder, replace_status) ## comment this line if you don't wish to make changes to database.
                    messagebox.showinfo('Info', 'Process completed! You can close this program now.')
                    print("Doing stuff with folder", folder, replace_status)
                except:
                    messagebox.showerror('Internal Error', 'Error: Algorithm error, please contact UCD MSBA team.')

            else:
                messagebox.showerror('Error', folder_check)
        except ValueError:
            messagebox.showerror('Value Error', 'Error: Folder cannot be opened.')
        except FileNotFoundError:
            messagebox.showerror('File Not Found Error', 'Error: Folder cannot be found.')
    else:
        messagebox.showerror('Folder Not Found Error', 'Error: No folder selected.')


def checkbox_val():
    '''
    :return: if 'replace' box is chosen.
    '''
    is_replace = False
    print(check_var.get())
    if check_var.get():
        is_replace = check_var.get()
    else:
        pass
    return is_replace

def relative_to_assets(path: str) -> Path:
    '''
    path setting utils
    :param path:
    :return:
    '''
    return ASSETS_PATH / Path(path)

def resource_path(relative_path):
    """ Get absolute path to resource, works for dev and for PyInstaller """
    base_path = getattr(sys, '_MEIPASS', os.path.dirname(os.path.abspath(__file__)))
    return os.path.join(base_path, relative_path)

window = Tk()
window.title("Automatic Liquid Biopsy Reporting")
folderPath = StringVar()
window.geometry("937x486")
window.configure(bg = "#FFFFFF")

'''
Page layout
'''
canvas = Canvas(
    window,
    bg = "#FFFFFF",
    height = 956,
    width = 2182,
    bd = 0,
    highlightthickness = 0,
    relief = "ridge"
)

canvas.place(x = 0, y = 0)
canvas.create_rectangle(
    0.0,
    0.0,
    937.0,
    486.0,
    fill="#F9F7F7",
    outline="")

canvas.create_text(
    169.0,
    253.0,
    anchor="nw",
    text="Select a folder to upload:",
    fill="#000000",
    font=("Verdana", 16 * -1)
)

canvas.create_text(
    195.0,
    307.0,
    anchor="nw",
    text="Replace existing files:",
    fill="#000000",
    font=("Verdana", 16 * -1)
)
'''
Checkbox
'''
check_var = tkinter.IntVar()
checkbox = tkinter.Checkbutton(window,variable=check_var, onvalue=True, offvalue=False, bg="#F9F7F7", command=checkbox_val)
checkbox.place(
    x=400.0,
    y=307.0,
    width=54.0,
    height=20.0,

)
'''
Browse folder button
'''
button_image_1 = PhotoImage(
    #file=relative_to_assets("button_1.png"))
    file = resource_path("assets\\button_1.png"))
button_1 = Button(
    image=button_image_1,
    borderwidth=0,
    highlightthickness=0,
    command=getFolderPath,
    relief="flat"
)
button_1.place(
    x=606.0,
    y=242.0,
    width=164.0,
    height=55.0
)

'''
Upload button
'''
button_image_2 = PhotoImage(
    #file=relative_to_assets("button_2.png"))
    file = resource_path("assets\\button_2.png"))
button_2 = Button(
    image=button_image_2,
    borderwidth=0,
    highlightthickness=0,
    command=doStuff,
    relief="flat"
)
button_2.place(
    x=388.0,
    y=351.0,
    width=164.0,
    height=68.0
)

canvas.create_text(
    169.0,
    136.0,
    anchor="nw",
    text="This program will upload Genialis Liquid Biospy reports, identify driving mutations and match patients.\nCheck Tableau updates after clicking Upload",
    fill="#000000",
    font=("Verdana", 13 * -1)
)

'''
Show selected folder path
'''
entry_image_1 = PhotoImage(
    #file=relative_to_assets("entry_1.png"))
    file = resource_path("assets\\entry_1.png"))
entry_bg_1 = canvas.create_image(
    506.0,
    264.5,
    image=entry_image_1
)
entry_1 = Entry(
    bd=0,
    bg="#FFFFFF",
    highlightthickness=0,
    textvariable=folderPath
)
entry_1.place(
    x=394.0,
    y=250.0,
    width=212.0,
    height=27.0
)

canvas.create_rectangle(
    0.0,
    0.0,
    937.0,
    104.0,
    fill="#039392",
    outline="")

canvas.create_text(
    181.0,
    33.0,
    anchor="nw",
    text="Genomic Dashboard - Files Upload",
    fill="#FFFFFF",
    font=("Verdana", 32 * -1)
)
#window.resizable(False, False) #Allow resize
window.mainloop()