import sys
sys.path.append('pythonfiles/')
from homepage import HomePage

try:
    # for Python2
    from Tkinter import *
    from ttk import *
except ImportError:
    # for Python3
    from tkinter import *
    from tkinter.ttk import *

if __name__ == "__main__":
    root = Tk()
    root.title("ScEpi Search")
    root.resizable(0,0)
    w, h = root.winfo_screenwidth(), root.winfo_screenheight()
    root.geometry("%dx%d+0+0" % (w, h))
    app = HomePage(root)
    app.AddContents()
    root.mainloop()
    exit()
