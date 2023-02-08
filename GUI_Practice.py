#Practice from CodeCamp Use Kinter

#Practice One Create Hello World

#Import tkinter materials
from tkinter import *

#Initialize the widget packages
root = Tk()

#Create a label widget
mylabel = Label(root, text = "Hello World!")

#Place it into a screen near closest space
mylabel.pack()

#Event loop (a while loop)
root.mainloop()
