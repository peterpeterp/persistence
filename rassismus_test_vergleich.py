import Tkinter as tk
import os, sys, subprocess, ConfigParser,time, string
import numpy as np
from random import random,randrange
from time import sleep










def spiel():
	ini()
	global hoehe,breite,skal,canvas,imgtk,root2
	root2 = tk.Tk()  
	root2.geometry('745x295+800+400')
	root2.grid()


	canvas = tk.Canvas(root2, width=breite/skal, height=hoehe/skal)
	canvas.create_rectangle(0,0,200,100,fill="red")

	#canvas.bind("<Button-1>",lambda event:links(event,1))
	#canvas.bind("<Button-3>",lambda event:rechts(event,0))
	canvas.grid(row=2,rowspan=5,column=1,columnspan=12)
	root2.mainloop() 

spiel()