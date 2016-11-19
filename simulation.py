from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *
import numpy as np
import matplotlib.pyplot as plt
import sys
import math
import threading
import time
shoulder = 30
elbow = 30
dotShoulder = 0
dotElbow = 0
i = 0
j = 0
result=[]
result2=[]
f = open("trainSample.txt")
f2 = open("trainSample.txt")
def init():
    newThread = threading.Thread(target = darwPlot, args = ())
    newThread.start()
    glClearColor(0.0, 0.0, 0.0, 0.0)
    gluOrtho2D(-5.0, 5.0, -5.0, 5.0)
    glLoadIdentity()
    glShadeModel(GL_FLAT)

def drawcoordinate():

    glPointSize(3.0)
    glColor3f(1.0, 1.0, 0.0)
    glBegin(GL_LINES)
    glVertex2f(-5.0, 0.0)
    glVertex2f(5.0, 0.0)
    glVertex2f(0.0, 5.0)
    glVertex2f(0.0, -5.0)
    glEnd()

def UpdateAngle():
    global shoulder
    global elbow
    global i
    line = f.readline()
    time.sleep(0.2)
    if line:
        result.append(map(float,line.split(' ')))
        shoulder = 180*result[i][0]/math.pi
        elbow = 180*result[i][1]/math.pi
        i = i+1
        glutPostRedisplay()  
    
def darwPlot():
    global j
    global i
    plt.figure(1)
    plt.axis([0, 200, -3,3])
    ax1 = plt.subplot(211)
    ax2 = plt.subplot(212)
    line = f2.readline()
    while line:
        result2.append(map(float,line.split(' ')))
        plt.sca(ax1)
        plt.scatter(j, result2[j][0])
        plt.sca(ax2)
        plt.scatter(j, result2[j][1])
        plt.pause(0.001)
        line = f2.readline()
        j = j+1
    plt.show()

def display():
    
    glClear (GL_COLOR_BUFFER_BIT);  
    drawcoordinate()
    glColor3f(1.0, 0.0, 1.0)
    glPushMatrix() 
    glTranslatef (0.0, 0.0, 0.0) 
    glRotatef (shoulder, 0.0, 0.0, 1.0)  
    glTranslatef (1.0, 0.0, 0.0)  
    glPushMatrix() 
    glScalef (2.0, 0.1, 0.1)  
    glutSolidCube (1.0) 
    glPopMatrix()  
    glColor3f(0.0, 1.0, 1.0)
    glTranslatef (1.0, 0.0, 0.0)  
    glRotatef (elbow, 0.0, 0.0, 0.4)  
    glTranslatef (1.0, 0.0, 0.0)  
    glPushMatrix()  
    glScalef (2.0, 0.1, 0.1)
    glutSolidCube (1.0)  
    glPopMatrix() 
    glPopMatrix()
    glutSwapBuffers()


def reshape(w, h):
    glViewport (0, 0, w, h);   
    glMatrixMode (GL_PROJECTION);  
    glLoadIdentity ();  
    gluPerspective(65.0, w/h, 1.0, 20.0);  
    glMatrixMode(GL_MODELVIEW);  
    glLoadIdentity();  
    glTranslatef (0.0, 0.0, -8.0);  

def keyboard(key, x, y):
    pass
def main():
    glutInit(sys.argv)
    glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA)
    glutInitWindowPosition(100,100)
    glutInitWindowSize(500,500)
    glutCreateWindow("Arm")
    glutDisplayFunc(display)
    glutIdleFunc(UpdateAngle)
    glutReshapeFunc(reshape) 
    glutKeyboardFunc(keyboard)
    init()
    glutMainLoop()
main()