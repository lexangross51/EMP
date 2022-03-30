import json
import pygame
from pygame.locals import *
from OpenGL.GL import *
from OpenGL.GLU import *

def input(elems, faces, points):
    with open("mesh\mesh.json", "r") as file:
        data = json.load(file)

    for nodes in data['elems']:
        elems.append(nodes)

    for point in data['points']:
        points.append(point)

    for face in data['faces']:
        faces.append(face)

def draw_mesh(elems, points):
    verticies = []
    edges = []

    for elem in elems:
        for i in range(len(elem)):
            verticies.append(points[i])




verticies = (
    (1, -1, -1),
    (1, 1, -1),
    (-1, 1, -1),
    (-1, -1, -1),
    (1, -1, 1),
    (1, 1, 1),
    (-1, -1, 1),
    (-1, 1, 1)
    )
edges = (
    (0,1),
    (0,3),
    (0,4),
    (2,1),
    (2,3),
    (2,7),
    (6,3),
    (6,4),
    (6,7),
    (5,1),
    (5,4),
    (5,7)
    )
def Cube():
    glBegin(GL_LINES)
    for edge in edges:
        for vertex in edge:
            glVertex3fv(verticies[vertex])
    glEnd()


def main():
    elems = []
    points = []
    faces = []


    input(elems, faces, points)

    pygame.init()
    display = (800,600)
    pygame.display.set_mode(display, DOUBLEBUF|OPENGL)
    gluPerspective(45, (display[0]/display[1]), 0.1, 50.0)
    glTranslatef(0.0,0.0, -10)
    while True:
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                pygame.quit()
                quit()
        glRotatef(1, 3, 1, 1)
        glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT)

        draw_mesh(elems, points)

        pygame.display.flip()
        pygame.time.wait(10)


if __name__ == "__main__":
    main()