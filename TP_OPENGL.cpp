#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>

#include <GL/freeglut.h>

#include "Point.h"
#include "Vector.h"

// define window
#define WIDTH 500
#define HEIGHT 500
#define RED 0.25
#define GREEN 0.25
#define BLUE 0.25
#define ALPHA 1
// for leaving the program
#define KEY_ESC 27

using namespace std;

void init_scene();
void render_scene();
GLvoid initGL();
GLvoid window_display();
GLvoid window_reshape(GLsizei width, GLsizei height);
GLvoid window_key(unsigned char key, int x, int y);
void SphereClavier(int key, int x, int y);
void DisplayVoxel(Point centre, double length);

long nbMeridienParallele = 8;

void drawCylindre(Point** cylindre, long nbMeridien);
Point** facetCylindre(double rayon, double hauteur, long nbMeridien);
void drawCone(Point* Cone, long nbMeridien);
Point* facetCone(long nbMeridien);
void drawSphere(Point** Sphere, long nbMeridienParallele);
Point** facet_Sphere(long rayon);
void displayVoxel(Point centre, double length);
void StdDisplayVoxel(Point* tab);

int main(int argc, char** argv) {

  //Glut initialisation
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_RGBA);
  glutInitWindowSize(WIDTH, HEIGHT);
  glutInitWindowPosition(0, 0);
  glutCreateWindow("GLUT Testing");

  //OpenGL initialisation
  initGL();
  init_scene();


  glutDisplayFunc(&window_display);
  glutReshapeFunc(&window_reshape);
  glutKeyboardFunc(&window_key);
  glutSpecialFunc(SphereClavier);
  glutMainLoop();

  return 1;
}

GLvoid initGL() {
  glClearColor(RED, GREEN, BLUE, ALPHA);
}

// Stock var for future scene, if needed. Leave empty if no var used.
void init_scene() {}

GLvoid window_display() {
  glClear(GL_COLOR_BUFFER_BIT);
  glLoadIdentity();

  //rendering the scene here : it's there for the functions
  render_scene();
  //putting the rendered scene on screen
  glFlush();
}

GLvoid window_reshape(GLsizei width, GLsizei height) {
  glViewport(0, 0, width, height);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  // size of the window
  glOrtho(-30.0, 30.0, -30.0, 30.0, -30.0, 30.0);

  //defining the gluLookAt with values... might have it changed later with variables
  double xCam, yCam, zCam;
  double xCentre, yCentre, zCentre;
  double xVue, yVue, zVue;
  xCam = 2.0;
  yCam = 1.0, zCam = 1.0;
  xCentre = 0.0;
  yCentre = 2.0, zCentre = 0.0;
  xVue = 1.0;
  yVue = 0.0, zVue = 2.0;
  gluLookAt(xCam, yCam, zCam, xCentre, yCentre, zCentre, xVue, yVue, zVue);

  //applying the gluLookAt
  glMatrixMode(GL_MODELVIEW);
}

  //keyboard events
GLvoid window_key(unsigned char key, int x, int y) {
  switch (key) {
    case KEY_ESC:
      exit(1);
      break;
    case 43:
      nbMeridienParallele++;
      glClear(GL_COLOR_BUFFER_BIT);
      glLoadIdentity();
      render_scene();
      glFlush();
      break;
    case 45:
      if (nbMeridienParallele/2 > 3) {
        nbMeridienParallele--;
      } else {
        cout << "Impossible de diminuer encore le nombre de méridiens et parallèles"
             << endl;
      }
      glClear(GL_COLOR_BUFFER_BIT);
      glLoadIdentity();
      render_scene();
      glFlush();
      break;

    default:
      printf("La touche %d n´est pas active.\n", key);
      break;
  }
}

void SphereClavier(int key, int x, int y) {
  if (key == GLUT_KEY_UP) {
    nbMeridienParallele++;
  } else if (key == GLUT_KEY_DOWN) {
    if (nbMeridienParallele > 3) {
      nbMeridienParallele--;  // Nombre réel : nbMeridiens*2
    } else {
      cout << "Impossible de dimunuer plus le nombre de méridiens et de parallèles"
           << endl;
    }
  }

  glClear(GL_COLOR_BUFFER_BIT);
  glLoadIdentity();
  render_scene();
  glFlush();
}

//That's the kicking : this function needs to call everyhting important
void render_scene() {
  glColor3f(0.0, 1.0, 1.0);

  // Drawing cylinder
  /*
      double rayon = 10.0;
      double hauteur = 20.0;
      long nbMeridien = 10;

      Point** Cylindre = facetCylindre(rayon, hauteur, nbMeridien);
      drawCylindre(Cylindre, nbMeridien);
    
   // Drawing Cone
  /*
        long nbMeridien = 10.0;

        Point* Cone = facetCone(nbMeridien);
        drawCone(Cone, nbMeridien);
  */

  // Drawing sphere
  /*
	  long rayon = 20;
	  Point** Sphere = facet_Sphere(rayon);
	  drawSphere(Sphere, nbMeridienParallele);
  */
}

Point** facet_Sphere(long rayon) 
{
  double teta = 2 * M_PI / (nbMeridienParallele);
  double phi = M_PI / (nbMeridienParallele + 1);
  double xU, yU, zU;

  // nb de parallele + les deux sommets
  Point** paralleles = new Point*[nbMeridienParallele + 2];
  // paralleles[j] represente un parrallele
  // paralleles[j][i] represente un point sur un parrallele
  Point poleNord = Point(0.0, 0.0, rayon);
  paralleles[0] = new Point[1];
  paralleles[0][0] = poleNord;

  int cpt = 0;
  for (int j = 1; j < nbMeridienParallele + 1; j++) {
    cpt = j;

    paralleles[j] = new Point[nbMeridienParallele + 1];
    double cosinus_Phi = cos(phi * j);  // cos prend des radians
    double sinus_Phi = sin(phi * j);

    for (int i = 0; i < nbMeridienParallele; i++) {
      double cosinus_Teta = cos(teta * i);
      double sinus_Teta = sin(teta * i);

      xU = rayon * cosinus_Teta * sinus_Phi;
      yU = rayon * sinus_Teta * sinus_Phi;
      zU = rayon * cosinus_Phi;
      paralleles[j][i] = Point(xU, yU, zU);
    }
    paralleles[j][nbMeridienParallele] = paralleles[j][0];
  }

  cpt++;
  Point poleSud = Point(0.0, 0.0, -rayon);
  paralleles[cpt] = new Point[1];
  paralleles[cpt][0] = poleSud;

  return paralleles;
}

void drawSphere(Point** Sphere, long nbMeridienParallele) {

  int cpt = 0;
  for (int j = 1; j < nbMeridienParallele + 1; j++) {
    cpt = j;
    drawCurve(Sphere[j], nbMeridienParallele + 1);

    // North of the sphere
    if (j == 1) {
      for (int i = 0; i < nbMeridienParallele; i++) {
        glBegin(GL_LINE_STRIP);
        glVertex3f(Sphere[0][0].getX(), Sphere[0][0].getY(),
                   Sphere[0][0].getZ());
        glVertex3f(Sphere[j][i].getX(), Sphere[j][i].getY(),
                   Sphere[j][i].getZ());
        glEnd();
      }
    }

    // Main body
    if (j < nbMeridienParallele) {
      for (int i = 0; i < nbMeridienParallele; i++) {
        glBegin(GL_LINE_STRIP);
        glVertex3f(Sphere[j][i].getX(), Sphere[j][i].getY(),
                   Sphere[j][i].getZ());
        glVertex3f(Sphere[j + 1][i].getX(), Sphere[j + 1][i].getY(),
                   Sphere[j + 1][i].getZ());
        glEnd();
      }
    }

    // South of the sphere
    if (j == nbMeridienParallele) {
      for (int i = 0; i < nbMeridienParallele; i++) {
        glBegin(GL_LINE_STRIP);
        glVertex3f(Sphere[j][i].getX(), Sphere[j][i].getY(),
                   Sphere[j][i].getZ());
        glVertex3f(Sphere[j + 1][0].getX(), Sphere[j + 1][0].getY(),
                   Sphere[j + 1][0].getZ());
        glEnd();
      }
    }
  }
  cpt++;
}

Point* facetCone(long nbMeridien) {
  double rayon_Base = 15.0;
  Point sommet(0.0, 0.0, 20.0);
  double hauteur = 20.0;

  Point* Cone = new Point[nbMeridien + 1];  // summet + other values
  Cone[0] = sommet;

  double tetaU = 2 * M_PI / nbMeridien;
  double xU, yU, zU;

  for (int i = 0; i < nbMeridien; i++) {
    double cosinus_Teta = cos(tetaU * i);  // cos needs radians
    double sinus_Teta = sin(tetaU * i);
    if ((tetaU * i) == M_PI) {
      cosinus_Teta = -1.0;
      sinus_Teta = 0.0;
    }
    xU = rayon_Base * cosinus_Teta;
    yU = rayon_Base * sinus_Teta;
    zU = sommet.getZ() - hauteur;

    Cone[i + 1] = Point(xU, yU, zU);
  }

  return Cone;
}

void drawCone(Point* Cone, long nbMeridien) {
  // glBegin(GL_POLYGON);
  glBegin(GL_LINE_STRIP);
  for (int i = 1; i <= nbMeridien; i++) {
    glVertex3f(Cone[i].getX(), Cone[i].getY(), Cone[i].getZ());
  }
  glVertex3f(Cone[1].getX(), Cone[1].getY(), Cone[1].getZ());
  glEnd();

  for (int i = 1; i <= nbMeridien; i++) {
    glBegin(GL_LINE_STRIP);
    glVertex3f(Cone[0].getX(), Cone[0].getY(), Cone[0].getZ());
    glVertex3f(Cone[i].getX(), Cone[i].getY(), Cone[i].getZ());
    glEnd();
  }
  glColor3f(1.0, 0.0, 0.0);
  Cone[0].affiche();
}

Point** facetCylindre(double rayon, double hauteur, long nbMeridien) {
  Point** Cercles = new Point*[2];
  Cercles[0] = new Point[nbMeridien];
  Cercles[1] = new Point[nbMeridien];

  double tetaU = 2 * M_PI / nbMeridien;
  double xU, yU, zU;
  for (int i = 0; i < nbMeridien; i++) {
    double cosinus_Teta = cos(tetaU * i);  // cos needs radians
    double sinus_Teta = sin(tetaU * i);
    if ((tetaU * i) == M_PI) {
      cosinus_Teta = -1.0;
      sinus_Teta = 0.0;
    }
    xU = rayon * cosinus_Teta;
    yU = rayon * sinus_Teta;
    zU = hauteur / 2;

    Point haut(xU, yU, zU);
    Cercles[0][i] = haut;
    // Cercles[0][i].affiche();

    Point bas(xU, yU, -zU);
    Cercles[1][i] = bas;
  }

  return Cercles;
}

void drawCylindre(Point** cylindre, long nbMeridien) {
  Point* CercleHaut = cylindre[0];
  Point* CercleBas = cylindre[1];

  glBegin(GL_POLYGON);
  for (int i = 0; i < nbMeridien; i++) {
    glVertex3f(CercleHaut[i].getX(), CercleHaut[i].getY(),
               CercleHaut[i].getZ());
    cout << CercleHaut[i].getX() << " " << CercleHaut[i].getY() << " " << CercleHaut[i].getZ() << endl;
  }
  glEnd();

  glColor3f(0.5, 0.5, 0.0);
  glBegin(GL_POLYGON);
  for (int i = 0; i < nbMeridien; i++) {
    glVertex3f(CercleBas[i].getX(), CercleBas[i].getY(), CercleBas[i].getZ());
  }
  glEnd();

  for (int i = 0; i < nbMeridien; i++) {
    glBegin(GL_LINE_STRIP);
    glVertex3f(CercleHaut[i].getX(), CercleHaut[i].getY(),
               CercleHaut[i].getZ());
    glVertex3f(CercleBas[i].getX(), CercleBas[i].getY(), CercleBas[i].getZ());
    glEnd();
  }
}

void StdDisplayVoxel(Point* tab)
{
  glBegin(GL_POINTS);
  for (int i = 0; i < 8; ++i)
  {
    glVertex3d(tab[i].getX(), tab[i].getY(), tab[i].getZ());
  }
  glEnd();

  /*glColor3f(1,0,1);
    glBegin(GL_QUADS);
        glVertex3d(tab[0].getX(), tab[0].getY(), tab[0].getZ());
        glVertex3d(tab[1].getX(), tab[1].getY(), tab[1].getZ());
        glVertex3d(tab[5].getX(), tab[5].getY(), tab[5].getZ());
        glVertex3d(tab[4].getX(), tab[4].getY(), tab[4].getZ());
  
        glVertex3d(tab[1].getX(), tab[1].getY(), tab[1].getZ());
        glVertex3d(tab[3].getX(), tab[3].getY(), tab[3].getZ());
        glVertex3d(tab[7].getX(), tab[7].getY(), tab[7].getZ());
        glVertex3d(tab[5].getX(), tab[5].getY(), tab[5].getZ());
   
        glVertex3d(tab[4].getX(), tab[4].getY(), tab[4].getZ());
        glVertex3d(tab[5].getX(), tab[5].getY(), tab[5].getZ());
        glVertex3d(tab[7].getX(), tab[7].getY(), tab[7].getZ());
        glVertex3d(tab[6].getX(), tab[6].getY(), tab[6].getZ());
  
        glVertex3d(tab[0].getX(), tab[0].getY(), tab[0].getZ());
        glVertex3d(tab[2].getX(), tab[2].getY(), tab[2].getZ());
        glVertex3d(tab[6].getX(), tab[6].getY(), tab[6].getZ());
        glVertex3d(tab[4].getX(), tab[4].getY(), tab[4].getZ());
   
        glVertex3d(tab[0].getX(), tab[0].getY(), tab[0].getZ());
        glVertex3d(tab[1].getX(), tab[1].getY(), tab[1].getZ());
        glVertex3d(tab[3].getX(), tab[3].getY(), tab[3].getZ());
        glVertex3d(tab[2].getX(), tab[2].getY(), tab[2].getZ());
   
        glVertex3d(tab[2].getX(), tab[2].getY(), tab[2].getZ());
        glVertex3d(tab[3].getX(), tab[3].getY(), tab[3].getZ());
        glVertex3d(tab[7].getX(), tab[7].getY(), tab[7].getZ());
        glVertex3d(tab[6].getX(), tab[6].getY(), tab[6].getZ());
    glEnd();*/

  glBegin(GL_LINES);
    glVertex3d(tab[0].getX(), tab[0].getY(), tab[0].getZ());
    glVertex3d(tab[1].getX(), tab[1].getY(), tab[1].getZ());

    glVertex3d(tab[0].getX(), tab[0].getY(), tab[0].getZ());
    glVertex3d(tab[2].getX(), tab[2].getY(), tab[2].getZ());

    glVertex3d(tab[0].getX(), tab[0].getY(), tab[0].getZ());
    glVertex3d(tab[4].getX(), tab[4].getY(), tab[4].getZ());

    glVertex3d(tab[1].getX(), tab[1].getY(), tab[1].getZ());
    glVertex3d(tab[3].getX(), tab[3].getY(), tab[3].getZ());

    glVertex3d(tab[1].getX(), tab[1].getY(), tab[1].getZ());
    glVertex3d(tab[5].getX(), tab[5].getY(), tab[5].getZ());

    glVertex3d(tab[2].getX(), tab[2].getY(), tab[2].getZ());
    glVertex3d(tab[3].getX(), tab[3].getY(), tab[3].getZ());

    glVertex3d(tab[2].getX(), tab[2].getY(), tab[2].getZ());
    glVertex3d(tab[6].getX(), tab[6].getY(), tab[6].getZ());

    glVertex3d(tab[3].getX(), tab[3].getY(), tab[3].getZ());
    glVertex3d(tab[7].getX(), tab[7].getY(), tab[7].getZ());

    glVertex3d(tab[4].getX(), tab[4].getY(), tab[4].getZ());
    glVertex3d(tab[5].getX(), tab[5].getY(), tab[5].getZ());

    glVertex3d(tab[4].getX(), tab[4].getY(), tab[4].getZ());
    glVertex3d(tab[6].getX(), tab[6].getY(), tab[6].getZ());

    glVertex3d(tab[5].getX(), tab[5].getY(), tab[5].getZ());
    glVertex3d(tab[7].getX(), tab[7].getY(), tab[7].getZ());

    glVertex3d(tab[6].getX(), tab[6].getY(), tab[6].getZ());
    glVertex3d(tab[7].getX(), tab[7].getY(), tab[7].getZ());
      glEnd();

}

Point ** createVoxel(Point* Tab, Point centre, double rayon)
{
  Point ** TabVox = new Point*[9];
  for (int i = 0; i < 9; ++i)
  {
    TabVox[i] = new Point[9];
  }

  Point milieu(centre.getX(), centre.getY(), centre.getZ());
  Point Sommet0(Tab[0].getX(),Tab[0].getY(), Tab[0].getZ());
  Point Sommet1(Tab[1].getX(),Tab[1].getY(), Tab[1].getZ());
  Point Sommet2(Tab[2].getX(),Tab[2].getY(), Tab[2].getZ());
  Point Sommet3(Tab[3].getX(),Tab[3].getY(), Tab[3].getZ());
  Point Sommet4(Tab[4].getX(),Tab[4].getY(), Tab[4].getZ());
  Point Sommet5(Tab[5].getX(),Tab[5].getY(), Tab[5].getZ());
  Point Sommet6(Tab[6].getX(),Tab[6].getY(), Tab[6].getZ());
  Point Sommet7(Tab[7].getX(),Tab[7].getY(), Tab[7].getZ());

  Point Arete01(Tab[0].getX(), Tab[0].getY() + rayon, Tab[0].getZ());
  Point Arete02(Tab[0].getX() + rayon, Tab[0].getY(), Tab[0].getZ());
  Point Arete04(Tab[0].getX(), Tab[0].getY(), Tab[0].getZ() - rayon);
  Point Arete13(Tab[1].getX() + rayon, Tab[1].getY(), Tab[1].getZ());
  Point Arete15(Tab[1].getX(), Tab[1].getY(), Tab[1].getZ() - rayon);
  Point Arete23(Tab[2].getX(), Tab[2].getY() + rayon, Tab[2].getZ());
  Point Arete26(Tab[2].getX(), Tab[2].getY(), Tab[2].getZ() - rayon);
  Point Arete37(Tab[3].getX(), Tab[3].getY(), Tab[3].getZ() - rayon);
  Point Arete45(Tab[4].getX(), Tab[4].getY() + rayon, Tab[4].getZ());
  Point Arete46(Tab[4].getX() + rayon, Tab[4].getY(), Tab[4].getZ());
  Point Arete57(Tab[5].getX() + rayon, Tab[5].getY(), Tab[5].getZ());
  Point Arete67(Tab[6].getX(), Tab[6].getY() + rayon, Tab[6].getZ());

  Point Face0132(centre.getX(), centre.getY(), centre.getZ() + rayon);
  Point Face1375(centre.getX(), centre.getY() + rayon, centre.getZ());
  Point Face2376(centre.getX() + rayon, centre.getY(), centre.getZ());
  Point Face0154(centre.getX() - rayon, centre.getY(), centre.getZ());
  Point Face4576(centre.getX(), centre.getY(), centre.getZ() - rayon);
  Point Face0264(centre.getX(), centre.getY() - rayon, centre.getZ());


TabVox[0][0] = Sommet0;
TabVox[0][1] = Arete01;
TabVox[0][2] = Arete02;
TabVox[0][3] = Face0132;
TabVox[0][4] = Arete04;
TabVox[0][5] = Face0154;
TabVox[0][6] = Face0264;
TabVox[0][7] = milieu;

TabVox[1][0] = Arete01;
TabVox[1][1] = Sommet1;
TabVox[1][2] = Face0132;
TabVox[1][3] = Arete13;
TabVox[1][4] = Face0154;
TabVox[1][5] = Arete15;
TabVox[1][6] = milieu;
TabVox[1][7] = Face1375;

TabVox[2][0] = Arete02;
TabVox[2][1] = Face0132;
TabVox[2][2] = Sommet2;
TabVox[2][3] = Arete23;
TabVox[2][4] = Face0264;
TabVox[2][5] = milieu;
TabVox[2][6] = Arete26;
TabVox[2][7] = Face2376;

TabVox[3][0] = Face0132;
TabVox[3][1] = Arete13;
TabVox[3][2] = Arete23;
TabVox[3][3] = Sommet3;
TabVox[3][4] = milieu;
TabVox[3][5] = Face1375;
TabVox[3][6] = Face2376;
TabVox[3][7] = Arete37;

TabVox[4][0] = Arete04;
TabVox[4][1] = Face0154;
TabVox[4][2] = Face0264;
TabVox[4][3] = milieu;
TabVox[4][4] = Sommet4;
TabVox[4][5] = Arete45;
TabVox[4][6] = Arete46;
TabVox[4][7] = Face4576;


TabVox[5][0] = Face0154;
TabVox[5][1] = Arete15;
TabVox[5][2] = milieu;
TabVox[5][3] = Face1375;
TabVox[5][4] = Arete45;
TabVox[5][5] = Sommet5;
TabVox[5][6] = Face4576;
TabVox[5][7] = Arete57;

TabVox[6][0] = Face0264;
TabVox[6][1] = milieu;
TabVox[6][2] = Arete26;
TabVox[6][3] = Face2376;
TabVox[6][4] = Arete46;
TabVox[6][5] = Face4576;
TabVox[6][6] = Sommet6;
TabVox[6][7] = Arete67;

TabVox[7][0] = milieu;
TabVox[7][1] = Face1375;
TabVox[7][2] = Face2376;
TabVox[7][3] = Arete37;
TabVox[7][4] = Face4576;
TabVox[7][5] = Arete57;
TabVox[7][6] = Arete67;
TabVox[7][7] = Sommet7;

TabVox[0][8] = Point(milieu.getX() - rayon/2, milieu.getY() - rayon/2, milieu.getZ() + rayon/2); //0
TabVox[1][8] = Point(milieu.getX() - rayon/2, milieu.getY() + rayon/2, milieu.getZ() + rayon/2); //1
TabVox[2][8] = Point(milieu.getX() + rayon/2, milieu.getY() - rayon/2, milieu.getZ() + rayon/2); //2 
TabVox[3][8] = Point(milieu.getX() + rayon/2, milieu.getY() + rayon/2, milieu.getZ() + rayon/2); //3
TabVox[4][8] = Point(milieu.getX() - rayon/2, milieu.getY() - rayon/2, milieu.getZ() - rayon/2); //4
TabVox[5][8] = Point(milieu.getX() - rayon/2, milieu.getY() + rayon/2, milieu.getZ() - rayon/2); //5
TabVox[6][8] = Point(milieu.getX() + rayon/2, milieu.getY() - rayon/2, milieu.getZ() - rayon/2); //6
TabVox[7][8] = Point(milieu.getX() + rayon/2, milieu.getY() + rayon/2, milieu.getZ() - rayon/2);

return TabVox;

}

double Distance(Point A, Point B)
{
  double Xbis = (A.getX() - B.getX()) * (A.getX() - B.getX());
  double Ybis = (A.getY() - B.getY()) * (A.getY() - B.getY());
  double Zbis = (A.getZ() - B.getZ()) * (A.getZ() - B.getZ());
  double somme = sqrt(Xbis + Ybis + Zbis);
  return somme;
}

void CreateSpherebyVoxel(Point centre, double length, double resolution)
{
  Point pts[8] = {
            Point(centre.getX() - length, centre.getY() - length, centre.getZ() + length), //0
            Point(centre.getX() - length, centre.getY() + length, centre.getZ() + length), //1
            Point(centre.getX() + length, centre.getY() - length, centre.getZ() + length), //2 
            Point(centre.getX() + length, centre.getY() + length, centre.getZ() + length), //3
            Point(centre.getX() - length, centre.getY() - length, centre.getZ() - length), //4
            Point(centre.getX() - length, centre.getY() + length, centre.getZ() - length), //5
            Point(centre.getX() + length, centre.getY() - length, centre.getZ() - length), //6
            Point(centre.getX() + length, centre.getY() + length, centre.getZ() - length) //7
          };
  Point** TabVoxel = createVoxel(pts,centre, length);
  for (int i = 0; i < 8; ++i)
  {
    GenerateSphere(TabVoxel[i], centre, length, resolution,length);
  }
}

void GenerateSphere(Point* TabVoxel, Point centre, double length, double resolution,  double currentR)
{
  int cpt = 0;
  for (int i = 0; i < 8; ++i)
  {
    if (Distance(TabVoxel[i], centre) <= length)
    {
      cpt++;
    }
  }
  if (cpt >= 8)
  {
    StdDisplayVoxel(TabVoxel);
  }
  else if(currentR >= resolution)
  {
    if (cpt > 0)
    {
      Point** SubVoxel = createVoxel(TabVoxel,TabVoxel[8], currentR/2);
      for (int i = 0; i < 8; ++i)
      {
        GenerateSphere(SubVoxel[i], centre, length, resolution, currentR/2);
      }
    }
  }
}

//Newton polygon calcul
int polyNewton(int n, int i) {
  return f(n) / (f(i) * f(n - i));
}

//factorielle, iterative style
int f(int n) {
  int facto = 1;
  for (int i = 1; i <= n; i++) {
    facto *= i;
  }
  return facto;
}
