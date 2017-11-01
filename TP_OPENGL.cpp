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

void drawCurve(Point* tabPointsOfCurve, long nbPoints);
void drawGrille(Point* grille, long nbU, long nbV);
Point* hermiteCubicCurve(const Point& p0,
                         const Point& p1,
                         const Vector& v0,
                         const Vector& v1,
                         long nbU);
Point* bezierCurveByBernstein(Point* tabControlPoint,
                              long nbControlPoint,
                              long nbU);
Point* BezierCurveByCasteljau(Point* tabControlPoint,
                              long nbControlPoint,
                              long nbU);
Point calculRecCasteljau(Point* tabControlPoint, long nbControlPoint, double u);
double bernstein(int i, int n, double u);
int polyNewton(int n, int i);
int f(int n);

Point** surfaceCylindrique(Point* courbeBezier,
                           Point* droite,
                           long nBezier,
                           long nbPointCourbe,
                           long nbPointU,
                           long nbPointV);
Point** surfaceReglee(Point* ptControlCourbeBezier1,
                      Point* ptControlCourbeBezier2,
                      long nbControlPoint1,
                      long nbControlPoint2,
                      long nbPointU,
                      long nbPointV);
Point* BezierSurfaceByCasteljau(Point** grilleControlPoint,
                                long nbControlePointU,
                                long nbU,
                                long nbControlePointV,
                                long nbV);
Point calculRecSurfaceCasteljau(Point** grilleControlPoint,
                                long nbControlePointU,
                                long nbControlePointV,
                                double u,
                                double v);
long nbMeridienParallele = 8;

void drawCylindre(Point** cylindre, long nbMeridien);
Point** facetCylindre(double rayon, double hauteur, long nbMeridien);
void drawCone(Point* Cone, long nbMeridien);
Point* facetCone(long nbMeridien);
void drawSphere(Point** Sphere, long nbMeridienParallele);
Point** facet_Sphere(long rayon);
void displayVoxel(Point centre, double length);
void StdDisplayVoxel(Point* tab);
Point ** createVoxel(Point* Tab, Point centre, double rayon);
void CreateSpherebyVoxel(Point centre, double length, double resolution);
void GenerateSphere(Point* TabVoxel, Point centre, double length, double resolution,  double currentR);
double Distance(Point A, Point B);
void GenerateCylinder(Point* TabVoxel, Point centre, Vector hauteur, double length, double resolution, double currentR);
void CreateCylinderbyVoxel(Point centre, Vector hauteur, double length, double resolution);

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
