#pragma once

#include <QGLWidget>
#include <QMouseEvent>
#include <QKeyEvent>
#include "Camera.h"
#include <QVector3D>
#include <vector>
#include "LSystem.h"

class MainWindow;

class GLWidget3D : public QGLWidget {
public:
	MainWindow* mainWin;
	Camera camera;
	QPoint lastPos;
	lsystem::LSystem lsystem;

public:
	GLWidget3D(MainWindow *parent);
	void drawScene();

protected:
	void initializeGL();
	void resizeGL(int width, int height);
	void paintGL();    
	void mousePressEvent(QMouseEvent *e);
	void mouseMoveEvent(QMouseEvent *e);
	void mouseReleaseEvent(QMouseEvent *e);

};

