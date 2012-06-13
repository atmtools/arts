#include <QtGui>

#include "mainwindow_qt.h"

MainWindow::MainWindow()
{
  QWidget *widget = new QWidget;
  QHBoxLayout *l = new QHBoxLayout;
  l->addWidget(new QLabel("Hello World"));
  
  widget->setLayout(l);
  widget->show();
  setCentralWidget(widget);
  
  createActions();
  createMenus();
  createStatusBar();
  
  readSettings();
  
  setUnifiedTitleAndToolBarOnMac(true);
}

void MainWindow::closeEvent(QCloseEvent *ev)
{
  writeSettings();
  ev->accept();
}

void MainWindow::about()
{
  QMessageBox::about(this, tr("About Application"),
                     tr("The <b>Application</b> example demonstrates how to "
                        "write modern GUI applications using Qt, with a menu bar, "
                        "toolbars, and a status bar."));
}

void MainWindow::createActions()
{
  exitAct = new QAction(tr("E&xit"), this);
  exitAct->setShortcuts(QKeySequence::Quit);
  exitAct->setStatusTip(tr("Exit the application"));
  connect(exitAct, SIGNAL(triggered()), this, SLOT(close()));
  
  aboutAct = new QAction(tr("&About"), this);
  aboutAct->setStatusTip(tr("Show the application's About box"));
  connect(aboutAct, SIGNAL(triggered()), this, SLOT(about()));
  
  aboutQtAct = new QAction(tr("About &Qt"), this);
  aboutQtAct->setStatusTip(tr("Show the Qt library's About box"));
  connect(aboutQtAct, SIGNAL(triggered()), qApp, SLOT(aboutQt()));
}

void MainWindow::createMenus()
{
  fileMenu = menuBar()->addMenu(tr("&File"));
  fileMenu->addAction(exitAct);

  menuBar()->addSeparator();

  helpMenu = menuBar()->addMenu(tr("&Help"));
  helpMenu->addAction(aboutAct);
  helpMenu->addAction(aboutQtAct);
}

void MainWindow::createStatusBar()
{
  statusBar()->showMessage(tr("Ready"));
}

void MainWindow::readSettings()
{
  QSettings settings("The ARTS Developers", "ARTS");
  QPoint p = settings.value("pos", QPoint(200, 200)).toPoint();
  QSize s = settings.value("size", QSize(400, 400)).toSize();
  resize(s);
  move(p);
}

void MainWindow::writeSettings()
{
  QSettings settings("Trolltech", "Application Example");
  settings.setValue("pos", pos());
  settings.setValue("size", size());
}

