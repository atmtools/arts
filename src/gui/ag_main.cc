#include <QApplication>

#include "ag_main.h"
#include "mainwindow_qt.h"

int ag_main(int argc, char *argv[])
{
  // Q_INIT_RESOURCE(application);
  
  QApplication app(argc, argv);
  app.setOrganizationName("The ARTS Developers");
  app.setApplicationName("ARTS - The Atmospheric Radiative Transfer Simulator");
  MainWindow mainWin;
  mainWin.show();
  return app.exec();
}
