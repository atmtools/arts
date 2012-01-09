#ifndef AG_MAINWINDOW_H
#define AG_MAINWINDOW_H

#include <QMainWindow>

QT_BEGIN_NAMESPACE
class QAction;
class QMenu;
class QPlainTextEdit;
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
  Q_OBJECT
  
public:
  MainWindow();
  
protected:
  void closeEvent(QCloseEvent *ev);
  
  private slots:
  void about();
  
private:
  void createActions();
  void createMenus();
  void createToolBars();
  void createStatusBar();
  void readSettings();
  void writeSettings();
  
  QMenu *fileMenu;
  QMenu *helpMenu;
  QAction *exitAct;
  QAction *aboutAct;
  QAction *aboutQtAct;
};

#endif /* AG_MAINWINDOW_H */

