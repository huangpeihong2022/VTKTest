#include "QQuickVTKInteractiveWidget.h"
#include "QQuickVTKRenderItem.h"
#include "QQuickVTKRenderWindow.h"
#include "vtkActor.h"
#include "vtkAppendPolyData.h"
#include "vtkCamera.h"
#include "vtkClipPolyData.h"
#include "vtkCommand.h"
#include "vtkConeSource.h"
#include "vtkGenericOpenGLRenderWindow.h"
#include "vtkGlyph3D.h"
#include "vtkImplicitPlaneRepresentation.h"
#include "vtkImplicitPlaneWidget2.h"
#include "vtkInteractorEventRecorder.h"
#include "vtkNew.h"
#include "vtkPNGWriter.h"
#include "vtkPlane.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkProperty.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRenderer.h"
#include "vtkSphereSource.h"
#include "vtkTestUtilities.h"
#include "vtkTesting.h"
#include "vtkWindowToImageFilter.h"

#include <QApplication>
#include <QDebug>
#include <QQmlApplicationEngine>
#include <QQuickWindow>
#include <QTimer>
#include <QUrl>
#include <qdebug.h>
#include <qlist.h>
#include <qobject.h>
#include <qqml.h>
#include <vtkAngleRepresentation.h>
#include <vtkAngleRepresentation2D.h>
#include <vtkAngleWidget.h>
#include <vtkLeaderActor2D.h>
#include <vtkProperty2D.h>
#include <vtkTextProperty.h>

int main(int argc, char *argv[]) {
    QQuickVTKRenderWindow::setupGraphicsBackend();
    qmlRegisterType<QQuickVTKRenderWindow>("QQuickVTKRenderWindow", 1, 0, "QQuickVTKRenderWindow");
    qmlRegisterType<QQuickVTKRenderItem>("QQuickVTKRenderItem", 1, 0, "QQuickVTKRenderItem");

    QGuiApplication app(argc, argv);

    QQmlApplicationEngine engine;
    // engine.addImportPath("G:/opensource/VTK/master/debug/lib/qml/Debug");
    engine.load(QUrl("qrc:///qml/main.qml"));

    QObject *topLevel = engine.rootObjects().value(0);
    QQuickWindow *window = qobject_cast<QQuickWindow *>(topLevel);

    window->show();

    // Fetch the QQuick window using the standard object name set up in the constructor
    QQuickVTKRenderItem *qquickvtkItem = dynamic_cast<QQuickVTKRenderItem *>(topLevel->findChild<QObject *>("ConeView"));

    // Create a cone pipeline and add it to the view
    vtkNew<vtkActor> actor;
    vtkNew<vtkPolyDataMapper> mapper;
    vtkNew<vtkConeSource> cone;
    mapper->SetInputConnection(cone->GetOutputPort());
    actor->SetMapper(mapper);
    qquickvtkItem->renderer()->AddActor(actor);
    qquickvtkItem->renderer()->ResetCamera();
    qquickvtkItem->renderer()->SetBackground(0, 0, 0);
    qquickvtkItem->renderer()->SetBackgroundAlpha(1);
    qquickvtkItem->update();
    app.exec();
}