// import related modules
import QtQuick 2.15
import QtQuick.Controls 2.15
import QtQuick.Window 2.15

// import the VTK module
// import VTK 9.0
import QQuickVTKRenderWindow 1.0
import QQuickVTKRenderItem 1.0

// window containing the application
ApplicationWindow {
    // title of the application
    title: qsTr("VTK QtQuick App")
    visible: true
    width: 400
    height: 400
    objectName: "_appWindow"
    color: palette.window

    SystemPalette {
        id: palette
        colorGroup: SystemPalette.Active
    }

    // Instantiate the vtk render window
    QQuickVTKRenderWindow {
        id: vtkwindow
        width: 400
        height: 400
    }

    // add one or more vtk render items
    QQuickVTKRenderItem {
        objectName: "ConeView"
        x: 0
        y: 0
        width: 400
        height: 400
        // Provide the handle to the render window
        renderWindow: vtkwindow
    }
}