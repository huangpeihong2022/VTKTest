#include <iostream>
#include <sstream>

#include <vtkActor.h>
#include <vtkAxesActor.h>
#include <vtkSmartPointer.h>
#include <vtkVolume.h>
#include <vtkNew.h>
#include <vtkPolyDataMapper.h>
#include <vtkCylinderSource.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkDataSetMapper.h>
#include <vtkCompositePolyDataMapper.h>
#include <vtkTubeFilter.h>
#include <vtkLineSource.h>
#include <vtkConeSource.h>
#include <vtkCubeSource.h>
#include <vtkSphereSource.h>
#include <vtkLineSource.h>
#include <vtkParametricSpline.h>
#include <vtkParametricFunctionSource.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkProperty.h>
#include <vtkCamera.h>
#include <vtkLight.h>
#include <vtkDICOMImageReader.h>
#include <vtkXYPlotActor.h>
#include <vtkProperty2D.h>
#include <vtkTextProperty.h>
#include <vtkImageExtractComponents.h>
#include <vtkImageAccumulate.h>
#include <vtkImageData.h>
#include <vtkGPUVolumeRayCastMapper.h>
#include <vtkVolumeProperty.h>
#include <vtkPiecewiseFunction.h>
#include <vtkColorTransferFunction.h>
#include <vtkBarChartActor.h>
#include <vtkFieldData.h>
#include <VtkLegendBoxActor.h>
#include <QtCore/qmath.h> 
#include <QDebug> 

void TestCyLinder()
{
    // You can only create cylinders along the Y-axis of the world coordinate system
    vtkSmartPointer<vtkCylinderSource> cylinderSource =vtkSmartPointer<vtkCylinderSource>::New();
    cylinderSource->SetHeight(10.0);
    cylinderSource->SetCenter(0.0, 0.0, 0.0);
    cylinderSource->SetRadius(2.0);
    cylinderSource->SetResolution(50);

    //Create a directional cylinder
	vtkSmartPointer<vtkLineSource> lineSource = vtkSmartPointer<vtkLineSource>::New();
    //Set the direction of the long axis of the cylinder（x,y,z）
	lineSource->SetPoint1(5.0, 0.0, 0.0);
	lineSource->SetPoint2(-5.0, 0.0, 0.0);

    vtkSmartPointer<vtkTubeFilter> tubeFilter = vtkSmartPointer<vtkTubeFilter>::New();
	tubeFilter->SetInputConnection(lineSource->GetOutputPort());
	tubeFilter->SetRadius(2.0);
	tubeFilter->SetNumberOfSides(50);
	tubeFilter->CappingOn();

    //feasible VTKMapper
    vtkSmartPointer<vtkPolyDataMapper> cylinderSourcemapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    vtkSmartPointer<vtkPolyDataMapper> lineSourcemapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    //vtkNew<vtkDataSetMapper> mapper;
    //vtkNew<vtkCompositePolyDataMapper> mapper;
    cylinderSourcemapper->SetInputConnection(cylinderSource->GetOutputPort());
    lineSourcemapper->SetInputConnection(tubeFilter->GetOutputPort());

    vtkSmartPointer<vtkActor> cylinderSourceactor = vtkSmartPointer<vtkActor>::New();
    vtkSmartPointer<vtkActor> lineSourceactor = vtkSmartPointer<vtkActor>::New();
    cylinderSourceactor->SetMapper(cylinderSourcemapper);
    lineSourceactor->SetMapper(lineSourcemapper);

    vtkSmartPointer<vtkRenderer> cylinderSourcerender = vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkRenderer> lineSourcerender = vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->AddRenderer(cylinderSourcerender);
    renderWindow->AddRenderer(lineSourcerender);
    renderWindow->SetSize(800,680);
    vtkSmartPointer<vtkRenderWindowInteractor> Interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    Interactor->SetRenderWindow(renderWindow);
    Interactor->Initialize();

    cylinderSourcerender->AddActor(cylinderSourceactor);
    cylinderSourcerender->ResetCamera();
    cylinderSourcerender->SetBackground(0,0,0);
	cylinderSourcerender->SetViewport(0.0, 0.0, 0.5, 1.0);

    lineSourcerender->AddActor(lineSourceactor);
    lineSourcerender->ResetCamera();
    lineSourcerender->SetBackground(0,0,0);
    lineSourcerender->SetViewport(0.5, 0.0, 1.0, 1.0 );

    renderWindow->Render();
    Interactor->Start();
}

void TestCone()
{
    vtkSmartPointer<vtkConeSource> cone = vtkSmartPointer<vtkConeSource>::New();
    cone->SetHeight(2.5);
    cone->SetRadius(1.0);
    cone->SetResolution(100);
    cone->SetDirection(1,0,0);

    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    //vtkSmartPointer<vtkDataSetMapper> mapper = vtkSmartPointer<vtkDataSetMapper>::New();
    //vtkSmartPointer<vtkCompositePolyDataMapper> mapper = vtkSmartPointer<vtkCompositePolyDataMapper>::New();
    mapper->SetInputConnection(cone->GetOutputPort());

    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);

    vtkSmartPointer<vtkLight> light = vtkSmartPointer<vtkLight>::New();
    light->SetColor(1.0,0.0,0.0);

    vtkSmartPointer<vtkRenderer> render = vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkRenderWindow> renderwindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderwindow->AddRenderer(render);
    renderwindow->SetSize(800,600);

    vtkSmartPointer<vtkRenderWindowInteractor> interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    interactor->SetRenderWindow(renderwindow);
    interactor->Initialize();

    render->AddActor(actor);
    //render->AddLight(light);
    render->ResetCamera();
    render->SetBackground(0,0,0);

    renderwindow->Render();
    
    //The camera is rotated 360 degrees
	// for (int i=0; i<360; i++)
	// {
	// 	renderwindow->Render();
	// 	//Rotate the camera 1 degree in the scene
	// 	render->GetActiveCamera()->Azimuth(1);
	// }
    interactor->Start();
}

void TestCube()
{
    vtkSmartPointer<vtkCubeSource> cube = vtkSmartPointer<vtkCubeSource>::New();
    cube->SetXLength(2);
    cube->SetYLength(2);
    cube->SetZLength(2);

    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(cube->GetOutputPort());

    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);

    vtkSmartPointer<vtkRenderer> render = vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkRenderWindow> renderwindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderwindow->AddRenderer(render);
    renderwindow->SetSize(800,600);

    vtkSmartPointer<vtkRenderWindowInteractor> interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    interactor->SetRenderWindow(renderwindow);
    interactor->Initialize();

    render->AddActor(actor);
    render->ResetCamera();
    render->SetBackground(0,0,0);

    renderwindow->Render();
    interactor->Start();
}

void TestSphere()
{
    vtkSmartPointer<vtkSphereSource> sphere = vtkSmartPointer<vtkSphereSource>::New();
    sphere->SetCenter(0,0,0);
    sphere->SetRadius(4.0);
    sphere->SetPhiResolution(50);
    sphere->SetThetaResolution(50);

    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(sphere->GetOutputPort());

    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);

    vtkSmartPointer<vtkRenderer> render = vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkRenderWindow> renderwindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderwindow->AddRenderer(render);
    renderwindow->SetSize(800,600);

    vtkSmartPointer<vtkRenderWindowInteractor> interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    interactor->SetRenderWindow(renderwindow);
    interactor->Initialize();

    render->AddActor(actor);
    render->ResetCamera();
    render->SetBackground(0,0,0);

    renderwindow->Render();
    interactor->Start();
}

void TestLine()
{
    vtkSmartPointer<vtkLineSource> sphere = vtkSmartPointer<vtkLineSource>::New();
    sphere->SetPoint1(0,0,0);
    sphere->SetPoint2(5,0,0);

    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(sphere->GetOutputPort());

    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);

    vtkSmartPointer<vtkRenderer> render = vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkRenderWindow> renderwindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderwindow->AddRenderer(render);
    renderwindow->SetSize(800,600);

    vtkSmartPointer<vtkRenderWindowInteractor> interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    interactor->SetRenderWindow(renderwindow);
    interactor->Initialize();

    render->AddActor(actor);
    render->ResetCamera();
    render->SetBackground(0,0,0);

    renderwindow->Render();
    interactor->Start();
}

void TestCurve()
{
    double p0[3] = { 1.0, 0.0, 0.0 };
    double p1[3] = { 0.0, 1.0, 0.0 };
    double p2[3] = { 0.0, 0.0, 1.0 };
    double p3[3] = { 1.0, 2.0, 3.0 };

    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->InsertNextPoint(p0);
    points->InsertNextPoint(p1);
    points->InsertNextPoint(p2);
    points->InsertNextPoint(p3);

    vtkSmartPointer<vtkParametricSpline> spline = vtkSmartPointer<vtkParametricSpline>::New();
    spline->SetPoints(points);
        
    vtkSmartPointer<vtkParametricFunctionSource> functionSource = vtkSmartPointer<vtkParametricFunctionSource>::New();
    functionSource->SetParametricFunction(spline);
    functionSource->Update();

    vtkSmartPointer<vtkPolyDataMapper> functionMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    functionMapper->SetInputConnection(functionSource->GetOutputPort());

    vtkSmartPointer<vtkActor> functionActor = vtkSmartPointer<vtkActor>::New();
    functionActor->SetMapper(functionMapper);

    vtkSmartPointer<vtkPolyData> resultPolydata = vtkSmartPointer<vtkPolyData>::New();
    resultPolydata->SetPoints(points);

    vtkSmartPointer<vtkVertexGlyphFilter> resultGlyphFilter = vtkSmartPointer<vtkVertexGlyphFilter>::New();
    resultGlyphFilter->AddInputData(resultPolydata);
    resultGlyphFilter->Update();

    vtkSmartPointer<vtkPolyDataMapper> resultMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    resultMapper->SetInputConnection(resultGlyphFilter->GetOutputPort());
    vtkSmartPointer<vtkActor> resultActor = vtkSmartPointer<vtkActor>::New();
    resultActor->SetMapper(resultMapper);
    resultActor->GetProperty()->SetPointSize(5);//定义点的尺寸大小，这样点才能在画布上显示出来
    resultActor->GetProperty()->SetColor(1, 0.0, 0.0);
        
    vtkSmartPointer<vtkRenderer> render = vtkSmartPointer<vtkRenderer>::New();
    render->AddActor(functionActor);
    render->AddActor(resultActor);
    render->SetBackground(0.1, 0.2, 0.4);

    vtkSmartPointer<vtkRenderWindow> renderwindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderwindow->AddRenderer(render);
    renderwindow->SetSize(300, 300);

    vtkSmartPointer<vtkRenderWindowInteractor> interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    interactor->SetRenderWindow(renderwindow);
    interactor->Initialize();
    interactor->Start();
}

void TestCSYS()
{
    vtkSmartPointer<vtkAxesActor> actor = vtkSmartPointer<vtkAxesActor>::New();
    actor->SetPosition(0, 0, 0);
    actor->SetTotalLength(2, 2, 2);
    actor->SetShaftType(0);
    actor->SetAxisLabels(2);
    actor->SetCylinderRadius(0.02);

    vtkSmartPointer<vtkRenderer> render = vtkSmartPointer<vtkRenderer>::New();
    render->AddActor(actor);
    render->SetBackground(0, 0, 0);

    vtkSmartPointer<vtkRenderWindow> renderwindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderwindow->AddRenderer(render);
    renderwindow->SetSize(300, 300);

    vtkSmartPointer<vtkRenderWindowInteractor> interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    interactor->SetRenderWindow(renderwindow);
    interactor->Initialize();
    interactor->Start();
}

void TestGrid2D()
{
    vtkSmartPointer<vtkPoints> _linePoints = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> _lineCellArray = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPolyData> _linePolyData  = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkDataSetMapper> _linePolyMapper = vtkSmartPointer<vtkDataSetMapper>::New();
    double rangeX[2] = {0,5};
    double rangeY[2] = {0,5};
    double gridSizeX = 1, gridSizeY = 1;

    // 平行Y方向的直线
    for(double gridX = rangeX[0]; gridX < rangeX[1] + (gridSizeX / 2.0); gridX += gridSizeX)
    {
        double lineStart[3] = {gridX, rangeY[0], 0.0};
        double lineEnd[3]   = {gridX, rangeY[1], 0.0};
 
        vtkIdType pointIdStart = _linePoints->InsertNextPoint(lineStart);
        vtkIdType pointIdEnd   = _linePoints->InsertNextPoint(lineEnd);
 
        vtkIdType singleLineCell[2] = {pointIdStart, pointIdEnd};
        _lineCellArray->InsertNextCell(2, singleLineCell);
    }
 
    // 平行X方向的直线
    for(double gridY = rangeY[0]; gridY < rangeY[1] + (gridSizeY / 2.0); gridY += gridSizeY)
    {
        double lineStart[3] = {rangeX[0], gridY, 0.0};
        double lineEnd[3]   = {rangeX[1], gridY, 0.0};
 
        vtkIdType pointIdStart = _linePoints->InsertNextPoint(lineStart);
        vtkIdType pointIdEnd   = _linePoints->InsertNextPoint(lineEnd);
 
        vtkIdType singleLineCell[2] = {pointIdStart, pointIdEnd};
        _lineCellArray->InsertNextCell(2, singleLineCell);
    }
 
    _linePolyData->SetLines(_lineCellArray);
    _linePolyData->SetPoints(_linePoints);
 
    _linePolyMapper->SetInputData(_linePolyData);

    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(_linePolyMapper);

    vtkSmartPointer<vtkRenderer> render = vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkRenderWindow> renderwindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderwindow->AddRenderer(render);
    renderwindow->SetSize(800,600);

    vtkSmartPointer<vtkRenderWindowInteractor> interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    interactor->SetRenderWindow(renderwindow);
    interactor->Initialize();

    render->AddActor(actor);
    render->ResetCamera();
    render->SetBackground(0,0,0);

    renderwindow->Render();
    interactor->Start();
}

void TestGrid3D()
{
    vtkSmartPointer<vtkPoints> _linePoints = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> _lineCellArray = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPolyData> _linePolyData  = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkDataSetMapper> _linePolyMapper = vtkSmartPointer<vtkDataSetMapper>::New();
    double rangeX[2] = {0,5};
    double rangeY[2] = {0,5};
    double rangeZ[2] = {0,5};
    double gridSizeX = 1, gridSizeY = 1, gridSizeZ = 1;
    // Z方向切片
    for(double gridZ = rangeZ[0]; gridZ < rangeZ[1] + (gridSizeZ / 2.0); gridZ += gridSizeZ)
    {
        // 平行Y方向的直线
        for(double gridX = rangeX[0]; gridX < rangeX[1] + (gridSizeX / 2.0); gridX += gridSizeX)
        {
            double lineStart[3] = {gridX, rangeY[0], gridZ};
            double lineEnd[3]   = {gridX, rangeY[1], gridZ};
 
            vtkIdType pointIdStart = _linePoints->InsertNextPoint(lineStart);
            vtkIdType pointIdEnd   = _linePoints->InsertNextPoint(lineEnd);
 
            vtkIdType singleLineCell[2] = {pointIdStart, pointIdEnd};
            _lineCellArray->InsertNextCell(2, singleLineCell);
        }
 
        // 平行X方向的直线
        for(double gridY = rangeY[0]; gridY < rangeY[1] + (gridSizeY / 2.0); gridY += gridSizeY)
        {
            double lineStart[3] = {rangeX[0], gridY, gridZ};
            double lineEnd[3]   = {rangeX[1], gridY, gridZ};
 
            vtkIdType pointIdStart = _linePoints->InsertNextPoint(lineStart);
            vtkIdType pointIdEnd   = _linePoints->InsertNextPoint(lineEnd);
 
            vtkIdType singleLineCell[2] = {pointIdStart, pointIdEnd};
            _lineCellArray->InsertNextCell(2, singleLineCell);
        }
    }
 
    for(double gridX = rangeX[0]; gridX < rangeX[1] + (gridSizeX / 2.0); gridX += gridSizeX)
    {
        for(double gridY = rangeY[0]; gridY < rangeY[1] + (gridSizeY / 2.0); gridY += gridSizeY)
        {
            double lineStart[3] = {gridX, gridY, rangeZ[0]};
            double lineEnd[3]   = {gridX, gridY, rangeZ[1]};
 
            vtkIdType pointIdStart = _linePoints->InsertNextPoint(lineStart);
            vtkIdType pointIdEnd   = _linePoints->InsertNextPoint(lineEnd);
 
            vtkIdType singleLineCell[2] = {pointIdStart, pointIdEnd};
            _lineCellArray->InsertNextCell(2, singleLineCell);
        }
    }
 
    _linePolyData->SetLines(_lineCellArray);
    _linePolyData->SetPoints(_linePoints);
 
    _linePolyMapper->SetInputData(_linePolyData);

    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(_linePolyMapper);

    vtkSmartPointer<vtkRenderer> render = vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkRenderWindow> renderwindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderwindow->AddRenderer(render);
    renderwindow->SetSize(800,600);

    vtkSmartPointer<vtkRenderWindowInteractor> interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    interactor->SetRenderWindow(renderwindow);
    interactor->Initialize();

    render->AddActor(actor);
    render->ResetCamera();
    render->SetBackground(0,0,0);

    renderwindow->Render();
    interactor->Start();

}

void TestHistogram()
{
    vtkSmartPointer<vtkDICOMImageReader>reader=vtkSmartPointer<vtkDICOMImageReader>::New();
    reader->SetDirectoryName("D:/spine_dcm/dcm_bsk_lumbar_case0016_ori");
    reader->SetDataByteOrderToLittleEndian();
    reader->Update();

    vtkSmartPointer<vtkGPUVolumeRayCastMapper> volumeMapper = vtkSmartPointer<vtkGPUVolumeRayCastMapper>::New();
    volumeMapper->SetInputData(reader->GetOutput());
    volumeMapper->SetSampleDistance(0.1);
    volumeMapper->SetAutoAdjustSampleDistances(0);
    volumeMapper->SetBlendModeToComposite();

    vtkSmartPointer<vtkVolumeProperty> volumeProperty = vtkSmartPointer<vtkVolumeProperty>::New();
    volumeProperty->SetInterpolationTypeToLinear();
    volumeProperty->ShadeOn(); //打开或者关闭阴影测试.
    volumeProperty->SetAmbient(0.4);
    volumeProperty->SetDiffuse(0.6);  //漫反射.
    volumeProperty->SetSpecular(0.2); //镜面反射.

    //设置不透明度.
    vtkSmartPointer<vtkPiecewiseFunction> compositeOpacity = vtkSmartPointer<vtkPiecewiseFunction>::New();
    compositeOpacity->AddPoint(-3024, 0, 0.5, 0.0);
    compositeOpacity->AddPoint(-16, 0, .49, .61);
    compositeOpacity->AddPoint(641, .72, .5, 0.0);
    compositeOpacity->AddPoint(3071, .71, 0.5, 0.0);
    volumeProperty->SetScalarOpacity(compositeOpacity); //设置不透明度传输函数.

    //设置梯度不透明属性.
    vtkSmartPointer<vtkPiecewiseFunction> volumeGradientOpacity = vtkSmartPointer<vtkPiecewiseFunction>::New();
    volumeGradientOpacity->AddPoint(10, 0.0);
    volumeGradientOpacity->AddPoint(90, 0.5);
    volumeGradientOpacity->AddPoint(100, 1.0);
    volumeProperty->SetGradientOpacity(volumeGradientOpacity); //设置梯度不透明度效果对比.

    //设置颜色属性.
    vtkSmartPointer<vtkColorTransferFunction> color = vtkSmartPointer<vtkColorTransferFunction>::New();
    color->AddRGBPoint(-3024, 0, 0, 0, 0.5, 0.0);
    color->AddRGBPoint(-16, 0.73, 0.25, 0.30, 0.49, .61);
    color->AddRGBPoint(641, .90, .82, .56, .5, 0.0);
    color->AddRGBPoint(3071, 1, 1, 1, .5, 0.0);
    volumeProperty->SetColor(color);

    vtkSmartPointer<vtkVolume> volume = vtkSmartPointer<vtkVolume>::New();
    volume->SetMapper(volumeMapper);
    volume->SetProperty(volumeProperty);

    vtkSmartPointer<vtkRenderer> Render = vtkSmartPointer<vtkRenderer>::New();
    Render->AddVolume(volume);
    Render->ResetCamera();
    Render->SetBackground(0, 0, 0);
    Render->SetViewport(0, 0, 0.5, 1);

    int numComponents = reader->GetOutput()->GetNumberOfScalarComponents();

    vtkSmartPointer<vtkXYPlotActor> plot = vtkSmartPointer<vtkXYPlotActor>::New();
    plot->ExchangeAxesOff();
    plot->SetLabelFormat( "%g" );
    plot->SetXTitle( "CT" );
    plot->SetYTitle( "Red" );
    plot->SetXValuesToValue();
    plot->GetProperty()->SetColor(0.0, 0.0, 0.0);
    plot->GetAxisLabelTextProperty()->SetColor(0.0, 0.0, 0.0);
    plot->GetAxisTitleTextProperty()->SetColor(0.0, 0.0, 0.0);
    
    double colors[3][3] = {
        { 1, 0, 0 },    
        { 0, 1, 0 },
        { 0, 0, 1 }
    };

    const char* labels[3] = { "Red", "Green", "Blue" };
    
    int xmax = 0;   //最大横坐标
    int ymax = 0;   //最大纵坐标
    
    for( int i = 0; i < numComponents; ++i )
    {
        //彩色图像不能直接计算直方图，因此需要先提取每个通道图像
        vtkSmartPointer<vtkImageExtractComponents> extract =
                vtkSmartPointer<vtkImageExtractComponents>::New();
        extract->SetInputConnection( reader->GetOutputPort() );
        extract->SetComponents( i );
        extract->Update();
        double range[2];
        extract->GetOutput()->GetScalarRange( range );

        //直方图间隔的个数为：最大灰度值减去最小灰度值，再减1
        int extent = static_cast<int>(range[1])-static_cast<int>(range[0])-1;
        
        vtkSmartPointer<vtkImageAccumulate> histogram =
                vtkSmartPointer<vtkImageAccumulate>::New();
        histogram->SetInputConnection( extract->GetOutputPort() );
        histogram->SetComponentExtent( 0,extent, 0,0, 0,0);   //直方图间隔的个数
        histogram->SetComponentOrigin( range[0],0,0 );        //灰度起点为图像的最小灰度值
        histogram->SetComponentSpacing( 1,0,0);              //直方图的间隔取(1, 0, 0)，即每个灰度计算统计一个频率
        histogram->SetIgnoreZero( 1 );                        //在统计直方图时，像素值为0的像素不进行统计
        histogram->Update();
        
        if( range[1] > xmax )
        {
            xmax = range[1];
        }
        if( histogram->GetOutput()->GetScalarRange()[1] > ymax )
        {
            ymax = histogram->GetOutput()->GetScalarRange()[1];
        }
        plot->AddDataSetInput(histogram->GetOutput());
        plot->SetPlotColor(i,colors[i]);
        plot->SetPlotLabel(i,labels[i]); 

        plot->LegendOn();
    }
    //设置X轴和Y轴的数据范围
    plot->SetXRange(-3024, xmax);
    plot->SetYRange(0, 600000);

	vtkSmartPointer<vtkRenderer> Histogramrenderer = vtkSmartPointer<vtkRenderer>::New();
	Histogramrenderer->AddActor(plot);
	Histogramrenderer->SetBackground(1.0, 1.0, 1.0);
    Histogramrenderer->SetViewport(0.55,0.55,0.99,0.99);
    
    vtkSmartPointer<vtkRenderWindow> renderWindow =vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->AddRenderer(Render);
    renderWindow->AddRenderer(Histogramrenderer);
    renderWindow->SetSize(1500, 1000);
    renderWindow->Render();
    renderWindow->SetWindowName("ImageAccumulateExample2");
    
    vtkSmartPointer<vtkRenderWindowInteractor> interactor =vtkSmartPointer<vtkRenderWindowInteractor>::New();
    interactor->SetRenderWindow(renderWindow);
    interactor->Initialize();
    interactor->Start();
}

int main(int argc, char *argv[]) {
    //TestCyLinder();
    //TestCone();
    //TestCube();
    //TestSphere();
    //TestLine();
    //TestCurve();
    //TestCSYS();
    //TestGrid2D();
    //TestGrid3D();
    TestHistogram();

    return 0;
}