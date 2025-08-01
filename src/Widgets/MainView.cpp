#include "MainView.h"
#include "Application.h"

#include "util/PixelSelectionTool.h"
#include "util/Math.h"
#include "util/Exception.h"

#include <vector>

#include <QSize>
#include <QPainter>
#include <QDebug>
#include <QOpenGLFramebufferObject>
#include <QWindow>

#include <math.h>

namespace
{
    Bounds getDataBounds(const std::vector<Vector2f>& points)
    {
        Bounds bounds = Bounds::Max;

        for (const Vector2f& point : points)
        {
            bounds.setLeft(std::min(point.x, bounds.getLeft()));
            bounds.setRight(std::max(point.x, bounds.getRight()));
            bounds.setBottom(std::min(point.y, bounds.getBottom()));
            bounds.setTop(std::max(point.y, bounds.getTop()));
        }

        return bounds;
    }
}

MainView::MainView() :
    _densityRenderer(DensityRenderer::RenderMode::DENSITY),
    _backgroundColor(0, 0, 22, 255),
    _pointRenderer(this),
    _cellRenderer(this),
    //_pixelSelectionTool(this),
    _showRandomWalk(false),
    _showDirections(false),
    _pixelRatio(1.0)
{
    setContextMenuPolicy(Qt::CustomContextMenu);
    setAcceptDrops(true);
    setMouseTracking(true);
    setFocusPolicy(Qt::ClickFocus);

    _pointRenderer.setPointScaling(Relative);

    _pointRenderer.getNavigator().setZoomMarginScreen(25.f);
    _cellRenderer.getNavigator().setZoomMarginScreen(25.f);

    //// Configure pixel selection tool
    //_pixelSelectionTool.setEnabled(true);
    //_pixelSelectionTool.setMainColor(QColor(Qt::black));
    
    //QObject::connect(&_pixelSelectionTool, &PixelSelectionTool::shapeChanged, [this]() {
    //    if (isInitialized())
    //        update();
    //});

    QSurfaceFormat surfaceFormat;

    surfaceFormat.setRenderableType(QSurfaceFormat::OpenGL);
    surfaceFormat.setVersion(3, 3);
    surfaceFormat.setProfile(QSurfaceFormat::CoreProfile);
    surfaceFormat.setSwapBehavior(QSurfaceFormat::DoubleBuffer);
    surfaceFormat.setSamples(16);
    surfaceFormat.setStencilBufferSize(8);

#ifdef _DEBUG
    surfaceFormat.setOption(QSurfaceFormat::DebugContext);
#endif

    setFormat(surfaceFormat);

    // Call updatePixelRatio when the window is moved between hi and low dpi screens
    // e.g., from a laptop display to a projector
    // Wait with the connection until we are sure that the window is created
    connect(this, &MainView::created, this, [this](){
        [[maybe_unused]] auto windowID = this->window()->winId(); // This is needed to produce a valid windowHandle on some systems

        QWindow* winHandle = windowHandle();

        // On some systems we might need to use a different windowHandle
        if(!winHandle)
        {
            const QWidget* nativeParent = nativeParentWidget();
            winHandle = nativeParent->windowHandle();
        }

        if(winHandle == nullptr)
        {
            qDebug() << "MainView: Not connecting updatePixelRatio - could not get window handle";
            return;
        }

        QObject::connect(winHandle, &QWindow::screenChanged, this, &MainView::updatePixelRatio, Qt::UniqueConnection);
    });
}

bool MainView::isInitialized()
{
    return _isInitialized;
}

MainView::RenderMode MainView::getRenderMode() const
{
    return _renderMode;
}

void MainView::setRenderMode(const RenderMode& renderMode)
{
    if (renderMode == _renderMode)
        return;

    _renderMode = renderMode;

    emit renderModeChanged(_renderMode);

    switch (_renderMode)
    {
        case MainView::SCATTERPLOT:
            break;
        
        case MainView::DENSITY:
            computeDensity();
            break;

        case MainView::LANDSCAPE:
            computeDensity();
            break;

        default:
            break;
    }

    update();
}

//MainView::ColoringMode MainView::getColoringMode() const
//{
//    return _coloringMode;
//}

//void MainView::setColoringMode(const ColoringMode& coloringMode)
//{
//    if (coloringMode == _coloringMode)
//        return;
//
//    _coloringMode = coloringMode;
//
//    emit coloringModeChanged(_coloringMode);
//}

//PixelSelectionTool& MainView::getPixelSelectionTool()
//{
//    return _pixelSelectionTool;
//}

void MainView::computeDensity()
{
    emit densityComputationStarted();

    _densityRenderer.computeDensity();

    emit densityComputationEnded();

    update();
}

Matrix3f createProjectionMatrix(Bounds bounds)
{
    Matrix3f m;
    m.setIdentity();
    m[0] = 2 / bounds.getWidth();
    m[4] = 2 / bounds.getHeight();
    m[6] = -((bounds.getRight() + bounds.getLeft()) / bounds.getWidth());
    m[7] = -((bounds.getTop() + bounds.getBottom()) / bounds.getHeight());
    return m;
}

// Positions need to be passed as a pointer as we need to store them locally in order
// to be able to find the subset of data that's part of a selection. If passed
// by reference then we can upload the data to the GPU, but not store it in the widget.
void MainView::setData(const std::vector<Vector2f>* points)
{
    _dataBounds = getDataBounds(*points);

    const auto dataBoundsRect = QRectF(QPointF(_dataBounds.getLeft(), _dataBounds.getBottom()), QSizeF(_dataBounds.getWidth(), _dataBounds.getHeight()));

    // Pass bounds and data to renderers
    _pointRenderer.setDataBounds(dataBoundsRect);
    _densityRenderer.setDataBounds(dataBoundsRect);
    _cellRenderer.setDataBounds(dataBoundsRect);

    _pointRenderer.getNavigator().resetView(true);
    _cellRenderer.getNavigator().resetView(true);

    _pointRenderer.setData(*points);
    _densityRenderer.setData(points);
    //_cellRenderer.setData(*points); qDebug() << "Cell renderer done";

    switch (_renderMode)
    {
        case MainView::SCATTERPLOT:
            break;
        
        case MainView::DENSITY:
        case MainView::LANDSCAPE:
        {
            _densityRenderer.computeDensity();
            break;
        }

        default:
            break;
    }

    _pointRenderer.setSelectionOutlineColor(Vector3f(1, 0, 0));

    update();
}

QColor MainView::getBackgroundColor()
{
    return _backgroundColor;
}

void MainView::setBackgroundColor(QColor color)
{
    _backgroundColor = color;

    update();
}

void MainView::setHighlights(const std::vector<char>& highlights, const std::int32_t& numSelectedPoints)
{
    _pointRenderer.setHighlights(highlights, numSelectedPoints);

    update();
}

void MainView::setScalars(const std::vector<float>& scalars)
{
    _pointRenderer.setColorChannelScalars(scalars);
    _cellRenderer.setColorChannelScalars(scalars);
    
    update();
}

void MainView::setColors(const std::vector<Vector3f>& colors)
{
    _pointRenderer.setColors(colors);
    //_pointRenderer.setScalarEffect(None);

    update();
}

void MainView::setPointSize(float pointSize)
{
    _pointRenderer.setPointSize(pointSize);

    update();
}

void MainView::setPointSizeScalars(const std::vector<float>& pointSizeScalars)
{
    _pointRenderer.setSizeChannelScalars(pointSizeScalars);
    _pointRenderer.setPointSize(*std::max_element(pointSizeScalars.begin(), pointSizeScalars.end()));

    update();
}

void MainView::setPointOpacityScalars(const std::vector<float>& pointOpacityScalars)
{
    _pointRenderer.setOpacityChannelScalars(pointOpacityScalars);

    update();
}

void MainView::setPointScaling(mv::gui::PointScaling scalingMode)
{
    _pointRenderer.setPointScaling(scalingMode);
}

void MainView::setScalarEffect(PointEffect effect)
{
    _pointRenderer.setScalarEffect(effect);
    _cellRenderer.setScalarEffect(ScalarEffect::Color);

    update();
}

void MainView::setSigma(const float sigma)
{
    _densityRenderer.setSigma(sigma);

    update();
}

mv::Vector3f MainView::getColorMapRange() const
{
    switch (_renderMode) {
        case SCATTERPLOT:
            return _pointRenderer.getColorMapRange();

        case LANDSCAPE:
            return _densityRenderer.getColorMapRange();

        default:
            break;
    }
    
    return Vector3f();
}

void MainView::setColorMapRange(const float& min, const float& max)
{
    switch (_renderMode) {
        case SCATTERPLOT:
        {
            _pointRenderer.setColorMapRange(min, max);
            break;
        }

        case LANDSCAPE:
        {
            _densityRenderer.setColorMapRange(min, max);
            break;
        }
        case CELL:
        {
            _cellRenderer.setColorMapRange(min, max);
            break;
        }

        default:
            break;
    }

    update();
}

void MainView::createScreenshot(std::int32_t width, std::int32_t height, const QString& fileName, const QColor& backgroundColor)
{
    // Exit if the viewer is not initialized
    if (!_isInitialized)
        return;

    // Exit prematurely if the file name is invalid
    if (fileName.isEmpty())
        return;

    makeCurrent();

    try {

        // Use custom FBO format
        QOpenGLFramebufferObjectFormat fboFormat;
        
        fboFormat.setTextureTarget(GL_TEXTURE_2D);
        fboFormat.setInternalTextureFormat(GL_RGB);

        QOpenGLFramebufferObject fbo(width, height, fboFormat);

        // Bind the FBO and render into it when successfully bound
        if (fbo.bind()) {

            // Clear the widget to the background color
            glClearColor(backgroundColor.redF(), backgroundColor.greenF(), backgroundColor.blueF(), backgroundColor.alphaF());
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

            // Reset the blending function
            glEnable(GL_BLEND);
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

            // Resize OpenGL to intended screenshot size
            resizeGL(width, height);

            switch (_renderMode)
            {
                case SCATTERPLOT:
                {
                    _pointRenderer.setPointScaling(Relative);
                    _pointRenderer.render();
                    _pointRenderer.setPointScaling(Absolute);

                    break;
                }

                case DENSITY:
                case LANDSCAPE:
                    _densityRenderer.setRenderMode(_renderMode == DENSITY ? DensityRenderer::DENSITY : DensityRenderer::LANDSCAPE);
                    _densityRenderer.render();
                    break;
            }

            // Save FBO image to disk
            //fbo.toImage(false, QImage::Format_RGB32).convertToFormat(QImage::Format_RGB32).save(fileName);
            //fbo.toImage(false, QImage::Format_ARGB32).save(fileName);

            QImage fboImage(fbo.toImage());
            QImage image(fboImage.constBits(), fboImage.width(), fboImage.height(), QImage::Format_ARGB32);

            image.save(fileName);

            // Resize OpenGL back to original OpenGL widget size
            resizeGL(this->width(), this->height());

            fbo.release();
        }
    }
    catch (std::exception& e)
    {
        exceptionMessageBox("Rendering failed", e);
    }
    catch (...) {
        exceptionMessageBox("Rendering failed");
    }
}

void MainView::setRandomWalks(const std::vector<std::vector<Vector2f>>& randomWalks)
{
    _randomWalks = randomWalks;
}

void MainView::setDirections(const std::vector<Vector2f>& directions)
{
    _directions = directions;
}

void MainView::initializeGL()
{
    initializeOpenGLFunctions();

#ifdef SCATTER_PLOT_WIDGET_VERBOSE
    qDebug() << "Initializing scatterplot widget with context: " << context();

    std::string versionString = std::string((const char*) glGetString(GL_VERSION));

    qDebug() << versionString.c_str();
#endif

    connect(context(), &QOpenGLContext::aboutToBeDestroyed, this, &MainView::cleanup);

    // Initialize renderers
    _pointRenderer.init();
    _densityRenderer.init();
    _cellRenderer.init();

    // Set a default color map for both renderers
    _pointRenderer.setScalarEffect(PointEffect::Color);

    // OpenGL is initialized
    _isInitialized = true;

    // Initialize the point and density renderer with a color map
    setColorMap(_colorMapImage);

    emit initialized();
}

void MainView::resizeGL(int w, int h)
{
    // we need this here as we do not have the screen yet to get the actual devicePixelRatio when the view is created
    //_pixelRatio = devicePixelRatio();

    // Pixelration tells us how many pixels map to a point
    // That is needed as macOS calculates in points and we do in pixels
    // On macOS high dpi displays pixel ration is 2
    //w *= _pixelRatio;
    //h *= _pixelRatio;

    _windowSize.setWidth(w);
    _windowSize.setHeight(h);

    _pointRenderer.resize(QSize(w, h));
    _densityRenderer.resize(QSize(w, h));
    _cellRenderer.resize(QSize(w, h));

    // Set matrix for normalizing from pixel coordinates to [0, 1]
    toNormalisedCoordinates = Matrix3f(1.0f / w, 0, 0, 1.0f / h, 0, 0);

    // Take the smallest dimensions in order to calculate the aspect ratio
    int size = w < h ? w : h;

    float wAspect = (float)w / size;
    float hAspect = (float)h / size;
    float wDiff = ((wAspect - 1) / 2.0);
    float hDiff = ((hAspect - 1) / 2.0);

    toIsotropicCoordinates = Matrix3f(wAspect, 0, 0, hAspect, -wDiff, -hDiff);
}

void MainView::paintGL()
{
    try {
        QPainter painter;

        // Begin mixed OpenGL/native painting
        if (!painter.begin(this))
            throw std::runtime_error("Unable to begin painting");

        // Draw layers with OpenGL
        painter.beginNativePainting();
        {
            // Bind the framebuffer belonging to the widget
            //glBindFramebuffer(GL_FRAMEBUFFER, defaultFramebufferObject());

            // Clear the widget to the background color
            glClearColor(_backgroundColor.redF(), _backgroundColor.greenF(), _backgroundColor.blueF(), _backgroundColor.alphaF());
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

            // Reset the blending function
            glEnable(GL_BLEND);
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

            switch (_renderMode)
            {
                case SCATTERPLOT:
                    _pointRenderer.render();
                    break;

                case DENSITY:
                case LANDSCAPE:
                    _densityRenderer.setRenderMode(_renderMode == DENSITY ? DensityRenderer::DENSITY : DensityRenderer::LANDSCAPE);
                    _densityRenderer.render();
                    break;
                case CELL:
                    _cellRenderer.render();
                    break;
            }

        }
        painter.endNativePainting();
        
        QFont font = painter.font();
        font.setPointSize(24);
        painter.setFont(font);
        painter.setPen(QPen(Qt::white));
        painter.drawText(14, 40, _projectionName);

        font.setPointSize(16);
        painter.setFont(font);
        painter.drawText(14, 70, _coloredBy);

        painter.drawText(14, 100, "Section " + _clusterName);

        // Draw random walks
        // [Bounds to -1, 1]
        Matrix3f orthoM = createProjectionMatrix(_dataBounds);
        
        int size = std::min(width(), height());
        Matrix3f toScreen;
        toScreen.setIdentity();
        toScreen[0] = size / 2;
        toScreen[4] = size / 2;
        toScreen[6] = width() / 2;
        toScreen[7] = height() / 2;

        Matrix3f invM;
        invM.setIdentity();
        invM[0] = 1;
        invM[4] = -1;

        const auto screenPoint  = _pointRenderer.getWorldPositionToScreenPoint(QVector3D(_currentPoint.x, _currentPoint.y, 0.f));
        const auto cp           = Vector2f(screenPoint.x(), screenPoint.y());

        if (_showFilterCircles)
        {
            // Render peak filter circles
            const auto screenPointInnerRadius = _pointRenderer.getWorldPositionToScreenPoint(QVector3D(_currentPoint.x + _radii.x, _currentPoint.y, 0.f));
            const auto screenPointOuterRadius = _pointRenderer.getWorldPositionToScreenPoint(QVector3D(_currentPoint.x + _radii.y, _currentPoint.y, 0.f));

            const auto inner_r = (screenPointInnerRadius - screenPoint).x();
            const auto outer_r = (screenPointOuterRadius - screenPoint).x();

            painter.setPen(QPen(QColor(255, 0, 0, 255)));
            painter.drawEllipse(QPointF(cp.x, cp.y), inner_r, inner_r);
            painter.drawEllipse(QPointF(cp.x, cp.y), outer_r, outer_r);
        }

        // Render selection dot
        painter.setBrush(Qt::red);
        painter.drawEllipse(QPointF(cp.x, cp.y), 5, 5);
        painter.setBrush(Qt::BrushStyle::NoBrush);

        if (_showRandomWalk)
        {
            painter.setPen(QPen(QColor(0, 0, 0, 255)));
            for (int i = 0; i < _randomWalks.size(); i++)
            {
                if (_randomWalks[i].size() < 2) continue;
                for (int j = 0; j < _randomWalks[i].size() - 1; j++)
                {
                    const Vector2f& p1 = toScreen * invM * orthoM * _randomWalks[i][j];
                    const Vector2f& p2 = toScreen * invM * orthoM * _randomWalks[i][j + 1];

                    painter.drawLine(p1.x, p1.y, p2.x, p2.y);
                }
            }
        }

        if (_showDirections)
        {
            for (int i = 0; i < _directions.size(); i += 2)
            {
                Vector2f p = _directions[i];
                Vector2f d = _directions[i+1];
                //const Vector2f& p1 = toScreen * invM * orthoM * p;
                const Vector2f& pt = toScreen * invM * orthoM * p;

                painter.drawLine(pt.x - d.x * 5, pt.y + d.y * 5, pt.x + d.x * 5, pt.y - d.y * 5);
            }
        }

        // Draw metadata
        int x = 10;
        int y = 160;
        painter.setFont(font);
        painter.setPen(QPen(Qt::white));
        for (QString s : _metadataList)
        {
            painter.drawText(x, y, s);
            y += 20;
        }

        //// Draw the pixel selection tool overlays if the pixel selection tool is enabled
        //if (_pixelSelectionTool.isEnabled()) {
        //    painter.drawPixmap(rect(), _pixelSelectionTool.getAreaPixmap());
        //    painter.drawPixmap(rect(), _pixelSelectionTool.getShapePixmap());
        //}
        
        painter.end();
    }
    catch (std::exception& e)
    {
        exceptionMessageBox("Rendering failed", e);
    }
    catch (...) {
        exceptionMessageBox("Rendering failed");
    }
}

void MainView::cleanup()
{
    qDebug() << "Deleting scatterplot widget, performing clean up...";
    _isInitialized = false;

    makeCurrent();
    _pointRenderer.destroy();
    _densityRenderer.destroy();
    _cellRenderer.destroy();
}

void MainView::setColorMap(const QImage& colorMapImage)
{
    _colorMapImage = colorMapImage;

    // Do not update color maps of the renderers when OpenGL is not initialized
    if (!_isInitialized)
        return;

    // Activate OpenGL context
    makeCurrent();

    // Apply color maps to renderers
    _pointRenderer.setColormap(_colorMapImage);
    _densityRenderer.setColormap(_colorMapImage);
    _cellRenderer.setColormap(_colorMapImage);

    // Render
    update();
}

void MainView::updatePixelRatio()
{
    float pixelRatio = devicePixelRatio();

#ifdef SCATTER_PLOT_WIDGET_VERBOSE
    qDebug() << "Window moved to screen " << window()->screen() << ".";
    qDebug() << "Pixelratio before was " << _pixelRatio << ". New pixelratio is: " << pixelRatio << ".";
#endif // SCATTER_PLOT_WIDGET_VERBOSE

    // we only update if the ratio actually changed
    if (_pixelRatio != pixelRatio)
    {
        _pixelRatio = pixelRatio;
        resizeGL(width(), height());
        update();
    }
}

MainView::~MainView()
{
    disconnect(QOpenGLWidget::context(), &QOpenGLContext::aboutToBeDestroyed, this, &MainView::cleanup);
    cleanup();
}
