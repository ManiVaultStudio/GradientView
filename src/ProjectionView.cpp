#include "ProjectionView.h"

#include "util/Exception.h"

#include <vector>
#include <algorithm>

#include <QSize>
#include <QPainter>

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

ProjectionView::ProjectionView() :
    _pointRenderer()
{
    //setContextMenuPolicy(Qt::CustomContextMenu);
    //setAcceptDrops(true);
    //setMouseTracking(true);
    //setFocusPolicy(Qt::ClickFocus);

    _pointRenderer.setPointScaling(Absolute);
    _pointRenderer.setPointSize(5);

    QSurfaceFormat surfaceFormat;

    surfaceFormat.setRenderableType(QSurfaceFormat::OpenGL);

#ifdef __APPLE__
    // Ask for an OpenGL 3.3 Core Context as the default
    surfaceFormat.setVersion(3, 3);
    surfaceFormat.setProfile(QSurfaceFormat::CoreProfile);
    surfaceFormat.setSwapBehavior(QSurfaceFormat::DoubleBuffer);
    //QSurfaceFormat::setDefaultFormat(defaultFormat);
#else
    // Ask for an OpenGL 4.3 Core Context as the default
    surfaceFormat.setVersion(4, 3);
    surfaceFormat.setProfile(QSurfaceFormat::CoreProfile);
    surfaceFormat.setSwapBehavior(QSurfaceFormat::DoubleBuffer);
#endif

#ifdef _DEBUG
    surfaceFormat.setOption(QSurfaceFormat::DebugContext);
#endif

    surfaceFormat.setSamples(16);

    setFormat(surfaceFormat);
}

ProjectionView::~ProjectionView()
{
    disconnect(QOpenGLWidget::context(), &QOpenGLContext::aboutToBeDestroyed, this, &ProjectionView::cleanup);
    cleanup();
}

bool ProjectionView::isInitialized()
{
    return _isInitialized;
}

void ProjectionView::setData(const std::vector<Vector2f>* points)
{
    auto dataBounds = getDataBounds(*points);

    dataBounds.ensureMinimumSize(1e-07f, 1e-07f);
    dataBounds.makeSquare();
    dataBounds.expand(0.1f);

    _dataBounds = dataBounds;
    qDebug() << _dataBounds.getLeft() << _dataBounds.getRight() << _dataBounds.getTop();

    // Pass bounds and data to renderer
    _pointRenderer.setBounds(_dataBounds);
    _pointRenderer.setData(*points);

    _pointRenderer.setSelectionOutlineColor(Vector3f(1, 0, 0));
    _pointRenderer.setAlpha(0.5f);
    //_pointRenderer.setPointScaling(PointScaling::Relative);

    update();
}

void ProjectionView::setScalars(const Eigen::Block<Eigen::MatrixXf, -1, 1, true>& scalars, int selectedPoint)
{
    float minV = *std::min_element(scalars.begin(), scalars.end());
    float maxV = *std::max_element(scalars.begin(), scalars.end());

    std::vector<Vector3f> colors(scalars.size());
    for (int i = 0; i < scalars.size(); i++)
    {
        float dimValue = scalars[i] / (maxV - minV);

        colors[i] = (i == selectedPoint) ? Vector3f(1, 0, 0) : Vector3f(1 - dimValue, 1 - dimValue, 1 - dimValue);
    }

    _pointRenderer.setColors(colors);
    _pointRenderer.setScalarEffect(None);

    update();
}

void ProjectionView::setColors(const std::vector<Vector3f>& colors)
{
    _pointRenderer.setColors(colors);
    _pointRenderer.setScalarEffect(None);

    update();
}

void ProjectionView::setProjectionName(QString name)
{
    _projectionName = name;
}

void ProjectionView::setPointSizeScalars(const std::vector<float>& pointSizeScalars)
{
    _pointRenderer.setSizeChannelScalars(pointSizeScalars);
    _pointRenderer.setPointSize(*std::max_element(pointSizeScalars.begin(), pointSizeScalars.end()));

    update();
}

void ProjectionView::setPointOpacityScalars(const std::vector<float>& pointOpacityScalars)
{
    _pointRenderer.setOpacityChannelScalars(pointOpacityScalars);

    update();
}

void ProjectionView::initializeGL()
{
    initializeOpenGLFunctions();

    connect(context(), &QOpenGLContext::aboutToBeDestroyed, this, &ProjectionView::cleanup);

    // Initialize renderers
    _pointRenderer.init();

    // Set a default color map for both renderers
    _pointRenderer.setScalarEffect(PointEffect::None);

    // OpenGL is initialized
    _isInitialized = true;

    emit initialized();
}

void ProjectionView::resizeGL(int w, int h)
{
    _windowSize.setWidth(w);
    _windowSize.setHeight(h);

    _pointRenderer.resize(QSize(w, h));

    // Set matrix for normalizing from pixel coordinates to [0, 1]
    //toNormalisedCoordinates = Matrix3f(1.0f / w, 0, 0, 1.0f / h, 0, 0);

    // Take the smallest dimensions in order to calculate the aspect ratio
    int size = w < h ? w : h;

    float wAspect = (float)w / size;
    float hAspect = (float)h / size;
    float wDiff = ((wAspect - 1) / 2.0);
    float hDiff = ((hAspect - 1) / 2.0);

    //toIsotropicCoordinates = Matrix3f(wAspect, 0, 0, hAspect, -wDiff, -hDiff);
}

void ProjectionView::paintGL()
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
            glClearColor(1, 1, 1, 1);
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

            // Reset the blending function
            glEnable(GL_BLEND);
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

            _pointRenderer.render();
        }
        QFont font = painter.font();
        font.setPointSize(font.pointSize() * 2);
        painter.setFont(font);
        painter.drawText(30, 30, _projectionName);

        painter.endNativePainting();

        painter.end();
    }
    catch (std::exception& e)
    {
        hdps::util::exceptionMessageBox("Rendering failed", e);
    }
    catch (...) {
        hdps::util::exceptionMessageBox("Rendering failed");
    }
}

void ProjectionView::cleanup()
{
    qDebug() << "Deleting projection view, performing clean up...";
    _isInitialized = false;

    makeCurrent();
    _pointRenderer.destroy();
}
