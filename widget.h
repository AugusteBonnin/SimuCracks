#ifndef WIDGET_H
#define WIDGET_H

#include <QWidget>

QT_BEGIN_NAMESPACE
namespace Ui { class Widget; }
QT_END_NAMESPACE

typedef struct Edge {
    double t0,t1;
    double ax,bx,cx,dx,ay,by,cy,dy;
    int next;
    int age;

}Edge ;

typedef struct {
    int edge;
} Face ;

class Widget : public QWidget
{
    Q_OBJECT

public:
    Widget(QWidget *parent = nullptr);
    ~Widget();

private:
    Ui::Widget *ui;
    QVector<Edge> edges;
    QVector<Face> faces;
    double gaussianValue();
    void computeEdge(Edge &e, double x0, double y0, double dx0, double dy0, double x1, double y1, double dx1, double dy1);
    QPointF  computePoint(const Edge &e, double t);
    QPointF computeDerivative(const Edge &e, double t);
    int iteration;
    QVector<int> facesToSubdivide;

    double computeMinDist(const QPointF &p1, const Edge &e2, double &d2);
private slots:
    void subdivide();

};
#endif // WIDGET_H
