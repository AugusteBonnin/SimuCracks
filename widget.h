#ifndef WIDGET_H
#define WIDGET_H

#include <QWidget>

QT_BEGIN_NAMESPACE
namespace Ui { class Widget; }
QT_END_NAMESPACE

typedef struct Edge {
    double t0,t1;
    QPointF A,B,C,D;
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
    QPointF  computePoint(const Edge &e, double t);
    QPointF computeDerivative(const Edge &e, double t);
    int iteration;
    QVector<int> facesToSubdivide;

    double computeMinDist(const QPointF &p1, const Edge &e2, double &d2);
private slots:
    void subdivide();

};
#endif // WIDGET_H
