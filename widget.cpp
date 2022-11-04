#include "widget.h"
#include "ui_widget.h"

#include <QPainter>
#include <QTimer>
#include <QVector2D>

Widget::Widget(QWidget *parent)
    : QWidget(parent)
    , ui(new Ui::Widget)
{
    iteration = 0 ;
    ui->setupUi(this);
    Edge e1,e2,e3,e4;
    e1.next=1;
    e2.next=2;
    e3.next=3;
    e4.next=0;
    e1.age = e2.age = e3.age = e4.age = 0 ;
    e1.t0 = e2.t0 = e3.t0 = e4.t0 = 0 ;
    e1.t1 = e2.t1 = e3.t1 = e4.t1 = 1 ;
    computeEdge(e1,0,0,0,511,0,511,0,511);
    computeEdge(e2,0,511,511,0,511,511,511,0);
    computeEdge(e3,511,511,0,-511,511,0,0,-511);
    computeEdge(e4,511,0,-511,0,0,0,-511,0);
    edges << e1;
    edges << e2;
    edges << e3;
    edges << e4;
    Face f;
    f.edge=0;
    faces << f;
    srand(time(NULL));
    facesToSubdivide << 0 ;
    subdivide();

}

Widget::~Widget()
{
    delete ui;
}

double Widget::gaussianValue()
{
    double v;
    do {
        v = rand()/(double)RAND_MAX;
    } while ((1.0f/((1/24.0)*sqrt(2*M_PI)))*exp(-(v-.5)*(v-.5)/(2*(1/24.0)*(1/24.0)))<rand()/(double)RAND_MAX);
    return v ;
}

void Widget::computeEdge(Edge & e,double x0,double y0,double dx0,double dy0,double x1,double y1,double dx1,double dy1)
{
    e.dx=x0;
    e.cx=dx0;
    e.ax=dx1+dx0-2*x1+2*x0;
    e.bx=-2*dx0-dx1+3*x1-3*x0;
    e.dy=y0;
    e.cy=dy0;
    e.ay=dy1+dy0-2*y1+2*y0;
    e.by=-2*dy0-dy1+3*y1-3*y0;
}

QPointF Widget::computePoint(const Edge & e,double t)
{
    QPointF res(e.dx+t*(e.cx+t*(e.bx+t*e.ax)),e.dy+t*(e.cy+t*(e.by+t*e.ay)));
    return res;
}

QPointF Widget::computeDerivative(const Edge & e,double t)
{
    QPointF res(e.cx+t*(2*e.bx+t*3*e.ax),e.cy+t*(2*e.by+t*3*e.ay));
    return res;
}

double Widget::computeMinDist(const QPointF & p1, const Edge & e2, double & d2)
{
double res = INFINITY;
for (double t = e2.t0 ; t <= e2.t1 ; t+= (e2.t1-e2.t0)*.01)
{
    QPointF p2 = computePoint(e2,t);
    QLineF l(p1,p2);
    double d = l.length();
    if (res>d)
    {
        res=d;
        d2=t;
    }
}
return res;
}

void Widget::subdivide()
{
    ++iteration;
    QVector<int> newFacesToSubdivide;
    for (int i = 0 ; i< facesToSubdivide.count();++i)
    {
        Face & f = faces[facesToSubdivide[i]];
        int e_idx = f.edge;
        double max = 0 ;
        int emax1_idx;
        do {
            Edge & e = edges[e_idx];
            QLineF l(computePoint(e,e.t0),computePoint(e,e.t1));
            double d = l.length()*(.95+.1*(rand()/(double)RAND_MAX));
            if (d>max)
            {
                max = d ;
                emax1_idx = e_idx ;
            }
            e_idx = e.next ;
        } while (e_idx!=f.edge);
        if (max < 80)
            continue;
        int emax2_idx = edges[edges[emax1_idx].next].next;
        /*
        max = 0 ;
        do {
            Edge & e = edges[e_idx];
            QLineF l(computePoint(e,e.t0),computePoint(e,e.t1));
            double d = l.length()*(.95+.1*(rand()/(double)RAND_MAX));
            if ((d>max)&&(e_idx!=emax1_idx))
            {
                max = d ;
                emax2_idx = e_idx ;
            }
            e_idx = e.next ;
        } while (e_idx!=f.edge);
        */
        /*
         * ______________
         * |             |
         * |             |
         * |emax1        |emax2
         * |             |
         * |_____________|
         *
         * _______________+p4
         * |emax1         |emax2prime
         * |              |
         * +p1__e12_______+p2
         * |______e21_____|
         * |emax1prime    |emax2
         * +______________|
         * p3
          */
        Edge & emax1 = edges[emax1_idx];
        double d1 = emax1.t0+(emax1.t1-emax1.t0)*gaussianValue();
        QPointF p1 = computePoint(emax1,d1);
        int emax1prime_idx=edges.count();
        Edge emax1prime(emax1);
        emax1.t1=d1;
        emax1prime.t0=d1;
        edges << emax1prime;

        /*
        double min = INFINITY ;
        double d2,min_t;

        do {
            Edge & e = edges[e_idx];
            if (e_idx!=emax1_idx)
            {
                double dist = computeMinDist(p1,e,min_t);
                if (dist < min)
                {
                    min = dist;
                    d2 = min_t ;
                    emax2_idx = e_idx ;
                }
            }
            e_idx = e.next ;
        } while (e_idx!=f.edge);
*/

        Edge & emax2 = edges[emax2_idx];
        double d2;
        double dist = computeMinDist(p1,emax2,d2);
        d2 += (.2*(rand()/(double)RAND_MAX)-.1)*(emax2.t1-emax2.t0);
        if (d2<emax2.t0)
            d2 = emax2.t0;
        else if (d2>emax2.t1)
            d2 = emax2.t1;
        QPointF p2 = computePoint(emax2,d2);
        QLineF l(p1,p2);
        double d = l.length();
        //QPointF dp1 = computeDerivative(emax1,d1)/pow(1.4,(512-d)/128);
        //QPointF dp2 = computeDerivative(emax2,d2)/pow(1.4,(512-d)/128);
        //QPointF dp1 = computeDerivative(emax1,d1)/(4/(d/80));
        //QPointF dp2 = computeDerivative(emax2,d2)/(4/(d/80));
        QPointF dp1 = computeDerivative(emax1,d1)/((d>160)?1:4);
        QPointF dp2 = computeDerivative(emax2,d2)/((d>160)?1:4);
        int emax2prime_idx=edges.count();
        Edge emax2prime(emax2) ;
        emax2.t1=d2;
        emax2prime.t0 = d2 ;
        edges << emax2prime;

        int e12_idx=edges.count();
        Edge e12;
        e12.age=iteration;
        e12.next=emax2prime_idx;
        e12.t0=0 ;
        e12.t1=1;
        computeEdge(e12,p1.x(),p1.y(),dp1.y(),-dp1.x(),p2.x(),p2.y(),-dp2.y(),dp2.x());
        edges << e12;
        emax1.next=e12_idx;
        edges[emax1_idx]=emax1;

        int e21_idx=edges.count();
        Edge e21;
        e21.age=iteration;
        e21.next=emax1prime_idx;
        e21.t0 = 0 ;
        e21.t1 = 1 ;
        computeEdge(e21,p2.x(),p2.y(),dp2.y(),-dp2.x(),p1.x(),p1.y(),-dp1.y(),dp1.x());
        edges << e21;
        emax2.next=e21_idx;
        edges[emax2_idx]=emax2;

        f.edge=emax1_idx;
        faces[facesToSubdivide[i]]=f;
        newFacesToSubdivide << facesToSubdivide[i];
        int f2_idx = faces.count();
        Face f2;
        f2.edge=emax2_idx;
        faces << f2 ;
        newFacesToSubdivide << f2_idx ;

    }
    facesToSubdivide = newFacesToSubdivide ;

    QImage image(512,512,QImage::Format_ARGB32);
    image.fill(Qt::white);
    QPainter painter(&image);
    for (int i = 0 ; i < edges.count() ; ++i)
    {
        Edge & e = edges[i];
        QPen pen;
        pen.setWidth(iteration+1-e.age);
        painter.setPen(pen);
        QPolygonF poly;
        for (double j = e.t0 ; j <= e.t1 ; j += .01)
            poly << computePoint(e,j);
        painter.drawPolyline(poly);
    }
    ui->label->setPixmap(QPixmap::fromImage(image));
    image.save(QString("%1.jpg").arg(iteration));

    if (facesToSubdivide.count())
    //if (iteration<6)
        QTimer::singleShot(2000,this,SLOT(subdivide()));
}
