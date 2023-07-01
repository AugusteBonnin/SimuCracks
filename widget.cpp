#include "widget.h"
#include "ui_widget.h"

#include <QPainter>
#include <QTimer>
#include <QVector2D>

Widget::Widget(QWidget *parent)
    : QWidget(parent)
    , ui(new Ui::Widget)
{
    srand(time(NULL));
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
    e1.A = QPointF(0,0);
    e1.B = QPointF(0,128);
    e1.C = QPointF(0,384);
    e1.D = QPointF(0,511);
    e2.A = QPointF(0,511);
    e2.B = QPointF(128,511);
    e2.C = QPointF(384,511);
    e2.D = QPointF(511,511);
    e3.A = QPointF(511,511);
    e3.B = QPointF(511,384);
    e3.C = QPointF(511,128);
    e3.D = QPointF(511,0);
    e4.A = QPointF(511,0);
    e4.B = QPointF(384,0);
    e4.C = QPointF(128,0);
    e4.D = QPointF(0,0);

    edges << e1;
    edges << e2;
    edges << e3;
    edges << e4;

    Face f;
    f.edge=0;
    faces << f;

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

QPointF Widget::computePoint(const Edge & e,double t)
{
    QPointF res((1-t)*(1-t)*(1-t)*e.A);
    res += 3*t*(1-t)*(1-t)*e.B;
    res += 3*t*t*(1-t)*e.C;
    res += t*t*t*e.D;
    return res;
}

QPointF Widget::computeDerivative(const Edge & e,double t)
{
    QPointF res(((-2+2*t)+(4*t-3*t*t))*e.A);
    res += (1-(-4*t)+(3*t*t))*e.B;
    res += (t*(2-3*t))*e.B;
    res += 3*t*t*e.B;
    return res;
}

double Widget::computeAnisotropicDist(const QPointF & O,const QVector2D & tangent,const QVector2D & normal,const QPointF & P)
{
    //const double factor=.025;
    double factor = (iteration<4)?1:.025;;
    QVector2D  OP(P-O);
    double dpx = QVector2D::dotProduct(tangent,OP);
    double dpy = factor*QVector2D::dotProduct(normal,OP);
    return dpx*dpx+dpy*dpy;
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

void Widget::computeEdgesPerFaceHistogram()
{
    uint16_t count[5];
    for (int i = 0 ; i < 5 ; ++i)
        count[i]=0;
    for (int i = 0 ; i < faces.count() ; ++i)
    {
        uint16_t total = 0 ;
        Face f = faces[i];
        int edge = f.edge;
        do {
            ++total ;
            edge = edges[edge].next;
        } while (edge != f .edge);
        switch (total)
        {
        case 3 :
            count[0]++;
            break;
        case 4 :
            count[1]++;
            break;
        case 5:
            count[2]++;
            break ;
        case 6 :
            count[3]++;
            break;
        default:
            count[4]++;
            break;

        }
    }

    QString text;
    text += QString("Triangles : %1 (%2%%)\n").arg(count[0]).arg(count[0]*100.0f/faces.count());
    text += QString("Quads : %1 (%2%%)\n").arg(count[1]).arg(count[1]*100.0f/faces.count());
    text += QString("Pentas : %1 (%2%%)\n").arg(count[2]).arg(count[2]*100.0f/faces.count());
    text += QString("Hexas : %1 (%2%%)\n").arg(count[3]).arg(count[3]*100.0f/faces.count());
    text += QString("Heptas & + : %1 (%2%%)\n").arg(count[4]).arg(count[4]*100.0f/faces.count());
    text += QString("Total : %1 (%2%%)\n").arg(faces.count()).arg(faces.count()*100.0f/faces.count());
    ui->label_2->setText(text);
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
        int emax2_idx ;//= edges[edges[emax1_idx].next].next;
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
        QPointF p1prime = computePoint(emax1,d1+(emax1.t1-emax1.t0)*.01);
        int emax1prime_idx=edges.count();
        Edge emax1prime(emax1);
        emax1.t1=d1;
        emax1prime.t0=d1;
        edges << emax1prime;


        double min = INFINITY ;
        double d2,min_t;

        QVector2D tangent(p1prime-p1);
        tangent.normalize();
        QVector2D normal(-tangent.y(),tangent.x());

        do {
            Edge & e = edges[e_idx];
            if (e_idx!=emax1_idx)
            {
                min_t = e.t0+(e.t1-e.t0)*gaussianValue();
                QPointF p = computePoint(e,min_t);
                double dist = computeAnisotropicDist(p1,tangent,normal,p);
                if (dist < min)
                {
                    min = dist;
                    d2 = min_t ;
                    emax2_idx = e_idx ;
                }
            }
            e_idx = e.next ;
        } while (e_idx!=f.edge);


        Edge & emax2 = edges[emax2_idx];
        /*double d2;
        double dist = computeMinDist(p1,emax2,d2);
        d2 += (.3*(rand()/(double)RAND_MAX)-.15)*(emax2.t1-emax2.t0);
        if (d2<emax2.t0)
            d2 = emax2.t0;
        else if (d2>emax2.t1)
            d2 = emax2.t1;
            */
        QPointF p2 = computePoint(emax2,d2);
        QPointF p2prime = computePoint(emax2,d2+(emax2.t1-emax2.t0)*.01);
        QLineF l(p1,p2);
        double d = l.length();
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
        QLineF l1(p1,p1prime);
        l1 = l1.normalVector();
        l1.setLength(d*.6);
        QLineF l2(p2,p2prime);
        l2 = l2.normalVector();
        l2.setLength(d/4);
        e12.A = p1 ;
        e12.B = l1.p2();
        e12.C = l2.p2();
        e12.D = p2 ;
        edges << e12;
        emax1.next=e12_idx;
        edges[emax1_idx]=emax1;

        int e21_idx=edges.count();
        Edge e21;
        e21.age=iteration;
        e21.next=emax1prime_idx;
        e21.t0 = 0 ;
        e21.t1 = 1 ;
        e21.D = p1 ;
        e21.C = l1.p2();
        e21.B = l2.p2();
        e21.A = p2 ;
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

    computeEdgesPerFaceHistogram();

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

    QImage widget_image(this->size(),QImage::Format_ARGB32);
    QPainter painter2(&widget_image);
    this->render(&painter2);
    widget_image.save(QString("%1.jpg").arg(iteration));

    if (facesToSubdivide.count())
    //if (iteration<6)
        QTimer::singleShot(2000,this,SLOT(subdivide()));
}
