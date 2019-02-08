#include <iostream>
#include <fstream>
#include <cassert>
#include <iomanip>
#include <cmath>
#include <stack>
#include <queue>
#include <algorithm>
#include "bitmap_image.hpp"
using namespace std;

#define pi (2*acos(0.0))
#define epsilon (1.0e-6)

class homogeneous_point
{
public:
    double x, y, z, w;

    // set the three coordinates, set w to 1
    homogeneous_point(double x, double y, double z)
    {
        this->x = x;
        this->y = y;
        this->z = z;
        this->w = 1;
    }

    /*
    default constructor. does nothing. allows declarations like below:
        matrix m;
    therefore, usage is dangerous
    */
    homogeneous_point() {
    }

    // constructs a homogeneous point with given coordinates. forces w to be 1.0
    // if w is zero, raises error
    homogeneous_point(double x, double y, double z, double w)
    {
        assert (w != 0);
        this->x = x/w;
        this->y = y/w;
        this->z = z/w;
        this->w = 1;
    }

    // adds two points. returns a point forcing w to be 1.0
    homogeneous_point operator+ (const homogeneous_point& point)
    {
        double x = this->x + point.x;
        double y = this->y + point.y;
        double z = this->z + point.z;
        double w = this->w + point.w;
        homogeneous_point p(x, y, z, w);
        return p;
    }

    // subtracts one point from another. returns a point forcing w to be 1.0
    homogeneous_point operator- (const homogeneous_point& point)
    {
        double x = this->x - point.x;
        double y = this->y - point.y;
        double z = this->z - point.z;
        double w = this->w - point.w;
        homogeneous_point p(x, y, z, w);
    }

    //----------------------
    /*
    homogeneous_point operator= (const homogeneous_point& point)
    {
        double x = point.x;
        double y = point.y;
        double z = point.z;
        double w = point.w;
        homogeneous_point p(x, y, z, w);
    }
    */

    // Print the coordinates of a point. exists for testing purpose.
    void print()
    {
        cout << "Point: " << endl;
        cout << x << " " << y << " " << z << " " << w << endl;
    }

};


class Vector
{
public:
    double x, y, z;

    // constructs a vector with given components
    Vector(double x, double y, double z) {
        this->x = x;
        this->y = y;
        this->z = z;
    }

    // keeps the direction same. recalculates the vector to be unit.
    void normalize()
    {
        double r = sqrt(x*x + y*y + z*z);
        x = x / r;
        y = y / r;
        z = z / r;
    }

    // add two vectors
    Vector operator+(const Vector& v)
    {
        Vector v1(x+v.x, y+v.y, z+v.z);
        return v1;
    }

    // subtract one vector from another
    Vector operator-(const Vector& v)
    {
        Vector v1(x-v.x, y-v.y, z-v.z);
        return v1;
    }

    // scale a vector with a given coefficient
    Vector operator* (double m)
    {
        Vector v(x*m, y*m, z*m);
        return v;
    }

    // get the dot product of two vectors
    static double dot(Vector a, Vector b)
    {
        return a.x*b.x + a.y*b.y + a.z*b.z;
    }

    // get the cross product of two vectors
    static Vector cross(Vector a, Vector b)
    {
        Vector v(a.y*b.z - a.z*b.y, b.x*a.z - b.z*a.x, a.x*b.y - a.y*b.x);
        return v;
    }

    // print a vector. only for testing purposes.
    void print ()
    {
        cout << "Vector" << endl;
        cout << x << " " << y << " " << z << endl;
    }
};


/*
The matrices are forced to be 4x4. This is because in this assignment, we will deal with points in triangles.
Maximum # of points that we will deal with at once is 3. And all the standard matrices are 4x4 (i.e. scale, translation, rotation etc.)
*/
class matrix
{
public:
    double values[4][4];
    int num_rows, num_cols;

    // only set the number of rows and cols
    matrix(int rows, int cols)
    {
        assert (rows <= 4 && cols <= 4);
        num_rows = rows;
        num_cols = cols;
    }

    // prepare an nxn square matrix
    matrix(int n)
    {
        assert (n <= 4);
        num_rows = num_cols = n;
    }

    // prepare and return an identity matrix of size nxn
    static matrix make_identity(int n)
    {
        assert (n <= 4);
        matrix m(n);
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                if (i == j)
                    m.values[i][j] = 1;
                else
                    m.values[i][j] = 0;
            }
        }
        return m;
    }

    // print the matrix. exists for testing purposes
    void print()
    {
        cout << "Matrix:" << endl;
        for (int i = 0; i < num_rows; i++)
        {
            for (int j = 0; j < num_cols; j++)
            {
                cout << values[i][j] << "\t";
            }
            cout << endl;
        }
    }

    // add the two matrices. Raise error if dimension mismatches
    matrix operator+ (const matrix& m)
    {
        assert (this->num_rows == m.num_rows);
        assert (this->num_cols == m.num_cols);

        matrix m1(num_rows, num_cols);
        for (int i = 0; i < num_rows; i++)
        {
            for (int j = 0; j < num_cols; j++)
            {
                m1.values[i][j] = values[i][j] + m.values[i][j];
            }
        }
        return m1;
    }

    // subtract a matrix from another. raise error if dimension mismatches
    matrix operator- (const matrix& m)
    {
        assert (this->num_rows == m.num_rows);
        assert (this->num_cols == m.num_cols);

        matrix m1(num_rows, num_cols);
        for (int i = 0; i < num_rows; i++)
        {
            for (int j = 0; j < num_cols; j++)
            {
                m1.values[i][j] = values[i][j] - m.values[i][j];
            }
        }
        return m1;
    }

    // multiply two matrices. allows statements like m1 = m2 * m3; raises error is dimension mismatches
    matrix operator* (const matrix& m)
    {
        assert (this->num_cols == m.num_rows);
        matrix m1(this->num_rows, m.num_cols);

        for (int i = 0; i < m1.num_rows; i++) {
            for (int j = 0; j < m1.num_cols; j++) {
                double val = 0;
                for (int k = 0; k < this->num_cols; k++) {
                    val += this->values[i][k] * m.values[k][j];
                }
                m1.values[i][j] = val;
            }
        }
        return m1;
    }

    // multiply a matrix with a constant
    matrix operator* (double m)
    {
        matrix m1(this->num_rows, this->num_cols);
        for (int i = 0; i < num_rows; i++) {
            for (int j = 0; j < num_cols; j++) {
                m1.values[i][j] = m * this->values[i][j];
            }
        }
        return m1;
    }

    // multiply a 4x4 matrix with a homogeneous point and return the resulting point.
    // usage: homogeneous_point p = m * p1;
    // here, m is a 4x4 matrix, intended to be the transformation matrix
    // p1 is the point on which the transformation is being made
    // p is the resulting homogeneous point
    homogeneous_point operator* (const homogeneous_point& p)
    {
        assert (this->num_rows == this->num_cols && this->num_rows == 4);

        matrix m(4, 1);
        m.values[0][0] = p.x;
        m.values[1][0] = p.y;
        m.values[2][0] = p.z;
        m.values[3][0] = p.w;

        matrix m1 = (*this)*m;
        homogeneous_point p1(m1.values[0][0], m1.values[1][0], m1.values[2][0], m1.values[3][0]);
        return p1;
    }

    // return the transpose of a matrix
    matrix transpose()
    {
        matrix m(num_cols, num_rows);
        for (int i = 0; i < num_rows; i++) {
            for (int j = 0; j < num_cols; j++) {
                m.values[j][i] = values[i][j];
            }
        }
        return m;
    }

    //------------------------
    void setIdx(int a, int b, double x){
        assert (a < 4 && b < 4);
        values[a][b]=x;
    }

    static matrix getScaleMatrix(double sx, double sy, double sz){ 
        matrix m=matrix::make_identity(4);
        m.setIdx(0,0,sx); 
        m.setIdx(1,1,sy);
        m.setIdx(2,2,sz);
        return m;
    }

    static matrix getTranslationMatrix(double dx, double dy, double dz){ 
        matrix m=matrix::make_identity(4);
        m.setIdx(0,3,dx); 
        m.setIdx(1,3,dy);
        m.setIdx(2,3,dz);
        return m;
    }

    static matrix getRotationMatrix(double theta,double kx, double ky, double kz){ 
        matrix m=matrix::make_identity(4);
        double c=cos(theta*pi/180);
        double s=sin(theta*pi/180);
        double norm=sqrt(kx*kx+ky*ky+kz*kz);
        kx=kx/norm;
        ky=ky/norm;
        kz=kz/norm;
        m.setIdx(0,0,kx*kx*(1-c)+c);
        m.setIdx(0,1,kx*ky*(1-c)-kz*s);
        m.setIdx(0,2,kx*kz*(1-c)+ky*s);
        m.setIdx(1,0,kx*ky*(1-c)+kz*s);
        m.setIdx(1,1,ky*ky*(1-c)+c);
        m.setIdx(1,2,ky*kz*(1-c)-kx*s);
        m.setIdx(2,0,kx*kz*(1-c)-ky*s);
        m.setIdx(2,1,ky*kz*(1-c)+kx*s);
        m.setIdx(2,2,kz*kz*(1-c)+c);
        return m;
    }

};

/*
A simple class to hold the color components, r, g, b of a certain shade.
*/
class color {
public:
    double r, g, b;
    color(double r, double g, double b) {
        this->r = r;
        this->g = g;
        this->b = b;
    }
    color() {
    }

    void printColor(){
        cout << "("<<r<<","<<g<<","<<b<<")"<<endl;
    } 
};


//-------Triangle-------------

class Triangle{
public:
    color c;
    homogeneous_point p1,p2,p3;
    bool clipped;

    Triangle(double ax,double ay,double az,double bx,double by,double bz,
            double cx,double cy,double cz,double r,double g,double b):p1(ax,ay,az),p2(bx,by,bz),
            p3(cx,cy,cz),c(r,g,b){
        clipped=false;
    }
    /*
    Triangle(homogeneous_point a,homogeneous_point b,homogeneous_point c,color clr):
        p1(a.x,a.y,a.z),p2(b.x,b.y,b.z),p3(c.x,c.y,c.z),c(clr.r,clr.g,clr.b){
        clipped=false;
    }*/
};


bool comparePoint(homogeneous_point h1, homogeneous_point h2){
    return h1.y<h2.y;
}


double eye_x, eye_y, eye_z;
double look_x, look_y, look_z;
double up_x, up_y, up_z;
double fov_x, fov_y, aspectRatio, near, far;
color backgroud;
int screen_x, screen_y;

int triangleSide=3;

vector<Triangle> triangles;

int getYIdx(double y){
    int idx=screen_y-1-int(0.5*(y+1)*screen_y);
    return idx;
}

int getXidx(double x){
    int idx=int(0.5*(x+1)*screen_x);
    return idx;
}

void scan_convert() {
    ifstream stage3;
    stage3.open("stage3.txt");

    color** pixels = new color*[screen_x];
    double** zs = new double*[screen_x];
    for (int i = 0; i < screen_x; i++) {
        pixels[i] = new color [screen_y];
        for (int j = 0; j < screen_y; j++) {
            pixels[i][j] = backgroud;
        }
        zs[i] = new double [screen_y];
        for (int j = 0; j < screen_y; j++) {
            zs[i][j] = +20; // a very large value intended as +INFINITY
        }
    }

    // perform scan conversion, populate the 2D array pixels
    // the array zs is the z-buffer.

    int sz=triangles.size(); 
    for(int i=0; i<sz; i++){ 
        if (triangles[i].clipped==false){ 
            Triangle t=triangles[i];
            homogeneous_point trPoints[triangleSide];
            trPoints[0]=t.p1;
            trPoints[1]=t.p2;
            trPoints[2]=t.p3;

            //for(int j=0;j<triangleSide;j++) trPoints[j].print();
            sort(trPoints, trPoints+triangleSide, comparePoint);
            //cout << endl;
            //for(int j=0;j<triangleSide;j++) trPoints[j].print();

            homogeneous_point hi=trPoints[2];
            homogeneous_point md=trPoints[1];
            homogeneous_point lo=trPoints[0];

            int hiIdx=0, loIdx=screen_y-1;
            if (hi.y<1) hiIdx=getYIdx(hi.y);
            if (lo.y>-1) loIdx=getYIdx(lo.y);
            //cout << "\nhiIdx:"<<hiIdx<<" loIdx:"<<loIdx<<endl;

            double gridY=2.0/screen_y, gridX=2.0/screen_x;
            double scanLine=1-(hiIdx*gridY+0.5*gridY);
            double za, zb, xa, xb;
            double xaIdx, xbIdx;
            double newz;

            //cout << gridY <<" "<<scanLine<< endl;
            for (int j=hiIdx; j<=loIdx; j++){ 
                if (scanLine >= md.y){
                    za=hi.z-(hi.z-md.z)*(hi.y-scanLine)/(hi.y-md.y);
                    zb=hi.z-(hi.z-lo.z)*(hi.y-scanLine)/(hi.y-lo.y);
                    xa=hi.x-((hi.x-md.x)*(hi.y-scanLine)/(hi.y-md.y));
                    xb=hi.x-(hi.x-lo.x)*(hi.y-scanLine)/(hi.y-lo.y);
                }
                else {
                    za=md.z-(md.z-lo.z)*(md.y-scanLine)/(md.y-lo.y);
                    zb=hi.z-(hi.z-lo.z)*(hi.y-scanLine)/(hi.y-lo.y);
                    xa=md.x-(md.x-lo.x)*(md.y-scanLine)/(md.y-lo.y);
                    xb=hi.x-(hi.x-lo.x)*(hi.y-scanLine)/(hi.y-lo.y); 
                }
                if (xa>xb){
                    swap(xa,xb);
                    swap(za,zb);
                }
                //cout <<"xa:"<<xa<<" xb:"<<xb<<" za:"<<za<<" zb:"<<zb<<endl;
                //cout <<"scanline y:"<< scanLine << endl;

                if (xa>=-1 || xb<=1){
                    xaIdx=0, xbIdx=screen_x-1; 
                    if(xa>=-1) xaIdx=getXidx(xa);
                    if(xb<=1) xbIdx=getXidx(xb);
                    double xscanLine=-1+xaIdx*gridX+0.5*gridX;
                    //cout << "xaidx:"<<xaIdx<<" xbidx:"<<xbIdx<<endl;

                    
                    for(int k=xaIdx; k<=xbIdx; k++){
                        newz=0;
                        if (k==xaIdx) newz=za;
                        else if (k==xbIdx) newz=zb;
                        else newz=zb-(zb-za)*(xb-xscanLine)/(xb-xa);

                        if (newz<zs[k][j]){ 
                            zs[k][j]=newz;
                            pixels[k][j]=t.c; 
                        }
                        xscanLine+=gridX;
                    }
                }
                scanLine-=gridY;
            }
        }
    }

    // the following code generates a bmp image. do not change this.
    bitmap_image image(screen_x, screen_y);
    for (int x = 0; x < screen_x; x++) {
        for (int y = 0; y < screen_y; y++) {
            image.set_pixel(x, y, pixels[x][y].r, pixels[x][y].g, pixels[x][y].b);
        }
    }
    image.save_image("out.bmp");

    // free the dynamically allocated memory

}

homogeneous_point getIntersectPoint(homogeneous_point p1,homogeneous_point p2, double zVal){
    Vector p(p1.x,p1.y,p1.z);
    Vector q(p2.x,p2.y,p2.z);
    Vector v=q-p;

    double t=(zVal-p.z)/v.z;
    homogeneous_point intersectPoint(p.x+t*v.x, p.y+t*v.y, p.z+t*v.z);
    return intersectPoint; 
}

void printVec(vector<homogeneous_point> v){
    for (int i=0;i<v.size();i++){
        v[i].print();
    }
    cout << endl;
}


void clip(){ 
    double nr=-near-.001;
    double fr=-far+.001;
    //double nr=-near;
    //double fr=-far;
    int sz=triangles.size();

    for(int i=0;i<sz;i++){
        if((triangles[i].p1.z>nr && triangles[i].p2.z>nr && triangles[i].p3.z>nr) ||
            (triangles[i].p1.z<fr && triangles[i].p2.z<fr && triangles[i].p3.z<fr)){
            triangles[i].clipped=true;
        }  
    }
     
    
    for(int i=0;i<sz;i++){
        if((triangles[i].p1.z<nr && triangles[i].p2.z<nr && triangles[i].p3.z<nr) &&
            (triangles[i].p1.z>fr && triangles[i].p2.z>fr && triangles[i].p3.z>fr)){
            //cout << "cont"<<endl;
            continue;
        } 
        if(triangles[i].clipped==false){
            Triangle t=triangles[i];
            homogeneous_point trPoints[triangleSide];
            trPoints[0]=t.p1;
            trPoints[1]=t.p2;
            trPoints[2]=t.p3;

            vector<homogeneous_point> outputLst;

            for(int j=0; j<triangleSide; j++){
                homogeneous_point a1=trPoints[j];
                homogeneous_point a2=trPoints[(j+1)%triangleSide];

                if (a1.z<=nr && a2.z<=nr) outputLst.push_back(a2);
                else if(a1.z<=nr && a2.z>nr) outputLst.push_back(getIntersectPoint(a1,a2,nr));
                else if(a1.z>nr && a2.z<=nr){
                    outputLst.push_back(getIntersectPoint(a1,a2,nr));
                    outputLst.push_back(a2);
                }
            } 

            vector<homogeneous_point> finalOutputLst;
            for(int j=0; j<outputLst.size(); j++){
                homogeneous_point a1=outputLst[j];
                homogeneous_point a2=outputLst[(j+1)%outputLst.size()];

                if (a1.z>=fr && a2.z>=fr) finalOutputLst.push_back(a2);
                else if(a1.z>=fr && a2.z<fr) finalOutputLst.push_back(getIntersectPoint(a1,a2,fr));
                else if(a1.z<fr && a2.z>=fr){
                    finalOutputLst.push_back(getIntersectPoint(a1,a2,fr));
                    finalOutputLst.push_back(a2);
                }
            } 
            
            homogeneous_point temp=finalOutputLst[0];
            //homogeneous_point temp(finalOutputLst[0].x,finalOutputLst[0].y,finalOutputLst[0].z);
            //cout << "i:" << i << endl;
            //temp.print(); 
            
            for (int j=1; j<finalOutputLst.size()-1;j++){
                homogeneous_point h1=finalOutputLst[j];
                homogeneous_point h2=finalOutputLst[j+1];

                Triangle newt(temp.x,temp.y,temp.z,h1.x,h1.y,h1.z,h2.x,h2.y,h2.z,t.c.r,t.c.g,t.c.b); 
                triangles.push_back(newt);
                triangles[i].clipped=true;
            }  
        }
    }
}


void stage3()
{
    if (near == far) return;
    ifstream stage2;
    ofstream stage3;
    stage2.open ("stage2.txt");
    stage3.open ("stage3.txt");
    stage3 << std::fixed;
    stage3 << std::setprecision(7);

    // process input from stage2 and write to stage3

    clip();

    fov_x=fov_y*aspectRatio;
    double t=near*tan(fov_y*0.5*pi/180);
    double r=near*tan(fov_x*0.5*pi/180);
    matrix projMatrix=matrix::make_identity(4);
    
    projMatrix.setIdx(0,0, near/r);
    projMatrix.setIdx(1,1, near/t);
    projMatrix.setIdx(2,2, (far+near)/(near-far));
    projMatrix.setIdx(2,3, 2*far*near/(near-far));
    projMatrix.setIdx(3,2, -1);
    projMatrix.setIdx(3,3, 0);

    for (int i=0; i<triangles.size(); i++){
        Triangle t=triangles[i];
        if (t.clipped==false){
            homogeneous_point h1=projMatrix*t.p1;
            homogeneous_point h2=projMatrix*t.p2;
            homogeneous_point h3=projMatrix*t.p3;

            stage3 << h1.x <<" "<< h1.y <<" "<< h1.z << endl;
            stage3 << h2.x <<" "<< h2.y <<" "<< h2.z << endl;
            stage3 << h3.x <<" "<< h3.y <<" "<< h3.z << endl;
            stage3 << endl;

            Triangle newt(h1.x,h1.y,h1.z,h2.x,h2.y,h2.z,h3.x,h3.y,h3.z,t.c.r,t.c.g,t.c.b);
            triangles[i]=newt;
        }
    }
    stage3.close();
    stage2.close();

}

void stage2()
{
    ifstream stage1;
    ofstream stage2;
    stage1.open ("stage1.txt");
    stage2.open ("stage2.txt");
    stage2 << std::fixed;
    stage2 << std::setprecision(7);

    // collect input from stage1 and process, write output to stage2
    Vector l(look_x-eye_x, look_y-eye_y, look_z-eye_z);
    l.normalize();
    Vector up(up_x,up_y,up_z);
    Vector r=Vector::cross(l,up);
    r.normalize();
    Vector u=Vector::cross(r,l);

    matrix T=matrix::getTranslationMatrix(-eye_x,-eye_y,-eye_z);
    matrix R=matrix::make_identity(4);
    R.setIdx(0,0,r.x);
    R.setIdx(0,1,r.y);
    R.setIdx(0,2,r.z);
    R.setIdx(1,0,u.x);
    R.setIdx(1,1,u.y);
    R.setIdx(1,2,u.z);
    R.setIdx(2,0,-l.x);
    R.setIdx(2,1,-l.y);
    R.setIdx(2,2,-l.z);
    matrix V=R*T;

    for (int i=0; i<triangles.size(); i++){
        Triangle t=triangles[i];
        homogeneous_point h1=V*t.p1;
        homogeneous_point h2=V*t.p2;
        homogeneous_point h3=V*t.p3;

        stage2 << h1.x <<" "<< h1.y <<" "<< h1.z << endl;
        stage2 << h2.x <<" "<< h2.y <<" "<< h2.z << endl;
        stage2 << h3.x <<" "<< h3.y <<" "<< h3.z << endl;
        stage2 << endl;

        Triangle newt(h1.x,h1.y,h1.z,h2.x,h2.y,h2.z,h3.x,h3.y,h3.z,t.c.r,t.c.g,t.c.b);
        triangles[i]=newt;
    }

    stage1.close();
    stage2.close();

}

void stage1()
{
    ifstream scene;
    ofstream stage1;
    scene.open ("scene.txt");
    stage1.open ("stage1.txt");
    stage1 << std::fixed;
    stage1 << std::setprecision(7);

    string command;

    scene >> eye_x >> eye_y >> eye_z;
    scene >> look_x >> look_y >> look_z;
    scene >> up_x >> up_y >> up_z;
    scene >> fov_y >> aspectRatio >> near >> far;
    scene >> screen_x >> screen_y;
    scene >> backgroud.r >> backgroud.g >> backgroud.b;

    // take other commands as input from scene in a loop
    // process accordingly
    // write to stage1
    double ax,ay,az,bx,by,bz,cx,cy,cz;
    int tempr,tempg,tempb;
    double px,py,pz,theta;
 
    stack<matrix> firstStack;
    firstStack.push(matrix::make_identity(4)); 
    
    while(true){
        scene >> command;

        if (command=="triangle"){
            scene >> ax >> ay >> az;
            scene >> bx >> by >> bz;
            scene >> cx >> cy >> cz;
            scene >> tempr >> tempg >> tempb;
 
            matrix transMatrx=firstStack.top(); 
            homogeneous_point h1(ax,ay,az);
            homogeneous_point h2(bx,by,bz);
            homogeneous_point h3(cx,cy,cz);

            h1=transMatrx*h1;
            h2=transMatrx*h2;
            h3=transMatrx*h3;

            Triangle t(h1.x,h1.y,h1.z,h2.x,h2.y,h2.z,h3.x,h3.y,h3.z,tempr,tempg,tempb);
            triangles.push_back(t);
            
            stage1 << h1.x <<" "<< h1.y <<" "<< h1.z << endl;
            stage1 << h2.x <<" "<< h2.y <<" "<< h2.z << endl;
            stage1 << h3.x <<" "<< h3.y <<" "<< h3.z << endl;
            stage1 << endl;
        }
        else if (command=="translate"){
            scene >> px >> py >> pz;
            matrix temp=matrix::getTranslationMatrix(px,py,pz);

            matrix topMatrx=firstStack.top();
            firstStack.pop();
            firstStack.push(topMatrx*temp); 
        }
        else if (command=="scale"){
            scene >> px >> py >> pz;
            matrix temp=matrix::getScaleMatrix(px,py,pz);

            matrix topMatrx=firstStack.top();
            firstStack.pop();
            firstStack.push(topMatrx*temp); 
        }
        else if (command=="rotate"){
            scene >> theta >> px >> py >> pz;
            matrix temp=matrix::getRotationMatrix(theta,px,py,pz); 

            matrix topMatrx=firstStack.top();
            firstStack.pop();
            firstStack.push(topMatrx*temp); 
        }
        else if (command=="push"){ 
            matrix topMatrx=firstStack.top();
            firstStack.push(topMatrx); 
        }
        else if (command=="pop"){ 
            if (firstStack.empty()==false) firstStack.pop(); 
        }
        else if (command=="end")
            break;

    }

    scene.close();
    stage1.close();

}

int main()
{
    cout << std::fixed;
    cout << std::setprecision(4);

    stage1();
    stage2();
    stage3();
    scan_convert();

    return 0;
}
