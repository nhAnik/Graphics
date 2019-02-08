#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include <windows.h>
#include <glut.h>
#include<vector>
#include<utility>
#include<fstream>
#include<string.h>
#include<queue>

#include "bitmap_image.hpp"

using namespace std;

typedef pair<int,int> pii;

#define pi (2*acos(0.0))


const int maxPixel=2000;
const int inf=50;
const int maxObjects=100;

double cameraHeight;
double cameraAngle;
double inc=0.5;

string fileName="input.txt";
bitmap_image b_img ("texture.bmp");


class point{
public:
	double x,y,z;

	point(double a, double b, double c){
        x=a; y=b; z=c;
	}

	point(){}

	void printPoint(){
        printf("%.2lf,%.2lf,%.2lf\n",x,y,z);
	}

	point operator+(const point& p){
	    point p1(x+p.x,y+p.y,z+p.z);
	    return p1;
	}

	point operator+(double m){
	    point p1(x+m,y+m,z+m);
	    return p1;
	}

	point operator-(const point& p){
	    point p1(x-p.x,y-p.y,z-p.z);
	    return p1;
	}

	point operator*(double m){
        point p1(m*x,m*y,m*z);
        return p1;
	}

	void normalize(){
        double magn=sqrt(x*x+y*y+z*z);
        x=x/magn; y=y/magn; z=z/magn;
	}
};
typedef pair<double,point> pdp;

class Color {
public:
    double r, g, b;
    Color(double r, double g, double b) {
        this->r = r;
        this->g = g;
        this->b = b;
    }
    Color() {
    }

    Color operator+(const Color& c1){
        Color c(r+c1.r,g+c1.g,b+c1.b);
        return c;
    }

    Color operator*(double m){
        Color c1(m*r,m*g,m*b);
        return c1;
    }

    void printColor(){
        printf("%.3lf,%.3lf,%.3lf\n",r,g,b);
	}
};

class Object{
public:
    point aPoint;
    Color c;
    double ambient,diffuse,specular,reflection;
    double shininess;
    double radius,width,height;
    int type;

    Object(point p,Color cl,double a,double d,double s,double r,double shini,double rad){
        aPoint=p;
        c=cl;
        ambient=a;diffuse=d;specular=s;reflection=r;
        shininess=shini;
        radius=rad;
        type=1;  //sphere
    }
    Object(point p,Color cl,double a,double d,double s,double r,double shini,double wid,double hi){
        aPoint=p;
        c=cl;
        ambient=a;diffuse=d;specular=s;reflection=r;
        shininess=shini;
        width=wid;height=hi;
        type=2;   //pyramid
    }
};

class Light{
public:
    point aPoint;
    double falloff;
    Light(point p,double f){aPoint=p; falloff=f;}
};

class SpotLight:public Light{
public:
    point look;
    double cutoff;
    SpotLight(point p,double f,point lk,double c):Light(p,f){
        look=lk;
        cutoff=c;
    }
};

point pos(100,100,0), u(0,0,1), r(-.707,.707,0), l(-.707,-.707,0);
point temp;

double nearPlane, farPlane, fovy, aspectRatio;
int checkerWidth,levOfRecursion, numOfPixels;
double checkerAmbient, checkerDiffuse, checkerReflection;
int numOfObjects, numOfLights, numOfSpotLights;

vector<Object> objectArr;
vector<Light> lightArr;
vector<SpotLight> spotLightArr;

point pointBuffer[maxPixel][maxPixel];
bool visited[maxPixel][maxPixel];

int textureMode;
Color **textureBuffer;
int texHeight, texWidth;

point getCrossProd(point a, point b){

    point temp;
    temp.x=a.y*b.z-a.z*b.y;
    temp.y=a.z*b.x-a.x*b.z;
    temp.z=a.x*b.y-a.y*b.x;
    return temp;
}

double getDotProd(point p1,point p2){
    double d=p1.x*p2.x+p1.y*p2.y+p1.z*p2.z;
    return d;
}

point getVec(point p1, point p2){
    point v(p2.x-p1.x,p2.y-p1.y,p2.z-p1.z);
    return v;
}

point getUnitVec(point p1, point p2){
    double a=p2.x-p1.x;
    double b=p2.y-p1.y;
    double c=p2.z-p1.z;
    double mag=sqrt(a*a+b*b+c*c);
    a=a/mag;
    b=b/mag;
    c=c/mag;
    point v(a,b,c);
    return v;
}

double getDistOfPoints(point p1,point p2){
    double d=sqrt((p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y)+(p1.z-p2.z)*(p1.z-p2.z));
    return d;
}


point rotateVec(point a, point p){
    a.x=a.x*cos(cameraAngle)+p.x*sin(cameraAngle);
    a.y=a.y*cos(cameraAngle)+p.y*sin(cameraAngle);
    a.z=a.z*cos(cameraAngle)+p.z*sin(cameraAngle);
    a.normalize();
    return a;
}

point rotateVecAng(point a, point p, double angle){
    a.x=a.x*cos(angle)+p.x*sin(angle);
    a.y=a.y*cos(angle)+p.y*sin(angle);
    a.z=a.z*cos(angle)+p.z*sin(angle);
    return a;
}


void readFile(){
    ifstream inFile;
    inFile.open(fileName.c_str());

    inFile >> nearPlane >> farPlane >> fovy >> aspectRatio;
    inFile >> levOfRecursion >> numOfPixels;
    inFile >> checkerWidth >> checkerAmbient >> checkerDiffuse >> checkerReflection;
    inFile >> numOfObjects;

    string objType;
    double xx,yy,zz,rr,gg,bb,a,d,s,r,rad,wid,hi,shini,falof,cutof,lx,ly,lz;

    for(int i=0;i<numOfObjects;i++){
        inFile >> objType;
        if (objType=="sphere"){
            inFile >> xx >> yy >> zz;
            point tempPoint(xx,yy,zz);
            inFile >> rad;
            inFile >> rr >> gg >> bb;
            Color tempClr(rr,gg,bb);
            inFile >> a >> d >> s >> r;
            inFile >> shini;
            Object obj(tempPoint,tempClr,a,d,s,r,shini,rad);
            objectArr.push_back(obj);

        }
        else if (objType=="pyramid"){
            inFile >> xx >> yy >> zz;
            point tempPoint(xx,yy,zz);
            inFile >> wid >> hi;
            inFile >> rr >> gg >> bb;
            Color tempClr(rr,gg,bb);
            inFile >> a >> d >> s >> r;
            inFile >> shini;
            Object obj(tempPoint,tempClr,a,d,s,r,shini,wid,hi);
            objectArr.push_back(obj);
        }
    }

    inFile >> numOfLights;
    for(int i=0;i<numOfLights;i++){
        inFile >> xx >> yy >> zz >> falof;
        point tempPoint(xx,yy,zz);
        Light tempLight(tempPoint,falof);
        lightArr.push_back(tempLight);
    }

    inFile >> numOfSpotLights;
    for(int i=0;i<numOfSpotLights;i++){
        inFile >> xx >> yy >> zz >> falof;
        inFile >> lx >> ly >> lz >> cutof;
        point tempPoint(xx,yy,zz);
        point tempLook(lx,ly,lz);
        SpotLight tempLight(tempPoint,falof,tempLook,cutof);
        spotLightArr.push_back(tempLight);
    }
    inFile.close();
    //numOfPixels=5;
}

void drawSquare(point p,double width,int color){
    if (color==0) glColor3f(0.0,0.0,0.0);
    else if (color==1) glColor3f(1.0,1.0,1.0);
	glBegin(GL_QUADS);{
		glVertex3f(p.x, p.y, p.z);
		glVertex3f(p.x+width, p.y, p.z);
		glVertex3f(p.x+width, p.y+width, p.z);
		glVertex3f(p.x, p.y+width, p.z);

	}glEnd();
}

void drawChecker(int checkerWidth){

    int color=0,color1;

    for(int i=-inf; i<inf; i++){
        color1=color;
        for(int j=-inf; j<inf; j++){
            point tempPoint(i*checkerWidth,j*checkerWidth,0);
            drawSquare(tempPoint,checkerWidth,color1);
            color1=1-color1;
        }
        color=1-color;
    }
}

void drawSphere(point center,double radius,Color c){
    int slices=24,stacks=20;
	struct point points[100][100];
	int i,j;
	double h,r;

	for(i=0;i<=stacks;i++)
	{
		h=radius*sin(((double)i/(double)stacks)*(pi/2));
		r=radius*cos(((double)i/(double)stacks)*(pi/2));
		for(j=0;j<=slices;j++)
		{
			points[i][j].x=center.x + r*cos(((double)j/(double)slices)*2*pi);
			points[i][j].y=center.y + r*sin(((double)j/(double)slices)*2*pi);
			points[i][j].z=h;
		}
	}
	for(i=0;i<stacks;i++)
	{
        glColor3f(c.r,c.g,c.b);
		for(j=0;j<slices;j++)
		{
			glBegin(GL_QUADS);{

				glVertex3f(points[i][j].x,     points[i][j].y,     center.z+points[i][j].z);
				glVertex3f(points[i][j+1].x,   points[i][j+1].y,   center.z+points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x, points[i+1][j+1].y, center.z+points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,   points[i+1][j].y,   center.z+points[i+1][j].z);

                glVertex3f(points[i][j].x,     points[i][j].y,     center.z-points[i][j].z);
				glVertex3f(points[i][j+1].x,   points[i][j+1].y,   center.z-points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x, points[i+1][j+1].y, center.z-points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,   points[i+1][j].y,   center.z-points[i+1][j].z);
			}glEnd();
		}
	}
}

void drawPyramid(point p,double width,double height,Color c){

    point topPoint(p.x+width/2,p.y+width/2,p.z+height);
    glColor3f(c.r,c.g,c.b);
	glBegin(GL_QUADS);{
		glVertex3f(p.x, p.y, p.z);
		glVertex3f(p.x+width, p.y, p.z);
		glVertex3f(p.x+width, p.y+width, p.z);
		glVertex3f(p.x, p.y+width, p.z);
	}glEnd();

    glBegin(GL_TRIANGLES);
    {
        glVertex3f(topPoint.x,topPoint.y,topPoint.z);
        glVertex3f(p.x, p.y, p.z);
		glVertex3f(p.x+width, p.y, p.z);
    }
    glEnd();
    glBegin(GL_TRIANGLES);
    {
        glVertex3f(topPoint.x,topPoint.y,topPoint.z);
        glVertex3f(p.x+width, p.y, p.z);
		glVertex3f(p.x+width, p.y+width, p.z);
    }
    glEnd();
    glBegin(GL_TRIANGLES);
    {
        glVertex3f(topPoint.x,topPoint.y,topPoint.z);
        glVertex3f(p.x+width, p.y+width, p.z);
		glVertex3f(p.x, p.y+width, p.z);
    }
    glEnd();
    glBegin(GL_TRIANGLES);
    {
        glVertex3f(topPoint.x,topPoint.y,topPoint.z);
        glVertex3f(p.x, p.y+width, p.z);
        glVertex3f(p.x, p.y, p.z);
    }
    glEnd();

}

void drawAll(){
    drawChecker(checkerWidth);
    Color lightColor(0.9,0.9,0.9);
    Color spotLightColor(0.2,0.9,0.9);

    for (int i=0;i<numOfObjects;i++){
        Object obj=objectArr[i];
        if(obj.type==1) drawSphere(obj.aPoint,obj.radius,obj.c);
        else if(obj.type==2) drawPyramid(obj.aPoint,obj.width,obj.height,obj.c);
    }
    for (int i=0;i<numOfLights;i++){
        Light lit=lightArr[i];
        drawSphere(lit.aPoint,5,lightColor);
    }
    for (int i=0;i<numOfSpotLights;i++){
        SpotLight sptlit=spotLightArr[i];
        drawSphere(sptlit.aPoint,10,spotLightColor);
    }

}

void genPointBuffer(){
    double screenY, screenX, fovx;
    double pixelWidY, pixelWidX;

    screenY=2*nearPlane*tan(fovy*pi/360.0);
    fovx=fovy*aspectRatio;
    screenX=2*nearPlane*tan(fovx*pi/360.0);

    point midPoint=pos+l*nearPlane;

    pixelWidY=screenY/numOfPixels;
    pixelWidX=screenX/numOfPixels;

    int dirX[]={1,-1, 0, 0};
    int dirY[]={0, 0, 1,-1};

    memset(visited, false, sizeof(visited));
    queue<pii> q;

    pii start(numOfPixels/2,numOfPixels/2);
    q.push(start);
    visited[start.first][start.second]=true;
    pointBuffer[start.first][start.second]=midPoint;

    while(!q.empty()){
        pii v=q.front();
        q.pop();
        for(int i=0;i<4;i++){
            pii temp(v.first+dirY[i], v.second+dirX[i]);
            if(visited[temp.first][temp.second]==false && (temp.first>=0 && temp.first<numOfPixels) &&
               (temp.second>=0 && temp.second<numOfPixels)){

                visited[temp.first][temp.second]=true;
                if (i==0) pointBuffer[temp.first][temp.second]=pointBuffer[v.first][v.second]+r*pixelWidX;
                else if (i==1) pointBuffer[temp.first][temp.second]=pointBuffer[v.first][v.second]+r*(-pixelWidX);
                else if (i==2) pointBuffer[temp.first][temp.second]=pointBuffer[v.first][v.second]+u*pixelWidY;
                else if (i==3) pointBuffer[temp.first][temp.second]=pointBuffer[v.first][v.second]+u*(-pixelWidY);
                q.push(temp);
            }
        }
    }
}


double getDeterminants(double arr[3][3]){
    double det1=arr[0][0]*(arr[1][1]*arr[2][2]-arr[1][2]*arr[2][1]);
    double det2=arr[0][1]*(arr[1][0]*arr[2][2]-arr[2][0]*arr[1][2]);
    double det3=arr[0][2]*(arr[1][0]*arr[2][1]-arr[1][1]*arr[2][0]);
    return det1-det2+det3;
}

double getTriangleT(point a,point b,point c,point p,point vec){
    double beta, gamma, t=-1, det;
    double A[3][3] ={{a.x-b.x, a.x-c.x, vec.x},  {a.y-b.y, a.y-c.y, vec.y},  {a.z-b.z, a.z-c.z, vec.z}};
    double A1[3][3]={{a.x-p.x, a.x-c.x, vec.x},  {a.y-p.y, a.y-c.y, vec.y},  {a.z-p.z, a.z-c.z, vec.z}};
    double A2[3][3]={{a.x-b.x, a.x-p.x, vec.x},  {a.y-b.y, a.y-p.y, vec.y},  {a.z-b.z, a.z-p.z, vec.z}};
    double A3[3][3]={{a.x-b.x, a.x-c.x, a.x-p.x},{a.y-b.y, a.y-c.y, a.y-p.y},{a.z-b.z, a.z-c.z, a.z-p.z}};

    det=getDeterminants(A);
    if (det!=0){
        beta=getDeterminants(A1)/det;
        gamma=getDeterminants(A2)/det;
        t=getDeterminants(A3)/det;

        if(beta+gamma<1 && beta>0 && gamma>0 && t>0) return t;
        else return -1;
    }


    return t;

}

double getRectangleT(point l,double width, point p,point vec){

    double t=(l.z-p.z)/vec.z;
    if (t<0) printf("l.z %lf p.z %lf and vec.z %lf",l.z,p.z,vec.z);
    point intersect=p+vec*t;
    //intersect.printPoint();
    if(intersect.x>=l.x && intersect.x<=(l.x+width) && intersect.y>=l.y && intersect.y<=(l.y+width) ) {
        //printf("four points %lf %lf %lf %lf\n",l.x,l.x+width,l.y,l.y+width);
        return t;
    }
    return -1;
}

point getNormalOnPlane(point a,point b,point c){
    point p1=getVec(c,b);
    point p2=getVec(c,a);
    point p=getCrossProd(p1,p2);
    p.normalize();
    return p;
}

pdp getPyramidT(point p,point vec,point lowerLeft,double width,double height){

    double tempT[6],minT=100000;
    int minIdx;
    point normal(0,0,0);

    point topPoint(lowerLeft.x+width/2, lowerLeft.y+width/2, lowerLeft.z+height);
    point upperLeft(lowerLeft.x,        lowerLeft.y+width,   lowerLeft.z);
    point lowerRight(lowerLeft.x+width, lowerLeft.y ,        lowerLeft.z);
    point upperRight(lowerLeft.x+width, lowerLeft.y+width,   lowerLeft.z);

    tempT[0]=getTriangleT(lowerLeft,  upperLeft,  topPoint,   p,vec);
    tempT[1]=getTriangleT(upperLeft,  upperRight, topPoint,   p,vec);
    tempT[2]=getTriangleT(upperRight, lowerRight, topPoint,   p,vec);
    tempT[3]=getTriangleT(lowerRight, lowerLeft,  topPoint,   p,vec);
    tempT[4]=getTriangleT(lowerLeft,  upperLeft,  upperRight, p,vec);
    tempT[5]=getTriangleT(lowerLeft,  upperRight, lowerRight, p,vec);
    //tempT[4]=getRectangleT(lowerLeft,width,p,vec);

    if(tempT[0]==-1 && tempT[1]==-1 && tempT[2]==-1 && tempT[3]==-1 && tempT[4]==-1 && tempT[5]==-1) {
        minT=-1;
        return pdp(minT,normal);
    }
    else {
        for(int i=0;i<6;i++){
            if(tempT[i]!=-1 && tempT[i]<minT){
                minT=tempT[i];
                minIdx=i;
            }
        }
    }
    if(minIdx==0) normal=getNormalOnPlane(upperLeft,lowerLeft,topPoint);
    else if(minIdx==1) normal=getNormalOnPlane(upperRight,upperLeft,topPoint);
    else if(minIdx==2) normal=getNormalOnPlane(lowerRight,upperRight,topPoint);
    else if(minIdx==3) normal=getNormalOnPlane(lowerLeft,lowerRight,topPoint);
    else if(minIdx==4 || minIdx==5) normal.z=-1;

    return pdp(minT,normal);
}

pdp getCircleT(point p,point vec,point center,double rad){
    double a,b,c,t1,t2,t=-1;
    point normal(0,0,0);
    point temp=p+center*(-1);
    a=1;
    b=2*(vec.x*temp.x+vec.y*temp.y+vec.z*temp.z);
    c=(temp.x*temp.x+temp.y*temp.y+temp.z*temp.z)-rad*rad;
    double discrim=b*b-4*a*c;
    if (discrim>=0){
        t1=(-b+sqrt(discrim))/(2.0*a);
        t2=(-b-sqrt(discrim))/(2.0*a);


        if(t1>0 && t2>0) t=min(t1,t2);
        else{
            if(t1>0) t=t1;
            else if(t2>0) t=t2;
            else t=-1;
        }
    }
    if (t!=-1) normal=getUnitVec(center,p+vec*t);

    return pdp(t,normal);
}

bool isIlluminated(point p,point vec){
    double tempT[maxObjects];
    for(int i=0;i<numOfObjects;i++){
        if(objectArr[i].type==1)
            tempT[i]=getCircleT(p,vec,objectArr[i].aPoint,objectArr[i].radius).first;

        else if(objectArr[i].type==2)
            tempT[i]=getPyramidT(p,vec,objectArr[i].aPoint,objectArr[i].width,objectArr[i].height).first;
    }

    bool allAreMinus=true;
    for(int i=0;i<numOfObjects;i++) if(tempT[i]!=-1) allAreMinus=false;

    return allAreMinus;
}

bool isIlluminatedSpotLight(point p,point vec,SpotLight spl){
    if(isIlluminated(p,vec)==false) return false;

    point srcToPoint=getUnitVec(spl.aPoint,p);
    point srcTolk=getUnitVec(spl.aPoint,spl.look);
    double angle=acos(getDotProd(srcToPoint,srcTolk));
    angle=angle*180.0/pi;
    if (angle>spl.cutoff) return false;
    return true;
}

Color getPixelColor(point p,point vec,int depth){  //return value 0 represents black and 1 represents white color

    Color black(0.0,0.0,0.0);
    Color white(1.0,1.0,1.0);

    //if(depth==0) return black;
    Color finalColor;
    int typeOfPixel;
    //typeOfPixel=-2 means blank -1 means checkerboard and another positive value means respective
    //index of objectArr

    double minT=100000,minIdx,tempT[maxObjects];
    point normalArr[maxObjects];
    point normalOnCheck(0,0,1);
    //stores the normal vector

    for(int i=0;i<numOfObjects;i++){
        if(objectArr[i].type==1){
            pdp tempPair=getCircleT(p,vec,objectArr[i].aPoint,objectArr[i].radius);
            tempT[i]=tempPair.first;
            normalArr[i]=tempPair.second;
        }

        else if(objectArr[i].type==2){ ;
            pdp tempPair=getPyramidT(p,vec,objectArr[i].aPoint,objectArr[i].width,objectArr[i].height);
            tempT[i]=tempPair.first;
            normalArr[i]=tempPair.second;
        }
    }

    bool allAreMinus=true;
    for(int i=0;i<numOfObjects;i++) if(tempT[i]!=-1) allAreMinus=false;

    if (allAreMinus) minT=-1;
    else {
        for(int i=0;i<numOfObjects;i++){
            if(tempT[i]!=-1 && tempT[i]<minT){
                minT=tempT[i];
                minIdx=i;
            }
        }
    }

    if(vec.z==0){
        if(minT==-1) {finalColor=black; typeOfPixel=-2;}
        else {finalColor=objectArr[minIdx].c; typeOfPixel=minIdx;}
    }
    else {
        double checkerT=-p.z/vec.z;
        point pointOnChecker=p+vec*checkerT;

        if(getDistOfPoints(p,pointOnChecker)>farPlane || checkerT<0) checkerT=-1;

        if(checkerT==-1){
            if(minT==-1) {finalColor=black; typeOfPixel=-2;}
            else {finalColor=objectArr[minIdx].c; typeOfPixel=minIdx;}
        }
        else {
            int y=int((pointOnChecker.y+inf*checkerWidth)/checkerWidth)-inf;
            int x=int((pointOnChecker.x+inf*checkerWidth)/checkerWidth)-inf;

            int colStart,colorType;
            Color texColor;
            if(textureMode==1){
                double texH=(pointOnChecker.y+inf*checkerWidth)-((y+inf)*checkerWidth);
                double texW=(pointOnChecker.x+inf*checkerWidth)-((x+inf)*checkerWidth);
                int hIdx=int(texH*texHeight/checkerWidth);
                int wIdx=int(texW*texWidth/checkerWidth);
                texColor=textureBuffer[wIdx][hIdx];
            }
            else{
                if(abs(x)%2==0) colStart=0;
                else colStart=1;

                if(abs(y)%2==0) colorType=colStart;
                else colorType=1-colStart;
            }

            if(minT==-1){
                minT=checkerT;
                if(textureMode==1) {finalColor=texColor; typeOfPixel=-1;}
                else if(colorType==0) {finalColor=black; typeOfPixel=-1;}
                else if(colorType==1) {finalColor=white; typeOfPixel=-1;}
            }
            else {
                if (checkerT<minT){
                    minT=checkerT;
                    if(textureMode==1) {finalColor=texColor; typeOfPixel=-1;}
                    else if(colorType==0) {finalColor=black; typeOfPixel=-1;}
                    else if(colorType==1) {finalColor=white; typeOfPixel=-1;}
                }
                else {finalColor=objectArr[minIdx].c; typeOfPixel=minIdx;}
            }
        }

    }

    if (typeOfPixel==-2) return black;

    Color lightedColor=finalColor;
    Color reflectedColor(0,0,0);
    point reflectionPoint,pixelPoint, lightPos;
    point intersectionPoint=p+vec*minT;
    point vecToSource;
    double fall;
    bool decision;

    double lambert=0.0, phong=0.0;
    for(int i=0;i<numOfLights+numOfSpotLights;i++){
        pixelPoint=intersectionPoint-vec*0.01;

        if (i<numOfLights){
            Light light=lightArr[i];
            lightPos=light.aPoint;
            fall=light.falloff;
            vecToSource=getUnitVec(pixelPoint,lightPos);
            decision=isIlluminated(pixelPoint,vecToSource);
        }
        else{
            SpotLight spotLight=spotLightArr[i-numOfLights];
            lightPos=spotLight.aPoint;
            fall=spotLight.falloff;
            vecToSource=getUnitVec(pixelPoint,lightPos);
            decision=isIlluminatedSpotLight(pixelPoint,vecToSource,spotLight);
        }
        point reflected(0,0,0);


        if (decision){
            double d=getDistOfPoints(pixelPoint,lightPos);
            double scaleFact=exp(-d*d*fall);
            if(typeOfPixel>=0){
                point normal=normalArr[typeOfPixel];
                lambert+=max(getDotProd(vecToSource,normal)*scaleFact,0.0);

                double factor=2*getDotProd(vecToSource,normal);
                reflected=normal*factor+vecToSource*(-1);
                reflected.normalize();
                phong+=max(pow(getDotProd(reflected,vec*(-1)),objectArr[typeOfPixel].shininess)*scaleFact, 0.0);
            }
            else if(typeOfPixel==-1){
                lambert+=max(getDotProd(vecToSource,normalOnCheck)*scaleFact, 0.0);
            }
        }
    }
    point reflectionRay;
    if(depth<levOfRecursion){
        if (typeOfPixel==-1)
            reflectionRay=vec+normalOnCheck*(-2*getDotProd(vec,normalOnCheck));
        else if(typeOfPixel>=0){
            point normal=normalArr[typeOfPixel];
            reflectionRay=vec+normal*(-2*getDotProd(vec,normal));
        }
        reflectionRay.normalize();
        reflectedColor=getPixelColor(intersectionPoint+reflectionRay*0.0001,reflectionRay,depth+1); ;
    }

    if(typeOfPixel==-1) lightedColor=finalColor*checkerAmbient
        +finalColor*(lambert*checkerDiffuse)
        +reflectedColor*checkerReflection;


    else if (typeOfPixel>=0) lightedColor=
         finalColor*objectArr[typeOfPixel].ambient
        +finalColor*(lambert*objectArr[typeOfPixel].diffuse)
        +finalColor*(phong*objectArr[typeOfPixel].specular)
        +reflectedColor*objectArr[typeOfPixel].reflection;


    return lightedColor;
}

void traceTheRay(){

    Color background(0.0,0.0,0.0);
    Color white(1.0,1.0,1.0);
    Color** pixels = new Color*[numOfPixels];
    for (int i=0;i<numOfPixels;i++) {
        pixels[i]=new Color[numOfPixels];
        for (int j=0;j<numOfPixels;j++)
            pixels[i][j]=background;
    }


    texHeight = b_img.height();
    texWidth = b_img.width();
    textureBuffer = new Color* [texWidth];
    for (int i = 0; i < texWidth; i++) {
        textureBuffer[i] = new Color [texHeight];
        for (int j = 0; j < texHeight; j++) {
            unsigned char r, g, b;
            b_img.get_pixel(i, j, r, g, b);
            Color c(r/255.0, g/255.0, b/255.0);
            textureBuffer[i][j] = c;
        }
    }

    genPointBuffer();

    for(int i=0;i<numOfPixels;i++){
        for(int j=0;j<numOfPixels;j++){
            point gridPoint=pointBuffer[i][j];
            point ray=getUnitVec(pos,gridPoint);
            pixels[j][numOfPixels-i]=getPixelColor(gridPoint,ray,0);
        }
    }

    bitmap_image image(numOfPixels,numOfPixels);
    for (int i=0;i<numOfPixels;i++) {
        for (int j=0;j<numOfPixels;j++) {
            image.set_pixel(i, j, min(pixels[i][j].r*255,255.0), min(pixels[i][j].g*255,255.0), min(pixels[i][j].b*255,255.0));
        }
    }
    image.save_image("out.bmp");

    printf("Ray tracing completed\n");
}

void keyboardListener(unsigned char key, int x,int y){
	switch(key){
        case '0':
            traceTheRay();
            break;

	    case '1':
			temp=getCrossProd(u,l);
			temp.normalize();
			l=rotateVec(l,temp);
			temp=getCrossProd(u,r);
			temp.normalize();
			r=rotateVec(r,temp);
			break;

		case '2':
			temp=getCrossProd(l,u);
			temp.normalize();
			l=rotateVec(l,temp);
			temp=getCrossProd(r,u);
			temp.normalize();
			r=rotateVec(r,temp);
			break;

        case '3':
            temp=getCrossProd(r,l);
            temp.normalize();
            l=rotateVec(l,temp);
            temp=getCrossProd(r,u);
            temp.normalize();
            u=rotateVec(u,temp);
            break;

        case '4':
            temp=getCrossProd(l,r);
            temp.normalize();
            l=rotateVec(l,temp);
            temp=getCrossProd(u,r);
            temp.normalize();
            u=rotateVec(u,temp);
            break;

        case '5': //counterclock
            temp=getCrossProd(l,u);
            temp.normalize();
            u=rotateVec(u,temp);
            temp=getCrossProd(l,r);
            temp.normalize();
            r=rotateVec(r,temp);
            break;

        case '6':
            temp=getCrossProd(u,l);
            temp.normalize();
            u=rotateVec(u,temp);
            temp=getCrossProd(r,l);
            temp.normalize();
            r=rotateVec(r,temp);
            break;

        case ' ':
            textureMode=1-textureMode;
            if (textureMode==1) printf("Texture Mode Enabled\n");
            else if (textureMode==0) printf("Texture Mode Disabled\n");

		default:
			break;
	}
}


void specialKeyListener(int key, int x,int y){
	switch(key){
	    case GLUT_KEY_UP:		// up arrow key
			pos.x+=l.x;
			pos.y+=l.y;
			pos.z+=l.z;
			break;
		case GLUT_KEY_DOWN:		//down arrow key
			pos.x-=l.x;
			pos.y-=l.y;
			pos.z-=l.z;
			break;


		case GLUT_KEY_RIGHT:
			pos.x+=r.x;
			pos.y+=r.y;
			pos.z+=r.z;
			break;
		case GLUT_KEY_LEFT:
			pos.x-=r.x;
			pos.y-=r.y;
			pos.z-=r.z;
			break;

		case GLUT_KEY_PAGE_UP:
		    pos.x+=u.x;
			pos.y+=u.y;
			pos.z+=u.z;
			break;
		case GLUT_KEY_PAGE_DOWN:
		    pos.x-=u.x;
			pos.y-=u.y;
			pos.z-=u.z;
			break;


		case GLUT_KEY_INSERT:
			break;

		case GLUT_KEY_HOME:
			break;
		case GLUT_KEY_END:
			break;

		default:
			break;
	}
}


void mouseListener(int button, int state, int x, int y){	//x, y is the x-y of the screen (2D)
	switch(button){
		case GLUT_LEFT_BUTTON:
			//.......
			break;

		case GLUT_RIGHT_BUTTON:
			//........
			break;

		case GLUT_MIDDLE_BUTTON:
			//........
			break;

		default:
			break;
	}
}



void display(){

	//clear the display
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(0,0,0,0);	//color black
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	/********************
	/ set-up camera here
	********************/
	//load the correct matrix -- MODEL-VIEW matrix
	glMatrixMode(GL_MODELVIEW);

	//initialize the matrix
	glLoadIdentity();

	//now give three info
	//1. where is the camera (viewer)?
	//2. where is the camera looking?
	//3. Which direction is the camera's UP direction?

	//gluLookAt(100,100,100,	0,0,0,	0,0,1);
	gluLookAt(pos.x, pos.y, pos.z, pos.x+l.x, pos.y+l.y, pos.z+l.z, u.x, u.y, u.z);
	//gluLookAt(0,0,200,	0,0,0,	0,1,0);


	//again select MODEL-VIEW
	glMatrixMode(GL_MODELVIEW);


	/****************************
	/ Add your objects from here
	****************************/
	//add objects
	//drawChecker();
    //drawSphere(center, 20.0, sphereColor);
    //drawPyramid(center,width,height,sphereColor);
    drawAll();

	//ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
	glutSwapBuffers();
}


void animate(){
	//codes for any changes in Models, Camera

	glutPostRedisplay();
}

void init(){
	//codes for initialization
	cameraHeight=150.0;
	cameraAngle=0.03;
	textureMode=0;

	//clear the screen
	glClearColor(0,0,0,0);

	readFile();

	/************************
	/ set-up projection here
	************************/
	//load the PROJECTION matrix
	glMatrixMode(GL_PROJECTION);

	//initialize the matrix
	glLoadIdentity();

	//give PERSPECTIVE parameters
	//gluPerspective(80,	1,	1,	1000.0);
	gluPerspective(fovy,aspectRatio,nearPlane,farPlane);
	//field of view in the Y (vertically)
	//aspect ratio that determines the field of view in the X direction (horizontally)
	//near distance
	//far distance
}

int main(int argc, char **argv){
	glutInit(&argc,argv);
	glutInitWindowSize(500, 500);
	glutInitWindowPosition(0, 0);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

	glutCreateWindow("My OpenGL Program");

	init();

	glEnable(GL_DEPTH_TEST);	//enable Depth Testing

	glutDisplayFunc(display);	//display callback function
	glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occuring)

	glutKeyboardFunc(keyboardListener);
	glutSpecialFunc(specialKeyListener);
	glutMouseFunc(mouseListener);

	glutMainLoop();		//The main loop of OpenGL

	return 0;
}
