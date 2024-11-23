#include<bits/stdc++.h>
using namespace std;
#define ld long double
#define ll long long
const ld eps = 1e-10; // 精度
const ld pi = acos(-1); //π

struct Point{
    ld x,y;
    ll id;
    Point(ld x=0, ld y=0):x(x),y(y){}
};
typedef Point Vector;
typedef pair<Point, Point> ppp;

bool eq(ld a, ld b){return fabs(a-b)<eps;} // 判断浮点数是否相等
ld sign(ld a){if(eq(a,0)) return 0; return a>0?1:-1;} // 符号
bool operator == (Point a, Point b){return eq(a.x,b.x) && eq(a.y,b.y);} // 判断点是否相等
bool operator < (Point a, Point b){return eq(a.x,b.x)?(sign(a.y-b.y)<0):(sign(a.x-b.x)<0);} // 默认按x排序
bool cmpy(Point a, Point b){return eq(a.y,b.y)?(sign(a.x-b.x)<0):(sign(a.y-b.y)<0);} // 按y坐标排序
ld lmax(ld a, ld b){return sign(a-b)>0?a:b;} // 精确求最大值

Vector operator + (Vector a, Vector b){return Vector(a.x+b.x,a.y+b.y);}
Vector operator - (Vector a, Vector b){return Vector(a.x-b.x,a.y-b.y);}
Vector operator * (Vector a, ld k){return Vector(a.x*k,a.y*k);}
Vector operator / (Vector a, ld k){return Vector(a.x/k,a.y/k);}

ld operator * (Vector a,Vector b){return a.x*b.x+a.y*b.y;} // 点乘
ld cross(Vector a,Vector b){return a.x*b.y-a.y*b.x;} // 叉乘
ld length(Vector a){return sqrt(a*a);} // 求向量长度
ld distance(Point a, Point b){return length(a-b);} // 求两点距离
ld distance(ppp a){return distance(a.first,a.second);} // 求两点距离

ld angle(Vector a, Vector b){ // 求向量a和b的夹角
    return acos((a*b)/length(a)/length(b));
}
ld rotation_angle(Vector a, Vector b){ // 求向量a逆时针旋转到向量b的角度
    ld res = atan2(cross(a,b),a*b);
    return sign(res)<0?res+2*pi:res;
}
Vector rotate(Vector a, ld rad){ // 逆时针旋转rad弧度
    return Vector(a.x*cos(rad)-a.y*sin(rad),a.x*sin(rad)+a.y*cos(rad));
}
Vector normal(Vector a){ // 求向量a的法向量
    return Vector(-a.y,a.x);
}
Vector projection(Vector a, Vector b){ // 求向量a在向量b上的投影向量
    return b*(a*b)/(b*b);
}
ld area(Point *p, int n){ // 凸多边形面积
    ld res = 0;
    for(int i=1;i<=n;i++) res += cross(p[i],p[i%n+1]);
    return abs(res)/2; // 避免转错方向
}


struct Line{
    Point u,v;
    Vector dir; // 方向向量
    Line():u(Point()),v(Point()),dir(Vector()){}
    Line(Point u, Point v):u(u),v(v),dir(v-u){} // 两点式
    Line(Point u, ld k):u(u),v(u+Point(1,k)),dir(v-u){} // 点斜式
    Line(ld k, ld b){
        u = Point(0,b);
        v = Point(1,k+b);
        dir = v-u;
    } // 斜截式
};
Point cross_point(Line a, Line b){ // 求两直线交点
    ld s1 = cross(a.v-a.u,b.u-a.u);
    ld s2 = cross(a.v-a.u,b.v-a.u);
    return (b.u*s2-b.v*s1)/(s2-s1);
} // 类似定比分点公式
ld distance(Point p, Line l){ // 点到直线距离
    return abs(cross(l.dir,p-l.u)/length(l.dir));
}
ld distance(Point p, Point A, Point B){ // 点到线段AB
    if(A==B) return distance(p,A);
    Vector v1 = B-A, v2 = p-A, v3 = p-B;
    if(sign(v1*v2)<0) return length(v2);
    if(sign(v1*v3)>0) return length(v3);
    return fabs(cross(v1,v2)/length(v1));
}
Point projection(Point p, Line l) {
    Vector v = l.dir;
    double ratio = (v*(p - l.u))/(v*v);
    return l.u + v*ratio;
}
Point reflection(Point p, Line l){ // 点关于直线的对称点
    return p+(projection(p,l)-p)*2;
}
bool in_line(Point p,Line l){ // 判断点p是否在直线上
    return eq(cross(l.v-l.u,p-l.u),0);
}
bool in_segment(Point p, Point A, Point B){ // 判断点p是否在线段AB上
    return eq(cross(A-p,B-p),0) && sign((A-p)*(B-p))<=0; // 在直线上 且 在AB之间
}
bool is_cross(Point a, Point b, Point c, Point d){ // 判断线段ab和cd是否相交
    return (sign(cross(b-a,c-a))*sign(cross(b-a,d-a))<0 && sign(cross(d-c,a-c))*sign(cross(d-c,b-c))<0) || in_segment(c,a,b) || in_segment(d,a,b) || in_segment(a,c,d) || in_segment(b,c,d);
} // 两个端点分别在另一条线段两侧，或者存在端点在另一条线段上
bool is_cross(Line a, Point A, Point B){ // 射线和线段是否相交(未验证)
    bool point_in_line = (eq(cross(A-a.u,a.dir),0) && sign((A-a.u)*a.dir)>=0) || (eq(cross(B-a.u,a.dir),0) && sign((B-a.u)*a.dir)>=0); // 线段端点在射线上
    bool is_cross = sign(cross(a.v-a.u,A-a.u))*sign(cross(a.v-a.u,B-a.u))<0 && sign((cross_point(a,Line(A,B))-a.u)*a.dir)>=0; // 存在交点
    return point_in_line || is_cross;
}
ld distance(Point a, Point b, Point c, Point d){ // 线段ab和cd距离
    if(is_cross(a,b,c,d)) return 0;
    return min(min(distance(a,c,d),distance(b,c,d)),min(distance(c,a,b),distance(d,a,b)));
}
bool in_polygon(Point p, Point *pos, int n) { // 判断点是否在多边形内（射线法）
    bool in = false;
    for(int i = 1; i <= n; i++) {
        Point A = pos[i], B = pos[i%n+1];
        if(in_segment(p, A, B)) return true; // 点在边上
        // 通过判断两端点是否在射线两侧来确定（一上一下）
        if(((sign(A.y - p.y) > 0) && (sign(B.y - p.y) <= 0)) || ((sign(B.y - p.y) > 0) && (sign(A.y - p.y) <= 0))) { // 注意，A和B都只有一次取等,若穿过顶点，只会计数一次；若与顶点相切，会计数2或0次
            Point cro = cross_point(Line(p, p+Point(1,0)), Line(A, B));
            if(sign(cro.x - p.x) > 0) in = !in;
        }
    }
    
    return in;
}

struct Circle{
    Point o;
    ld r;
    Circle(Point o=Point(), ld r=0):o(o),r(r){}
};
ld distance(Line l, Circle c){ // 直线到圆距离
    return fabs(distance(c.o,l)-c.r);
}
bool in_circle(Circle c,Point p){
    return sign(distance(c.o,p)-c.r)<=0;
}

// 圆和直线的交点
ppp cross_point(Circle c, Line l) {
    Point pro = projection(c.o, l);
    double d = sqrt(c.r*c.r - (pro-c.o)*(pro-c.o)); // 半弦长
    double len = length(l.dir);
    if (eq(d,0)) return {pro, pro};
    return {pro + l.dir * (d/len), pro - l.dir * (d/len)};
}
ppp cross_point(Circle c1, Circle c2){ // 两圆交点
    Vector v = c2.o-c1.o;
    ld d = length(v);
    ld cosA = (c1.r*c1.r+d*d-c2.r*c2.r)/(2*c1.r*d);
    ld angleA = acos(cosA);
    Vector v1 = rotate(v,angleA)*c1.r/d;
    Vector v2 = rotate(v,-angleA)*c1.r/d;
    return make_pair(c1.o+v1,c1.o+v2);
}

ppp tangent_point(Circle c, Point p){ // 点到圆的切点
    ld d_o = distance(c.o,p);
    ld d_t = sqrt(d_o*d_o-c.r*c.r);
    ld angleA = asin(c.r/d_o);
    Vector v = c.o-p;
    Vector v1 = rotate(v,angleA)*d_t/d_o;
    Vector v2 = rotate(v,-angleA)*d_t/d_o;
    return make_pair(p+v1,p+v2);
}

void tangent_point(Circle c1, Circle c2, Point *res, int &num){ // 两圆的公切点(只返回在c1上的切点)
    Vector v = c2.o-c1.o;
    ld d = length(v);
    int s1 = sign(d-(c1.r+c2.r)), s2 = sign(d-abs(c1.r-c2.r)); // 两个临界值
    num = 0;
    if(s2<0) return; // 内含
    if(s2==0){ // 内切
        if(c1.r<c2.r) swap(c1,c2);
        res[++num] = c1.o+(c2.o-c1.o)*c1.r/d;
        return;
    }
    ld angleA = acos((c1.r-c2.r)/d);
    Vector v1 = rotate(v,angleA)*c1.r/d;
    Vector v2 = rotate(v,-angleA)*c1.r/d;
    res[++num] = c1.o+v1;
    res[++num] = c1.o+v2;
    if(s1<0) return; // 相交
    if(s1==0){ // 外切
        res[++num] = c1.o+(c2.o-c1.o)*c1.r/d;
        return;
    }
    ld angleB = acos((c1.r+c2.r)/d);
    Vector v3 = rotate(v,angleB)*c1.r/d;
    Vector v4 = rotate(v,-angleB)*c1.r/d;
    res[++num] = c1.o+v3;
    res[++num] = c1.o+v4;
    return;
}

Circle inscribed_circle(Point a, Point b, Point c){ // 三角形内切圆
    ld A = distance(b,c), B = distance(a,c), C = distance(a,b);
    ld p = (A+B+C)/2; // 半周长
    ld r = sqrt((p-A)*(p-B)*(p-C)/p); // 海伦公式
    Point o = (a*A+b*B+c*C)/(A+B+C); // 奔驰定理
    return Circle(o,r);
}

Circle circum_circle(Point a, Point b, Point c) { // 三角形外接圆
    Vector n1 = normal(b-a); // 两条垂线
    Vector n2 = normal(c-a);
    Point m1 = (a+b)/2; // 两个中点
    Point m2 = (a+c)/2;
    Point o = cross_point(Line(m1, m1+n1), Line(m2, m2+n2));
    return Circle(o, distance(o, a));
}

ld CulArea(Point A, Point B, Circle C) { // 圆 与 圆心与两点连线 的面积交
    Vector OA=A-C.o, OB=B-C.o;
    Vector BA=A-B, BC=C.o-B;
    Vector AB=B-A, AC=C.o-A;
    ld DOA=length(OA), DOB=length(OB), DAB=length(AB), r=C.r;
    if(sign(cross(OA, OB))==0) return 0;
    if(sign(DOA-C.r)<0 && sign(DOB-C.r)<0) return cross(OA, OB)*0.5;
    else if(DOB<r && DOA>=r) {
        ld sqrt_arg=r*r*DAB*DAB-cross(BA, BC)*cross(BA, BC);
        if(sqrt_arg<0) sqrt_arg=0;
        ld x=(BA*BC+sqrt(sqrt_arg))/DAB;
        ld TS=cross(OA, OB)*0.5;
        return asin(TS*(1-x/DAB)*2/r/DOA)*r*r*0.5+TS*x/DAB;
    } 
    else if(DOB>=r && DOA<r) {
        ld sqrt_arg=r*r*DAB*DAB-cross(AB, AC)*cross(AB, AC);
        if(sqrt_arg<0) sqrt_arg=0;
        ld y=(AB*AC+sqrt(sqrt_arg))/DAB;
        ld TS=cross(OA, OB)*0.5;
        return asin(TS*(1-y/DAB)*2/r/DOB)*r*r*0.5+TS*y/DAB;
    } 
    else if(fabs(cross(OA, OB))>=r*DAB || AB*AC<=0 || BA*BC<=0){
        ld sin_arg=cross(OA, OB)/DOA/DOB;
        if(sin_arg>1) sin_arg=1;
        if(sin_arg<-1) sin_arg=-1;
        ld ang=asin(sin_arg);
        if(OA*OB<0) {
            if(cross(OA, OB)<0) return (-pi-ang)*r*r*0.5;
            else return (pi-ang)*r*r*0.5;
        } else return ang*r*r*0.5;
    } 
    else{
        ld sqrt_arg1=r*r*DAB*DAB-cross(BA, BC)*cross(BA, BC);
        ld sqrt_arg2=r*r*DAB*DAB-cross(AB, AC)*cross(AB, AC);
        if(sqrt_arg1<0) sqrt_arg1=0;
        if(sqrt_arg2<0) sqrt_arg2=0;
        ld x=(BA*BC+sqrt(sqrt_arg1))/DAB;
        ld y=(AB*AC+sqrt(sqrt_arg2))/DAB;
        ld TS=cross(OA, OB)*0.5;
        return (asin(TS*(1-x/DAB)*2/r/DOA)+asin(TS*(1-y/DAB)*2/r/DOB))*r*r*0.5+TS*((x+y)/DAB-1);
    }
}
ld inter_area(Circle c, Point *pos, int n){ // 圆与多边形的面积交
    ld res = 0;
    for(int i=1;i<=n;i++){
        res += CulArea(pos[i],pos[i%n+1],c);
    }
    return abs(res);
}

ld inter_area(Circle c1, Circle c2){ // 两圆面积交
    ld d = distance(c1.o,c2.o);
    if(sign(d-c1.r-c2.r)>=0) return 0; // 相离
    if(sign(d-abs(c1.r-c2.r))<=0) return pi*min(c1.r,c2.r)*min(c1.r,c2.r); // 内含
    ld angleA = acos((c1.r*c1.r+d*d-c2.r*c2.r)/(2*c1.r*d));
    ld angleB = acos((c2.r*c2.r+d*d-c1.r*c1.r)/(2*c2.r*d));
    return angleA*c1.r*c1.r+angleB*c2.r*c2.r-d*c1.r*sin(angleA);
}

std::ostream& operator<<(std::ostream& os, const Point& p) {
    os << "(" << p.x << ", " << p.y << ")";
    return os;
}
std::ostream& operator<<(std::ostream& os, const Line& l) {
    os << "Line(u: " << l.u << ", v: " << l.v << ", dir: " << l.dir << ")";
    return os;
}
std::ostream& operator<<(std::ostream& os, const Circle& c) {
    os << "Circle(o: " << c.o << ", r: " << c.r << ")";
    return os;
}

void convex(Point *pos, Point *res, int n, int &num){ // 求凸包
    num = 0;
    sort(pos+1,pos+n+1);
    for(int i=1;i<=n;i++){ // 下凸壳，逆时针转
        while(num>1 && sign(cross(res[num]-res[num-1],pos[i]-res[num-1]))<0) num--; // 改为<=0可以包含共线点
        res[++num] = pos[i];
    }
    int k = num;
    for(int i=n-1;i>=1;i--){ // 上凸壳，逆时针转
        while(num>k && sign(cross(res[num]-res[num-1],pos[i]-res[num-1]))<0) num--;
        res[++num] = pos[i];
    }
    if(n>1) num--;
}

ld diameter(Point *pos, int n){ // 求凸多边形直径
    int l = 1, r = 1;
    ld res = 0;
    for(int i=1;i<=n;i++){
        if(cmpy(pos[i],pos[l])) l = i;
        if(cmpy(pos[r],pos[i])) r = i;
    } // 一组对踵点
    res = lmax(res,distance(pos[l],pos[r]));
    for(int i=l,j=r;i!=r||j!=l;){ //没有回到起点，继续
        if(sign(cross(pos[i%n+1]-pos[i],pos[j%n+1]-pos[j]))<0) i = i%n+1; // 若能使叉积变大，继续
        else j = j%n+1;
        res = lmax(res,distance(pos[i],pos[j]));
    }
    return res;
}

ppp nps_re(Point *pos, int l,int r, ld &mindist) { // 递归部分
    int num = r-l+1;
    if(num<=1) return make_pair(Point(), Point(INFINITY,INFINITY));
    if(num==2) return make_pair(pos[l], pos[r]);

    int mid = (l+r)/2;
    ppp res = make_pair(Point(), Point());
    ppp resl = nps_re(pos,l,mid,mindist);
    ppp resr = nps_re(pos,mid+1,r,mindist);
    if(sign(distance(resl)-distance(resr))<0) res = resl;
    else res = resr;
    
    Point midp = pos[mid];
    Point *tmp = new Point[num+5];
    int numt = 0;
    for(int i=l;i<=r;i++){
        if(sign(abs(pos[i].x-midp.x)-mindist)<0) tmp[++numt] = pos[i];
    }
    sort(tmp+1,tmp+numt+1,cmpy);
    for(int i=1;i<=numt;i++){
        for(int j=i+1;j<=numt && tmp[j].y-tmp[i].y<mindist;j++){
            if(sign(distance(tmp[i],tmp[j])-mindist)<0){
                mindist = distance(tmp[i],tmp[j]);
                res = make_pair(tmp[i],tmp[j]);
            }
        }
    }

    delete[] tmp;
    return res;  // 返回最近点对
}
ppp nearest_points(Point *pos, int n){ // 最近点对
    sort(pos+1,pos+n+1); // 默认按x排序
    ld mindist = INFINITY;
    return nps_re(pos,1,n,mindist);
}

#define BOTTOM 0
#define LEFT 1
#define RIGHT 2
#define TOP 3
struct EndPoint{
    Point p;
    int id, type; //下左右上
    EndPoint(Point p=Point(), int id=0, int type=0):p(p),id(id),type(type){}

    bool operator < (const EndPoint &ep) const {
        //按y坐标升序排序
        if(p.y == ep.p.y) {
            return type < ep.type;//y相同时，按照下端点、左端点、右端点、上端点的顺序排列 
        } else return p.y < ep.p.y; 
	}
};
//线段相交问题：曼哈顿几何 EP开2倍大小
int manhattanIntersection(Line *S, EndPoint *EP, int n) {
	for(int i=1, k=0;i<=n;i++) {
		//调整端点p1、p2，保证左小右大
		if(S[i].u.y == S[i].v.y) {
			if(S[i].u.x > S[i].v.x) swap(S[i].u, S[i].v);
		} else if (S[i].u.y > S[i].v.y) swap(S[i].u, S[i].v);
		
		if(S[i].u.y == S[i].v.y) {
			EP[++k] = EndPoint(S[i].u, i, LEFT);
			EP[++k] = EndPoint(S[i].v, i, RIGHT);
		} else {
			EP[++k] = EndPoint(S[i].u, i, BOTTOM);
			EP[++k] = EndPoint(S[i].v, i, TOP);
		}
	}
	
	sort(EP+1,EP+2*n+1);//按端点的y坐标升序排列
	
	set<int> BT;//二叉搜索树
	BT.insert((int)1e9+10);//设置标记
	int cnt = 0;
	
	for(int i=1;i<=2*n;i++) {
		if(EP[i].type == TOP) {
			BT.erase(EP[i].p.x);//删除上端点 
		} else if(EP[i].type == BOTTOM) {
			BT.insert(EP[i].p.x);//添加下端点 
		} else if(EP[i].type == LEFT) {
			set<int>::iterator b = BT.lower_bound(S[EP[i].id].u.x);//O(log n)
			set<int>::iterator e = BT.upper_bound(S[EP[i].id].v.x);//O(log n)
			cnt += distance(b, e);//加上b和e的距离（点数） O(k) 
		}
	} 
	return cnt;
}

const int N = 1e5+10;
int n,num,q,r;
Point pos[N],res[N];
Circle c[5];
ld ans=0;
int main(){
    for(int i=1;i<=2;i++){
        cin>>c[i].o.x>>c[i].o.y>>c[i].r;
    }
    printf("%.10Lf\n",inter_area(c[1],c[2]));
    return 0;
}