#include <complex>
#include <math.h>
#include <fftw3.h>
#include <iostream>
#include <stdio.h>
using namespace std;

typedef complex<double> c;

class NotInAr{};
class TooBig{};
class TooSmall{};
class Ar{
public:
  int w,h;
  c*data; // row major order (id1,id2,id3,...) -> idn * in (id(n-1) * i(n-1) *(...
  Ar(char*filename){
    FILE*f=fopen(filename,"r");
    fscanf(f,"P5\n%d %d\n255\n",&w,&h);
    //if(w<10 || h<10)
    //  throw(TooSmall());
    //    if(w>800 || h>600)
    //  throw(TooBig());
    unsigned char buf[w*h];
    fread(buf,w,h,f);
    data=new c[w*h];
    for(int j=0;j<h;j++)
      for(int i=0;i<w;i++)
	(*this)(i,j)=buf[i+w*j];
    fclose(f);
  }
  Ar(int W,int H):w(W),h(H){
    data=new c[w*h];
  }
  // copy fouriertransform with correct zero padding
  Ar(Ar&a,int w,int h):w(w),h(h){
    //   if(a.w>w || a.h>h)
    //  throw NotInAr();
    int 
      W=min(a.w,w),
      H=min(a.h,h),
      DX=max(0,a.w-w),
      DY=max(0,a.h-h);
    data=new c[w*h];
    int i,j,dx=1,dy=1;
    // upper left
    for(j=0;j<H/2+dy;j++)
      for(i=0;i<W/2+dx;i++)
	(*this)(i,j)=a(i,j);
    // upper right
    for(j=0;j<H/2+dy;j++)
      for(i=W/2+dx;i<W;i++)
	(*this)(i-W+w,j)=a(i+DX,j);
    // lower left
    for(j=H/2+dy;j<H;j++)
      for(i=0;i<W/2+dx;i++)
	(*this)(i,j-H+h)=a(i,j+DY);
    // lower right
    for(j=H/2+dy;j<H;j++)
      for(i=W/2+dx;i<W;i++)
	(*this)(i-W+w,j-H+h)=a(i+DX,j+DY);
  }
  ~Ar(){
    //    delete [] data;
  }
  c&operator()(int i,int j){ // row major access
    if(i>=w || j>=h)
      throw NotInAr();
    return data[j+h*i];
  }
  c&operator[](int i){
    if(i<0 || i>w*h)
      throw(NotInAr());
    return data[i];
  }
  void output(char*filename){
    FILE*f=fopen(filename,"w");
    fprintf(f,"P5\n%d %d\n255\n",w,h);
    unsigned char buf[w*h]; // fill buf in column-major format
    for(int j=0;j<h;j++)
      for(int i=0;i<w;i++)
	buf[i+w*j]=(unsigned char)real((*this)(i,j));
    fwrite(buf,w,h,f);
    fclose(f);
  }
};

class FT{
public:
  Ar in,out;
  int dir;
  fftw_plan plan;
  
  FT(Ar&In,Ar&Out,int dir=FFTW_FORWARD):in(In),out(Out),dir(dir){
    if(out.w!=in.w || out.h!=in.h)
      throw NotInAr();
    plan=fftw_plan_dft_2d(in.w,in.h,
			  (fftw_complex*)in.data,
			  (fftw_complex*)out.data,
-			  dir,FFTW_ESTIMATE);
   }
  void operator()(){
    fftw_execute(plan);
  }
};

int
main(int argc,char**argv)
{
  if(argc!=3){
    printf("upscale a small image to screensize with fouriertransform.\n"
	   "usage: %s input.pgm output.pgm\n",argv[0]);
    return 0;
  }
  Ar in(argv[1]);
  int w=795, h=455,nw,nh;
  double 
    qw=w*1./in.w,
    qh=h*1./in.h;
  if(qh > qw){
    nw=int(qw*in.w+.5);
    nh=int(qw*in.h+.5);
  } else {
    nw=int(qh*in.w+.5);
    nh=int(qh*in.h+.5);
  }
  Ar kin(in.w,in.h);
  FT f(in,kin);
  
  f(); 
  Ar kbig(kin,nw,nh);
  Ar out(nw,nh);
  FT f2(kbig,out,FFTW_BACKWARD);
  f2();
  for(int i=0;i<nw*nh;i++)
    out[i]=round(real(out[i])*.75/(in.w*in.h)+42);
  out.output(argv[2]);
  return 0;
}

