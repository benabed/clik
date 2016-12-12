import itertools
import sys
import time

import scipy.integrate

from numpy import *
from numpy.linalg import *
from pylab import *
import slik_fast as sf


################################################################################
# class to hold stuff                                                          #
################################################################################

class stevenlike:
  def __init__(self,lmin=2,lmax=12,loadroot="teststore/",
               do_linear=True,use_offset_k2=False,
               lminlike=2,lmaxlike=10,regcl=True):
    if(lminlike<lmin):
      print 'error: lminlike less than lmin'
      assert(False)

    if(lmaxlike>lmax):
      print 'error: lmaxlike greater than lmax'
      assert(False)

    self.lmin = lmin
    self.lmax = lmax
    self.loadroot = loadroot
    self.do_linear=do_linear
    self.use_offset_k2=use_offset_k2
    if((not self.do_linear) & self.use_offset_k2):
      print 'Offset covariance only implemented for linear approximation'
      assert(False)
    self.lminlike=lminlike
    self.lmaxlike=lmaxlike

    self.regcleeval=0.
    self.regclbbval=0.
    if regcl:
      self.regcleeval=.1
      self.regclbbval=5.e-5

    relbase='lmin'+str(lmin)+'_lmax'+str(lmax)+'_calclike_07_Nov_2016.dat'
    self.yn_11=inpcovmat(loadroot+'yn_11'+relbase)
    self.yn_12=inpcovmat(loadroot+'yn_12'+relbase)
    self.yn_22=inpcovmat(loadroot+'yn_22'+relbase)
    self.yn_21=transpose(self.yn_12)
    self.f_11=inpcovmat(loadroot+'f_11'+relbase)
    self.f_22=inpcovmat(loadroot+'f_22'+relbase)
    self.ym_1=inpvec(loadroot+'ym_1'+relbase)
    self.ym_2=inpvec(loadroot+'ym_2'+relbase)
    y_1=inpvec(loadroot+'y_1'+relbase)
    y_2=inpvec(loadroot+'y_2'+relbase)
    self.b_l=inpvec(loadroot+'b_l'+relbase)

    self.yp_11=sum_vv(y_1,y_1,self)
    self.yp_22=sum_vv(y_2,y_2,self)
    self.yp_12=sum_vv(y_1,y_2,self)

    self.n_11_p=sum_mat(self.yn_11,self)
    self.n_22_p=sum_mat(self.yn_22,self)
    self.n_12_p=sum_mat(self.yn_12,self)

    #change that at some point
    clroot=loadroot
    fidcldatfile=clroot+'tau_scan_0_060_lensedCls.dat'
    self.ci=inv(d_symmean(self.f_11,self.f_22,self))
    self.fidcls=getclsmulti(fidcldatfile,(1,4,2),self)
    self.fidgrads_part= self.grad_s_part(self.ci,self.fidcls)

  def lkl(self,cls):
    #print >>sys.stderr,cls
    lkl=self.s_romb_part(cls)
    #print >>sys.stderr,"WOWOWO",lkl
    return lkl

  def __call__(self,cls):
    return self.lkl(cls)

  def s_romb_part(self,cls):
    lmin,lmax = self.lmin,self.lmax
    grads=zeros(5)
    dcls=(cls-self.fidcls)
    partdcls=dcls[3*(self.lminlike-lmin):3*(self.lmaxlike-lmin+1)].T
    grads[0]=dot(partdcls,self.fidgrads_part)
    grads[1]=dot(partdcls,self.grad_s_part(self.ci,self.fidcls+.25*dcls))
    grads[2]=dot(partdcls,self.grad_s_part(self.ci,self.fidcls+.50*dcls))
    grads[3]=dot(partdcls,self.grad_s_part(self.ci,self.fidcls+.75*dcls))
    grads[4]=dot(partdcls,self.grad_s_part(self.ci,cls))
    return scipy.integrate.romb(grads,.25,show=False)

  def grad_s_part(self,ci,cls):
    lmin,lmax,b_l = self.lmin,self.lmax,self.b_l
    
    #print 'generating pixel covmats...',
    cllong=makecllong_from_cls(cls,self)

    cov_11=make_m_slow(self.f_11, cllong, self.f_11, self.yn_11)
    cov_22=make_m_slow(self.f_22, cllong, self.f_22, self.yn_22)
    cov_12=make_m_slow(self.f_11, cllong, self.f_22, self.yn_12)

    nells=self.lmaxlike-self.lminlike+1

    nmin=3*(self.lminlike-self.lmin)
    nmax=3*(self.lmaxlike-self.lmin+1)
    
    pars=cls[nmin:nmax].copy()

    n1=shape(pars)[0]
    n2=n1*(n1+1)/2
    n=n1+n2

    #print 'generating cumulants...',
    #sys.stdout.flush()

    k1=pars
    
    if self.use_offset_k2:
      k2=dec_kap2(ci,symcov12_inc_offset(cov_11,cov_12,cov_22,ym_1,ym_2))[nmin:nmax,nmin:nmax]
    else:
      k2=dec_kap2(ci,symcov12_short(cov_11,cov_12,cov_22,self))[nmin:nmax,nmin:nmax]
    
    if self.do_linear:
        tmpdat=dec_kap1(ci,squash_vec(self.yp_12-self.n_12_p,lmin,lmax))[nmin:nmax]
        return -dot(dot(transpose(tmpdat-k1),inv(k2)),identity(n1))    
    
    # now that's the non linear part...

    #print "before k3"
    #t1 = time.time()
    k3=dec_kap3(ci,symkappa3_12(cov_11,cov_12,cov_22,self))[nmin:nmax,nmin:nmax,nmin:nmax].copy()
    #print "after k3",time.time()-t1
    
    #print "before k4"
    #t2 = time.time()
    udk4=symkappa4_12(cov_11,cov_12,cov_22,self)
    #print 'deconvolving k4...',  time.time()-t2
    
    #t3 = time.time()
    k4=dec_kap4(ci,udk4)[nmin:nmax,nmin:nmax,nmin:nmax,nmin:nmax].copy()
    #print 'done deconvolving k4...',time.time()-t3

    dk1=identity(n1)
    
    #print 'generating d_k2...',
    #t4 = time.time()
    uddk2=d_symcov12(cov_11,cov_12,cov_22,self.f_11,self.f_22,self)
    #print 'deconvolving d_k2...',  time.time()-t4
    
    #t5 = time.time()
    dk2=dec_derivkap2(ci,uddk2)[nmin:nmax,nmin:nmax,nmin:nmax].copy()
    #print 'done deconvolving d_k2...', time.time()-t5
    
    deriv=zeros((n,n1))
    diff=zeros(n)
    mat=zeros((n,n))

    #print 'filling matrix...',
    #t6 = time.time()

    mat[:n1,:n1]=k2

    p=n1
    for p1 in xrange(0,n1):
        for p2 in xrange(p1,n1):
            for q in xrange(n1):
                mat[p,q] = k3[p1,p2,q]
            p+=1

    for p in xrange(0,n1):
        q=n1
        for q1 in xrange(0,n1):
            for q2 in xrange(q1,n1):
                mat[p,q] = (
                    + k3[p,q1,q2]
                    + k2[p,q1]*k1[q2] + k2[p,q2]*k1[q1]
                    )
                q+=1

    
    p=n1
    for p1 in xrange(0,n1):
        for p2 in xrange(p1,n1):
            q=n1
            for q1 in xrange(0,n1):
                for q2 in xrange(q1,n1):
                    mat[p,q] = (
                        + k4[p1,p2,q1,q2]
                        + k3[p1,p2,q1]*k1[q2] + k3[p1,p2,q2]*k1[q1]
                        + k2[p1,q1]*k2[q2,p2] + k2[p1,q2]*k2[q1,p2]
                        )
                    q+=1
            p+=1

    
    #print time.time()-t6
    #print 'inverting matrix...',
    #t7 = time.time()
    matinv=inv(mat)
    #print time.time()-t7

    
    #print 'filling diff and deriv...',
    #t8=time.time()
    tk1=dec_kap1(ci,squash_vec(self.yp_12-self.n_12_p,lmin,lmax))[nmin:nmax].copy()

    for p in xrange(0,n1):
        diff[p]=tk1[p]-pars[p]
    p=n1
    for p1 in xrange(0,n1):
        for p2 in xrange(p1,n1):
            diff[p]=tk1[p1]*tk1[p2]-k2[p1,p2]-k1[p1]*k1[p2]
            p+=1

    deriv[:n1,:n1]=dk1
    p=n1
    for p1 in xrange(0,n1):
        for p2 in xrange(p1,n1):
            deriv[p,:]=dk2[p1,p2,:]
            p+=1

    #print 'computing likelihood...', time.time()-t8
    #print "total",time.time()-t1
    res = -dot(dot(transpose(diff),matinv),deriv)
    #print "RES",res
    return res

################################################################################
# reading files                                                                #
################################################################################

def inpvec(fname):
    return fromfile(fname,dtype=float32).astype(float)

def inpcovmat(fname):
    dat=fromfile(fname,dtype=float32).astype(float)
    n=int(sqrt(dat.shape))
    return reshape(dat,(n,n))

def getcls(filename, wantedcol,sih):
    lmin,lmax,b_l = sih.lmin,sih.lmax,sih.b_l
    inpdat=loadtxt(filename)
    ells=inpdat[lmin-2:lmax-1,0]
    dls=inpdat[lmin-2:lmax-1,wantedcol]
    return 2.0*math.pi*dls/ells/(ells+1.)*b_l[lmin:lmax+1]**2

#CHECK: Should this have the beam in or not?qqqq
def getclsmulti(filename, wantedcols,sih):
    lmin,lmax,b_l = sih.lmin,sih.lmax,sih.b_l
    inpdat=loadtxt(filename)
    ells=inpdat[lmin-2:lmax-1,0]
    ncols=len(wantedcols)
    cls=zeros((lmax+1-lmin)*ncols)
    wcol=0
    for col in wantedcols:
        cls[wcol::ncols]=2.0*math.pi*inpdat[lmin-2:lmax-1,col]/ells/(ells+1.)*b_l[lmin:lmax+1]**2
        wcol+=1
    return cls

def makecllong(filename,sih):
    lmin,lmax,b_l = sih.lmin,sih.lmax,sih.b_l
    inpdat=loadtxt(filename)
    inpells=inpdat[lmin-2:lmax-1,0]
    inpclee=inpdat[lmin-2:lmax-1,1] # i.e. TT
    inpclbb=inpdat[lmin-2:lmax-1,2] # i.e. EE
    inpcleb=inpdat[lmin-2:lmax-1,4]  # i.e. TE
    inpclbe=inpcleb
    outlen=2*((lmax+1)**2-lmin**2)
    cldiag=zeros((outlen,outlen))
    for l in range(lmin,lmax+1):
        for m in range(0,2*l+1):
            outpos=2*(l**2+m-lmin**2)
            cldiag[outpos,outpos]=2.0*math.pi*inpclee[l-lmin]/l/(l+1)*b_l[l]**2
            cldiag[outpos,outpos+1]=2.0*math.pi*inpcleb[l-lmin]/l/(l+1)*b_l[l]**2
            cldiag[outpos+1,outpos]=2.0*math.pi*inpclbe[l-lmin]/l/(l+1)*b_l[l]**2
            cldiag[outpos+1,outpos+1]=2.0*math.pi*inpclbb[l-lmin]/l/(l+1)*b_l[l]**2
    return cldiag

def makecllong_altlmax(filename,altlmax,sih):
    lmin,lmax,b_l = sih.lmin,sih.lmax,sih.b_l
    inpdat=loadtxt(filename)
    inpells=inpdat[lmin-2:altlmax-1,0]
    inpclee=inpdat[lmin-2:altlmax-1,1] # i.e. TT
    inpclbb=inpdat[lmin-2:altlmax-1,2] # i.e. EE
    inpcleb=inpdat[lmin-2:altlmax-1,4]  # i.e. TE
    inpclbe=inpcleb
    outlen=2*((altlmax+1)**2-lmin**2)
    cldiag=zeros((outlen,outlen))
    for l in range(lmin,altlmax+1):
        for m in range(0,2*l+1):
            outpos=2*(l**2+m-lmin**2)
            cldiag[outpos,outpos]=2.0*math.pi*inpclee[l-lmin]/l/(l+1)*b_l[l]**2
            cldiag[outpos,outpos+1]=2.0*math.pi*inpcleb[l-lmin]/l/(l+1)*b_l[l]**2
            cldiag[outpos+1,outpos]=2.0*math.pi*inpclbe[l-lmin]/l/(l+1)*b_l[l]**2
            cldiag[outpos+1,outpos+1]=2.0*math.pi*inpclbb[l-lmin]/l/(l+1)*b_l[l]**2
    return cldiag

def makecllong_regulate(filename,sih):
    lmin,lmax,b_l,regcleeval,regclbbval = sih.lmin,sih.lmax,sih.b_l,sih.regcleeval,sih.regclbbval
    inpdat=loadtxt(filename)
    inpells=inpdat[lmin-2:lmax-1,0]
    inpclee=inpdat[lmin-2:lmax-1,1] # i.e. TT
    inpclbb=inpdat[lmin-2:lmax-1,2] # i.e. EE
    inpcleb=inpdat[lmin-2:lmax-1,4]  # i.e. TE
    inpclbe=inpcleb
    outlen=2*((lmax+1)**2-lmin**2)
    cldiag=zeros((outlen,outlen))
    for l in range(lmin,lmax+1):
        for m in range(0,2*l+1):
            outpos=2*(l**2+m-lmin**2)
            cldiag[outpos,outpos]=2.0*math.pi*inpclee[l-lmin]/l/(l+1)*b_l[l]**2+regcleeval
            cldiag[outpos,outpos+1]=2.0*math.pi*inpcleb[l-lmin]/l/(l+1)*b_l[l]**2
            cldiag[outpos+1,outpos]=2.0*math.pi*inpclbe[l-lmin]/l/(l+1)*b_l[l]**2
            cldiag[outpos+1,outpos+1]=2.0*math.pi*inpclbb[l-lmin]/l/(l+1)*b_l[l]**2+regclbbval
    return cldiag
    
def makecllong_reshaped(filename,sih):
    lmin,lmax,b_l = sih.lmin,sih.lmax,sih.b_l
    inpdat=loadtxt(filename)
    inpells=inpdat[lmin-2:lmax-1,0]
    inpclee=inpdat[lmin-2:lmax-1,1] # i.e. TT
    inpclbb=inpdat[lmin-2:lmax-1,2] # i.e. EE
    inpcleb=zeros_like(inpclee) 
    inpclbe=zeros_like(inpclee)
    outlen=2*((lmax+1)**2-lmin**2)
    cldiag=zeros((outlen,outlen))
    for l in range(lmin,lmax+1):
        for m in range(0,2*l+1):
            outpos=2*(l**2+m-lmin**2)
            cldiag[outpos,outpos]=2.0*math.pi*inpclee[l-lmin]/l/(l+1)*b_l[l]**2+regcleeval
            cldiag[outpos,outpos+1]=2.0*math.pi*inpcleb[l-lmin]/l/(l+1)*b_l[l]**2
            cldiag[outpos+1,outpos]=2.0*math.pi*inpclbe[l-lmin]/l/(l+1)*b_l[l]**2
            cldiag[outpos+1,outpos+1]=2.0*math.pi*inpclbb[l-lmin]/l/(l+1)*b_l[l]**2+regclbbval
    return cldiag

def makecllong_reshaped_skip(filename,startl,sih):
    lmin,lmax,b_l = sih.lmin,sih.lmax,sih.b_l
    inpdat=loadtxt(filename)
    inpells=inpdat[lmin-2:lmax-1,0]
    inpclee=inpdat[lmin-2:lmax-1,1] # i.e. TT
    inpclbb=inpdat[lmin-2:lmax-1,2] # i.e. EE
    inpcleb=zeros_like(inpclee) 
    inpclbe=zeros_like(inpclee)
    outlen=2*((lmax+1)**2-lmin**2)
    cldiag=zeros((outlen,outlen))
    for l in range(startl,lmax+1):
        for m in range(0,2*l+1):
            outpos=2*(l**2+m-lmin**2)
            cldiag[outpos,outpos]=2.0*math.pi*inpclee[l-lmin]/l/(l+1)*b_l[l]**2+regcleeval
            cldiag[outpos,outpos+1]=2.0*math.pi*inpcleb[l-lmin]/l/(l+1)*b_l[l]**2
            cldiag[outpos+1,outpos]=2.0*math.pi*inpclbe[l-lmin]/l/(l+1)*b_l[l]**2
            cldiag[outpos+1,outpos+1]=2.0*math.pi*inpclbb[l-lmin]/l/(l+1)*b_l[l]**2+regclbbval
    return cldiag


def makecllong_from_cls(cls,sih):
    lmin,lmax,b_l = sih.lmin,sih.lmax,sih.b_l
    nells=lmax-lmin+1
    inpclee=cls[::3] # i.e. TT
    inpclbb=cls[2::3] # i.e. EE
    inpcleb=cls[1::3]  # i.e. TE
    inpclbe=inpcleb
    outlen=2*((lmax+1)**2-lmin**2)
    cldiag=zeros((outlen,outlen))
    for l in range(lmin,lmax+1):
        for m in range(0,2*l+1):
            outpos=2*(l**2+m-lmin**2)
            cldiag[outpos,outpos]=inpclee[l-lmin]*b_l[l]**2
            cldiag[outpos,outpos+1]=inpcleb[l-lmin]*b_l[l]**2
            cldiag[outpos+1,outpos]=inpclbe[l-lmin]*b_l[l]**2
            cldiag[outpos+1,outpos+1]=inpclbb[l-lmin]*b_l[l]**2
            
    return cldiag


################################################################################
# matrix utils                                                                 #
################################################################################

def mlnlike(d,m,cov):
    return 0.5*dot(transpose(d-m),dot(inv(cov),d-m))+0.5*slogdet(cov)[1]

def d_from_c(c,sih):
    ells=arange(sih.lmin,sih.lmax+1.0)
    return c*ells*(ells+1.0)/2.0/pi

def mul_rows(mat,d):
    mat2=mat.copy()
    for i in range(len(d)):
        mat2[i,:]*=d[i]
    return mat2

def mul_cols(mat,d):
    mat2=mat.copy()
    for i in range(len(d)):
        mat2[:,i]*=d[i]
    return mat2

def make_m_slow(mat1, clmat, mat2, mat3):
    return(dot(transpose(mat1),dot(clmat,mat2))+mat3)

def vecpart(p,v,lmin):
    pstart=p**2-lmin**2
    return v[2*pstart: 2*(pstart+2*p+1)].copy()

def sum_vv(v1,v2,sih):
    lmin,lmax = sih.lmin,sih.lmax
    means=zeros(4*(lmax-lmin+1))
    for l in range(lmin,lmax+1):
        s_ee=0.0
        s_eb=0.0
        s_be=0.0
        s_bb=0.0
        v1part=vecpart(l,v1,lmin)
        v2part=vecpart(l,v2,lmin)        
        for m in range(0,2*l+1):
            s_ee+=v1part[2*m]*v2part[2*m]
            s_eb+=v1part[2*m]*v2part[2*m+1]
            s_be+=v1part[2*m+1]*v2part[2*m]
            s_bb+=v1part[2*m+1]*v2part[2*m+1]
        means[4*(l-lmin)]=s_ee
        means[4*(l-lmin)+1]=s_eb
        means[4*(l-lmin)+2]=s_be
        means[4*(l-lmin)+3]=s_bb
    return means

def d_from_c_mat(m,sih):
    lmin,lmax = sih.lmin,sih.lmax
    ells=arange(lmin,lmax+1.0)
    scal=ells*(ells+1.0)/2.0/pi
    return mul_cols(mul_rows(m,scal),scal)

def matpart(p,q,mat,sih):
    lmin,lmax = sih.lmin,sih.lmax
    pstart=p**2-lmin**2
    qstart=q**2-lmin**2
    return(mat[2*pstart:2*(pstart+2*p+1),2*qstart:2*(qstart+2*q+1)])

def sum_mat(mat,sih):
    lmin,lmax = sih.lmin,sih.lmax
    means=zeros(4*(lmax-lmin+1))
    for l in range(lmin,lmax+1):
        s_ee=0.0
        s_eb=0.0
        s_be=0.0
        s_bb=0.0
        mpart=matpart(l,l,mat,sih)
        #        print mat.shape, mpart.shape
        for m in range(0,2*l+1):
            s_ee+=mpart[2*m,2*m]
            s_eb+=mpart[2*m,2*m+1]
            s_be+=mpart[2*m+1,2*m]
            s_bb+=mpart[2*m+1,2*m+1]
        means[4*(l-lmin)]=s_ee
        means[4*(l-lmin)+1]=s_eb
        means[4*(l-lmin)+2]=s_be
        means[4*(l-lmin)+3]=s_bb
    return means

def pairtrace(m1, m2):
   return sum(m1*(m2.transpose()))


def dec_kap1(ci,k1):
    n=k1.shape[0]
    k1t=k1.copy()
    dk1=zeros_like(k1)
    for i in xrange(n):
        for s in xrange(n):
            dk1[i]+=ci[i,s]*k1t[s]
    return dk1

def dec_kap2(ci,k2):
    n=k2.shape[0]

    k2t=k2.copy()
    dk2=zeros_like(k2)
    for p in range(n):
        for j in range(n):
            for s in range(n):
                dk2[p,j]+=ci[j,s]*k2t[p,s]

    k2t=dk2.copy()
    dk2=zeros_like(k2)
    for i in range(n):
        for j in range(n):
            for s in range(n):
                dk2[i,j]+=ci[i,s]*k2t[s,j]

    return dk2

dec_derivkap2 = sf.dec_derivkap2
def _dec_derivkap2(ci,k2):
    n=k2.shape[0]

    k2t=k2.copy()
    dk2=zeros_like(k2)
    for p in range(n):
        for j in range(n):
            for x in range(n):
                for s in range(n):
                    dk2[p,j,x]+=ci[j,s]*k2t[p,s,x]

    k2t=dk2.copy()
    dk2=zeros_like(k2)
    for i in range(n):
        for j in range(n):
            for x in range(n):
                for s in range(n):
                    dk2[i,j,x]+=ci[i,s]*k2t[s,j,x]

    return dk2

dec_kap3 = sf.dec_kap3
def _dec_kap3(ci,k3):
    n=k3.shape[0]

    k3t=k3.copy()
    dk3=zeros_like(k3)
    for p in range(n):
        for q in range(n):
            for k in range(n):
                for s in range(n):
                    dk3[p,q,k]+=ci[k,s]*k3t[p,q,s]
    
    k3t=dk3.copy()
    dk3=zeros_like(k3)
    for p in range(n):
        for j in range(n):
            for k in range(n):
                for s in range(n):
                    dk3[p,j,k]+=ci[j,s]*k3t[p,s,k]
    
    k3t=dk3.copy()
    dk3=zeros_like(k3)
    for i in range(n):
        for j in range(n):
            for k in range(n):
                for s in range(n):
                    dk3[i,j,k]+=ci[i,s]*k3t[s,j,k]

    return dk3

def dec_kap4(ci,k4):
    dk4=k4.copy()
    dk4=tensordot(ci,dk4,axes=([1],[3]))
    dk4=transpose(dk4,(1,2,3,0))
    dk4=tensordot(ci,dk4,axes=([1],[2]))
    dk4=transpose(dk4,(1,2,0,3))
    dk4=tensordot(ci,dk4,axes=([1],[1]))
    dk4=transpose(dk4,(1,0,2,3))
    dk4=tensordot(ci,dk4,axes=([1],[0]))
    return dk4

def squash_vec(v,lmin,lmax):
    y=zeros(3*(lmax-lmin+1))
    for p in range(0,(lmax-lmin+1)):
        y[3*p]=v[4*p]
        y[3*p+1]=v[4*p+1]+v[4*p+2]
        y[3*p+2]=v[4*p+3]
    return y


def grad_s_part(ci,cls,sih):
    lmin,lmax,b_l = sih.lmin,sih.lmax,sih.b_l
    
    #print 'generating pixel covmats...',
    cllong=makecllong_from_cls(cls,sih)

    cov_11=make_m_slow(sih.f_11, cllong, sih.f_11, sih.yn_11)
    cov_22=make_m_slow(sih.f_22, cllong, sih.f_22, sih.yn_22)
    cov_12=make_m_slow(sih.f_11, cllong, sih.f_22, sih.yn_12)

    nells=sih.lmaxlike-sih.lminlike+1

    nmin=3*(sih.lminlike-sih.lmin)
    nmax=3*(sih.lmaxlike-sih.lmin+1)
    
    pars=cls[nmin:nmax].copy()

    n1=shape(pars)[0]
    n2=n1*(n1+1)/2
    n=n1+n2

    #print 'generating cumulants...',
    #sys.stdout.flush()

    k1=pars.copy()
    
    t1 = time.time()
    k2=dec_kap2(ci,symcov12_short(cov_11,cov_12,cov_22,sih))[nmin:nmax,nmin:nmax].copy()
    print "linear case",time.time()-t1 
    if sih.do_linear:
        if sih.use_offset_k2:
            k2=dec_kap2(ci,symcov12_inc_offset(cov_11,cov_12,cov_22,ym_1,ym_2))[nmin:nmax,nmin:nmax].copy()
    # testing returning the "linear" result...
        tmpdat=dec_kap1(ci,squash_vec(sih.yp_12-sih.n_12_p,lmin,lmax))[nmin:nmax].copy()
    #tmpdat=fidcls
        return -dot(dot(transpose(tmpdat-k1),inv(k2)),identity(n1))    
    #tmpdiff=zeros_like(cls)
    #tmpdiff[2::3]=(ells+.5)*(-fidcls[2::3]/k1[2::3]**2+1./k1[2::3])
    #return tmpdiff

    #testing zeroing the higher cumulants...
    print "before k3"
    t1 = time.time()
    k3=dec_kap3(ci,symkappa3_12(cov_11,cov_12,cov_22,sih))[nmin:nmax,nmin:nmax,nmin:nmax].copy()
    print "after k3",time.time()-t1
    #udk4=symkappa4_12(cov_11,cov_12,cov_22)
    #k4=dec_kap4(ci,udk4)
    #sys.exit()
    print "before k4"
    t2 = time.time()
    udk4=symkappa4_12(cov_11,cov_12,cov_22,sih)
    print 'deconvolving k4...',  time.time()-t2
    sys.stdout.flush()
    t3 = time.time()
    k4=dec_kap4(ci,udk4)[nmin:nmax,nmin:nmax,nmin:nmax,nmin:nmax].copy()
    print 'done deconvolving k4...',time.time()-t3
    sys.stdout.flush()
    #k3=zeros((n1,n1,n1))
    #k4=zeros((n1,n1,n1,n1))

    #    print k1
    #print k2
    #print k3
    #print k4


    dk1=identity(n1)
    print 'generating d_k2...',
    sys.stdout.flush()
    t4 = time.time()
    uddk2=d_symcov12(cov_11,cov_12,cov_22,sih.f_11,sih.f_22,sih)
    print 'deconvolving d_k2...',  time.time()-t4
    sys.stdout.flush()
    t5 = time.time()
    dk2=dec_derivkap2(ci,uddk2)[nmin:nmax,nmin:nmax,nmin:nmax].copy()
    print 'done deconvolving d_k2...', time.time()-t5
    sys.stdout.flush()

    #print dk1
    #print dk2


    deriv=zeros((n,n1))
    diff=zeros(n)
    mat=zeros((n,n))

    print 'filling matrix...',
    sys.stdout.flush()
    t6 = time.time()

    mat[:n1,:n1]=k2

    p=n1
    for p1 in xrange(0,n1):
        for p2 in xrange(p1,n1):
            for q in xrange(n1):
                mat[p,q] = k3[p1,p2,q]
            p+=1

    for p in xrange(0,n1):
        q=n1
        for q1 in xrange(0,n1):
            for q2 in xrange(q1,n1):
                mat[p,q] = (
                    + k3[p,q1,q2]
                    + k2[p,q1]*k1[q2] + k2[p,q2]*k1[q1]
                    )
                q+=1

    #zeroing k3 for xx-xx block...
    #k3=zeros((n1,n1,n1))            
    p=n1
    for p1 in xrange(0,n1):
        for p2 in xrange(p1,n1):
            q=n1
            for q1 in xrange(0,n1):
                for q2 in xrange(q1,n1):
                    mat[p,q] = (
                        + k4[p1,p2,q1,q2]
                        + k3[p1,p2,q1]*k1[q2] + k3[p1,p2,q2]*k1[q1]
                        + k2[p1,q1]*k2[q2,p2] + k2[p1,q2]*k2[q1,p2]
                        )
                    q+=1
            p+=1

            #print mat

    print time.time()-t6
    print 'inverting matrix...',
    sys.stdout.flush()
    t7 = time.time()
    matinv=inv(mat)
    print time.time()-t7

    #print matinv

    #print dot(mat,matinv)

    print 'filling diff and deriv...',
    sys.stdout.flush()
    t8=time.time()
    tk1=dec_kap1(ci,squash_vec(sih.yp_12-sih.n_12_p,lmin,lmax))[nmin:nmax].copy()

    for p in xrange(0,n1):
        diff[p]=tk1[p]-pars[p]
    p=n1
    for p1 in xrange(0,n1):
        for p2 in xrange(p1,n1):
            diff[p]=tk1[p1]*tk1[p2]-k2[p1,p2]-k1[p1]*k1[p2]
            p+=1

    deriv[:n1,:n1]=dk1
    p=n1
    for p1 in xrange(0,n1):
        for p2 in xrange(p1,n1):
            deriv[p,:]=dk2[p1,p2,:]
            p+=1

            #assert(False)

    print 'computing likelihood...', time.time()-t8
    print "total",time.time()-t1
    sys.stdout.flush()
    res = -dot(dot(transpose(diff),matinv),deriv)
    print "RES",res
    return res

# for now assuming M_AB^T=M_BA...
def cov_12_34(mat23,mat41,mat24,mat31,sih):
    lmin,lmax = sih.lmin,sih.lmax
    cov=zeros((4*(lmax-lmin+1),4*(lmax-lmin+1)))
    for p in range(lmin,lmax+1):
        for q in range(lmin,lmax+1):
            m23=matpart(p,q,mat23,sih)
            m41=matpart(q,p,mat41,sih)
            m24=matpart(p,q,mat24,sih)
            m31=matpart(q,p,mat31,sih)
            for m in range(0,2*p+1):
                for n in range(0,2*q+1):
                    cov[4*(p-lmin)+0,4*(q-lmin)+0]+=m23[2*m+0,2*n+0]*m41[2*n+0,2*m+0]+m24[2*m+0,2*n+0]*m31[2*n+0,2*m+0]
                    cov[4*(p-lmin)+0,4*(q-lmin)+1]+=m23[2*m+0,2*n+0]*m41[2*n+1,2*m+0]+m24[2*m+0,2*n+1]*m31[2*n+0,2*m+0]
                    cov[4*(p-lmin)+0,4*(q-lmin)+2]+=m23[2*m+0,2*n+1]*m41[2*n+0,2*m+0]+m24[2*m+0,2*n+0]*m31[2*n+1,2*m+0]
                    cov[4*(p-lmin)+0,4*(q-lmin)+3]+=m23[2*m+0,2*n+1]*m41[2*n+1,2*m+0]+m24[2*m+0,2*n+1]*m31[2*n+1,2*m+0]
                    cov[4*(p-lmin)+1,4*(q-lmin)+0]+=m23[2*m+1,2*n+0]*m41[2*n+0,2*m+0]+m24[2*m+1,2*n+0]*m31[2*n+0,2*m+0]
                    cov[4*(p-lmin)+1,4*(q-lmin)+1]+=m23[2*m+1,2*n+0]*m41[2*n+1,2*m+0]+m24[2*m+1,2*n+1]*m31[2*n+0,2*m+0]
                    cov[4*(p-lmin)+1,4*(q-lmin)+2]+=m23[2*m+1,2*n+1]*m41[2*n+0,2*m+0]+m24[2*m+1,2*n+0]*m31[2*n+1,2*m+0]
                    cov[4*(p-lmin)+1,4*(q-lmin)+3]+=m23[2*m+1,2*n+1]*m41[2*n+1,2*m+0]+m24[2*m+1,2*n+1]*m31[2*n+1,2*m+0]
                    cov[4*(p-lmin)+2,4*(q-lmin)+0]+=m23[2*m+0,2*n+0]*m41[2*n+0,2*m+1]+m24[2*m+0,2*n+0]*m31[2*n+0,2*m+1]
                    cov[4*(p-lmin)+2,4*(q-lmin)+1]+=m23[2*m+0,2*n+0]*m41[2*n+1,2*m+1]+m24[2*m+0,2*n+1]*m31[2*n+0,2*m+1]
                    cov[4*(p-lmin)+2,4*(q-lmin)+2]+=m23[2*m+0,2*n+1]*m41[2*n+0,2*m+1]+m24[2*m+0,2*n+0]*m31[2*n+1,2*m+1]
                    cov[4*(p-lmin)+2,4*(q-lmin)+3]+=m23[2*m+0,2*n+1]*m41[2*n+1,2*m+1]+m24[2*m+0,2*n+1]*m31[2*n+1,2*m+1]
                    cov[4*(p-lmin)+3,4*(q-lmin)+0]+=m23[2*m+1,2*n+0]*m41[2*n+0,2*m+1]+m24[2*m+1,2*n+0]*m31[2*n+0,2*m+1]
                    cov[4*(p-lmin)+3,4*(q-lmin)+1]+=m23[2*m+1,2*n+0]*m41[2*n+1,2*m+1]+m24[2*m+1,2*n+1]*m31[2*n+0,2*m+1]
                    cov[4*(p-lmin)+3,4*(q-lmin)+2]+=m23[2*m+1,2*n+1]*m41[2*n+0,2*m+1]+m24[2*m+1,2*n+0]*m31[2*n+1,2*m+1]
                    cov[4*(p-lmin)+3,4*(q-lmin)+3]+=m23[2*m+1,2*n+1]*m41[2*n+1,2*m+1]+m24[2*m+1,2*n+1]*m31[2*n+1,2*m+1]
    return cov       

symcov12_short = sf.symcov12_short

def _symcov12_short(mat11,mat12,mat22,sih):
    lmin,lmax = sih.lmin,sih.lmax
    mat21=transpose(mat12).copy()
    cov=zeros((3*(lmax-lmin+1),3*(lmax-lmin+1)))
    for l in itertools.product(xrange(lmin,lmax+1),repeat=2):
        zz1=zeros((3,3))
        m23=matpart(l[0],l[1],mat21,sih)
        m41=matpart(l[1],l[0],mat21,sih)
        m24=matpart(l[0],l[1],mat22,sih)
        m31=matpart(l[1],l[0],mat11,sih)
        for i in itertools.product([0,1],repeat=4):
           #print pairtrace(m23[i[1]::2,i[2]::2],m41[i[3]::2,i[0]::2]) 
           #print zz1[i[0]+i[1],i[2]+i[3]]
           zz1[i[0]+i[1],i[2]+i[3]] += (
                + pairtrace(m23[i[1]::2,i[2]::2],m41[i[3]::2,i[0]::2])
                + pairtrace(m24[i[1]::2,i[3]::2],m31[i[2]::2,i[0]::2])
                )
        cov[3*(l[0]-lmin):3*(l[0]+1-lmin),3*(l[1]-lmin):3*(l[1]+1-lmin)]=zz1
    return cov

# WARNING: Provisionally adding a b_l'^2 to the coupling matrix...
# note an arbitrary tranpose in i[2] and i[3] to match cov formula...
def d_symmean(mat11_fid,mat22_fid,sih):
    lmin,lmax = sih.lmin,sih.lmax
    coup=zeros((3*(lmax-lmin+1),3*(lmax-lmin+1)))
    for l in itertools.product(range(lmin,lmax+1),repeat=2):
        m22=matpart(l[0],l[1],mat22_fid,sih)
        m11=matpart(l[1],l[0],mat11_fid,sih)
        zz1=zeros((3,3))
        for m in itertools.product(range(0,2*l[0]+1),range(0,2*l[1]+1)):
            for i in itertools.product([0,1],repeat=4):
                zz1[i[0]+i[1],i[2]+i[3]]+=m22[2*m[0]+i[1],2*m[1]+i[3]]*m11[2*m[1]+i[2],2*m[0]+i[0]]
        coup[3*(l[0]-lmin):3*(l[0]+1-lmin),3*(l[1]-lmin):3*(l[1]+1-lmin)]=zz1*sih.b_l[l[1]]**2
    return coup

# warning: checking order of i[] terms...
def symcov12(mat11,mat12,mat22,sih):
    lmin,lmax = sih.lmin,sih.lmax
    mat21=transpose(mat12).copy()
    cov=zeros((3*(lmax-lmin+1),3*(lmax-lmin+1)))
    for l in itertools.product(xrange(lmin,lmax+1),repeat=2):
        zz1=zeros((3,3))
        m23=matpart(l[0],l[1],mat21,sih)
        m41=matpart(l[1],l[0],mat21,sih)
        m24=matpart(l[0],l[1],mat22,sih)
        m31=matpart(l[1],l[0],mat11,sih)
        for m in itertools.product(xrange(0,2*l[0]+1),xrange(0,2*l[1]+1)):
            for i in itertools.product([0,1],repeat=4):
                zz1[i[0]+i[1],i[2]+i[3]] += (
                    + m23[2*m[0]+i[1],2*m[1]+i[2]]*m41[2*m[1]+i[3],2*m[0]+i[0]]
                    + m24[2*m[0]+i[1],2*m[1]+i[3]]*m31[2*m[1]+i[2],2*m[0]+i[0]]
                    )
        cov[3*(l[0]-lmin):3*(l[0]+1-lmin),3*(l[1]-lmin):3*(l[1]+1-lmin)]=zz1
    return cov

def symcov12_inc_offset(mat11,mat12,mat22,ym1,ym2,sih):
    lmin,lmax = sih.lmin,sih.lmax
    mat21=transpose(mat12).copy()
    cov=zeros((3*(lmax-lmin+1),3*(lmax-lmin+1)))
    for l in itertools.product(xrange(lmin,lmax+1),repeat=2):
        zz1=zeros((3,3))
        v1=vecpart(l[0],ym1,lmin)
        v2=vecpart(l[0],ym2,lmin)
        v3=vecpart(l[1],ym1,lmin)
        v4=vecpart(l[1],ym2,lmin)
        m23=matpart(l[0],l[1],mat21,sih)+outer(v2,v3)
        m41=matpart(l[1],l[0],mat21,sih)+outer(v4,v1)
        m24=matpart(l[0],l[1],mat22,sih)+outer(v2,v4)
        m31=matpart(l[1],l[0],mat11,sih)+outer(v3,v1)
        for m in itertools.product(xrange(0,2*l[0]+1),xrange(0,2*l[1]+1)):
            for i in itertools.product([0,1],repeat=4):
                zz1[i[0]+i[1],i[2]+i[3]] += (
                    + m23[2*m[0]+i[1],2*m[1]+i[2]]*m41[2*m[1]+i[3],2*m[0]+i[0]]
                    + m24[2*m[0]+i[1],2*m[1]+i[3]]*m31[2*m[1]+i[2],2*m[0]+i[0]]
                    - 2.*v1[2*m[0]+i[0]]*v2[2*m[0]+i[1]]*v3[2*m[1]+i[2]]*v4[2*m[1]+i[3]]
                    )
        cov[3*(l[0]-lmin):3*(l[0]+1-lmin),3*(l[1]-lmin):3*(l[1]+1-lmin)]=zz1
    return cov

symkappa3_12 = sf.symkappa3_12


def _symkappa3_12(mat11,mat12,mat22,sih):
    lmin,lmax = sih.lmin,sih.lmax
    mat21=transpose(mat12).copy()
    k3=zeros((3*(lmax-lmin+1),3*(lmax-lmin+1),3*(lmax-lmin+1)))
    outlen=2*((lmax+1)**2-lmin**2)
    m_11_21=zeros((2,2,outlen,outlen))
    m_11_22=zeros((2,2,outlen,outlen))
    m_21_21=zeros((2,2,outlen,outlen))
    m_21_22=zeros((2,2,outlen,outlen))

    t1 = time.time()
    for p in range(lmin,lmax+1):
        m_11_21[0,0,:,:]=dot(mat11[:,2*(p**2-lmin**2):2*(2*p+1+p**2-lmin**2):2],mat21[2*(p**2-lmin**2):2*(2*p+1+p**2-lmin**2):2,:])
        m_11_21[0,1,:,:]=dot(mat11[:,2*(p**2-lmin**2):2*(2*p+1+p**2-lmin**2):2],mat21[2*(p**2-lmin**2)+1:2*(2*p+1+p**2-lmin**2):2,:])
        m_11_21[1,0,:,:]=dot(mat11[:,2*(p**2-lmin**2)+1:2*(2*p+1+p**2-lmin**2):2],mat21[2*(p**2-lmin**2):2*(2*p+1+p**2-lmin**2):2,:])
        m_11_21[1,1,:,:]=dot(mat11[:,2*(p**2-lmin**2)+1:2*(2*p+1+p**2-lmin**2):2],mat21[2*(p**2-lmin**2)+1:2*(2*p+1+p**2-lmin**2):2,:])

        m_11_22[0,0,:,:]=dot(mat11[:,2*(p**2-lmin**2):2*(2*p+1+p**2-lmin**2):2],mat22[2*(p**2-lmin**2):2*(2*p+1+p**2-lmin**2):2,:])
        m_11_22[0,1,:,:]=dot(mat11[:,2*(p**2-lmin**2):2*(2*p+1+p**2-lmin**2):2],mat22[2*(p**2-lmin**2)+1:2*(2*p+1+p**2-lmin**2):2,:])
        m_11_22[1,0,:,:]=dot(mat11[:,2*(p**2-lmin**2)+1:2*(2*p+1+p**2-lmin**2):2],mat22[2*(p**2-lmin**2):2*(2*p+1+p**2-lmin**2):2,:])
        m_11_22[1,1,:,:]=dot(mat11[:,2*(p**2-lmin**2)+1:2*(2*p+1+p**2-lmin**2):2],mat22[2*(p**2-lmin**2)+1:2*(2*p+1+p**2-lmin**2):2,:])

        m_21_21[0,0,:,:]=dot(mat21[:,2*(p**2-lmin**2):2*(2*p+1+p**2-lmin**2):2],mat21[2*(p**2-lmin**2):2*(2*p+1+p**2-lmin**2):2,:])
        m_21_21[0,1,:,:]=dot(mat21[:,2*(p**2-lmin**2):2*(2*p+1+p**2-lmin**2):2],mat21[2*(p**2-lmin**2)+1:2*(2*p+1+p**2-lmin**2):2,:])
        m_21_21[1,0,:,:]=dot(mat21[:,2*(p**2-lmin**2)+1:2*(2*p+1+p**2-lmin**2):2],mat21[2*(p**2-lmin**2):2*(2*p+1+p**2-lmin**2):2,:])
        m_21_21[1,1,:,:]=dot(mat21[:,2*(p**2-lmin**2)+1:2*(2*p+1+p**2-lmin**2):2],mat21[2*(p**2-lmin**2)+1:2*(2*p+1+p**2-lmin**2):2,:])

        m_21_22[0,0,:,:]=dot(mat21[:,2*(p**2-lmin**2):2*(2*p+1+p**2-lmin**2):2],mat22[2*(p**2-lmin**2):2*(2*p+1+p**2-lmin**2):2,:])
        m_21_22[0,1,:,:]=dot(mat21[:,2*(p**2-lmin**2):2*(2*p+1+p**2-lmin**2):2],mat22[2*(p**2-lmin**2)+1:2*(2*p+1+p**2-lmin**2):2,:])
        m_21_22[1,0,:,:]=dot(mat21[:,2*(p**2-lmin**2)+1:2*(2*p+1+p**2-lmin**2):2],mat22[2*(p**2-lmin**2):2*(2*p+1+p**2-lmin**2):2,:])
        m_21_22[1,1,:,:]=dot(mat21[:,2*(p**2-lmin**2)+1:2*(2*p+1+p**2-lmin**2):2],mat22[2*(p**2-lmin**2)+1:2*(2*p+1+p**2-lmin**2):2,:])

        for q in range(p,lmax+1):
            for r in range(q,lmax+1):
                z=zeros((3,3,3))
                for inds in itertools.product([0,1],repeat=6):
                    biginds=tuple(array(inds[::2])+array(inds[1::2]))
                    z[biginds] += (
                        + pairtrace(matpart(r,q,m_21_21[inds[0],inds[1],:,:],sih)[inds[5]::2,inds[2]::2],matpart(q,r,mat21,sih)[inds[3]::2,inds[4]::2])
                        + pairtrace(matpart(r,q,m_11_21[inds[0],inds[1],:,:],sih)[inds[4]::2,inds[2]::2],matpart(q,r,mat22,sih)[inds[3]::2,inds[5]::2])
                        + pairtrace(matpart(r,q,m_21_22[inds[0],inds[1],:,:],sih)[inds[5]::2,inds[3]::2],matpart(q,r,mat11,sih)[inds[2]::2,inds[4]::2])
                        + pairtrace(matpart(r,q,m_11_22[inds[0],inds[1],:,:],sih)[inds[4]::2,inds[3]::2],matpart(q,r,mat12,sih)[inds[2]::2,inds[5]::2])
                        + pairtrace(matpart(q,r,m_21_21[inds[0],inds[1],:,:],sih)[inds[3]::2,inds[4]::2],matpart(r,q,mat21,sih)[inds[5]::2,inds[2]::2])
                        + pairtrace(matpart(q,r,m_11_21[inds[0],inds[1],:,:],sih)[inds[2]::2,inds[4]::2],matpart(r,q,mat22,sih)[inds[5]::2,inds[3]::2])
                        + pairtrace(matpart(q,r,m_21_22[inds[0],inds[1],:,:],sih)[inds[3]::2,inds[5]::2],matpart(r,q,mat11,sih)[inds[4]::2,inds[2]::2])
                        + pairtrace(matpart(q,r,m_11_22[inds[0],inds[1],:,:],sih)[inds[2]::2,inds[5]::2],matpart(r,q,mat12,sih)[inds[4]::2,inds[3]::2])
                        )
                k3[3*(p-lmin):3*(p+1-lmin),3*(q-lmin):3*(q+1-lmin),3*(r-lmin):3*(r+1-lmin)]=z     
            #print p,q,r
    t2 = time.time()
    print "s",t2-t1
    for ells in (itertools.product(range(lmin,lmax+1),repeat=3)):
        p,q,r=list(ells)
        sargs=argsort(ells)
        psort,qsort,rsort=sorted(list(ells))
        zz1 = k3[3*(psort-lmin):3*(psort+1-lmin),3*(qsort-lmin):3*(qsort+1-lmin),3*(rsort-lmin):3*(rsort+1-lmin)].copy()
        zz2 = zeros_like(zz1)
        for inds in itertools.product(range(0,3),repeat=3):
            tinds=list(inds)
            sinds=array(tinds)[sargs]
            zz2[tuple(tinds)] = zz1[tuple(sinds)]
        k3[3*(p-lmin):3*(p+1-lmin),3*(q-lmin):3*(q+1-lmin),3*(r-lmin):3*(r+1-lmin)] = zz2
    print "s",time.time()-t2
    #k3.tofile("sk3_slow.dat")
    return k3

symkappa4_12 = sf.symkappa4_12
# Probably the slowest one in the slow version. OPTIMIZE THIS ONE FIRST
import sys,time
def _symkappa4_12(mat11,mat12,mat22,sih):
    lmin,lmax = sih.lmin,sih.lmax
    mat21=transpose(mat12).copy()
    mats=zeros((2,2,mat11.shape[0],mat11.shape[1]))
    mats[0,0,:,:]=mat11
    mats[0,1,:,:]=mat12
    mats[1,0,:,:]=mat21
    mats[1,1,:,:]=mat22

    count=0
    
    k4=zeros((3*(lmax-lmin+1),3*(lmax-lmin+1),3*(lmax-lmin+1),3*(lmax-lmin+1)))
    outlen=2*((lmax+1)**2-lmin**2)
    # l, inside freqs, inside t/e, outside freqs, outside t/e by ::2 or 1::2  
    m=zeros((lmax-lmin+1,2,2,2,2,2,2,outlen,outlen))
    for p in range(lmin,lmax+1):
        for i in itertools.product([0,1],repeat=6):
            m1=mats[i[4],i[0],:,i[2]+2*(p**2-lmin**2):2*(2*p+1+p**2-lmin**2):2].copy()
            m2=mats[i[1],i[5],i[3]+2*(p**2-lmin**2):2*(2*p+1+p**2-lmin**2):2,:].copy()
            m[p-lmin,i[0],i[1],i[2],i[3],i[4],i[5],:,:]=dot(m1,m2)
        #print p

    for p in xrange(lmin,lmax+1):
        for q in xrange(p,lmax+1):
            for r in xrange(q,lmax+1):
                for s in xrange(r,lmax+1):
                    #count+=1
                    #continue
                    _t1 = time.time()
                    z=zeros((3,3,3,3))
                        #                        count=0
                    for changingmatinds in itertools.permutations(zip([q,r,s],range(2,8,2),range(3,8,2))):
                        for perms in itertools.product([1,-1],repeat=3):
                            terms=[(p,0,1)]+list(changingmatinds)
                            fullperms=[1]+list(perms)
                            freqinds=list()
                            lterms=list()
                            for i in range(len(terms)):
                                freqinds+=list(array(terms[i][1:3])[::fullperms[i]])
                                #freqinds+=terms[i][1:3]
                                #print terms[i][0],array(terms[i][1:3])[::fullperms[i]],terms[i][0],' ',
                                lterms+=[terms[i][0]]
                                lterms+=[terms[i][0]]
                            fis=freqinds[-1:]+freqinds[:-1]
                            fls=lterms[-1:]+lterms[:-1]
                            for inds1 in itertools.product([0,1],repeat=4):
                                mp1full=matpart(fls[0],fls[3],m[fls[1]-lmin,fis[1]%2,fis[2]%2,inds1[0],inds1[1],fis[0]%2,fis[3]%2,:,:],sih)
                                mp2full=matpart(fls[4],fls[7],m[fls[5]-lmin,fis[5]%2,fis[6]%2,inds1[2],inds1[3],fis[4]%2,fis[7]%2,:,:],sih)
                                for inds2 in itertools.product([0,1],repeat=4):
                                #count+=1
                                #biginds=tuple(2*array(inds[::2])+array(inds[1::2]))
                                #print fls[0],fls[3],fls[1]-lmin,fis[1]%2,fis[2]%2,biginds[0],biginds[1],fis[0]%2,fis[3]%2,biginds[2],biginds[3]
                                    mp1=mp1full[inds2[0]::2,inds2[1]::2]
                                    mp2=mp2full[inds2[2]::2,inds2[3]::2]
                                    biginds=(inds1[0]+inds1[1],inds2[0]+inds2[1],inds1[2]+inds1[3],inds2[2]+inds2[3])
                                    zval = pairtrace(mp1,mp2)
                                    #print "slow",zval
                                    z[biginds] += zval
                                    
                                #count+=1

                    k4[3*(p-lmin):3*(p+1-lmin),3*(q-lmin):3*(q+1-lmin),3*(r-lmin):3*(r+1-lmin),3*(s-lmin):3*(s+1-lmin)]=z     
                    t1 = time.time() - _t1
                    #print "(",p,q,r,s,")",t1

                    
    #k4.tofile("slik_slow_k4.dat")
    for ells in (itertools.product(range(lmin,lmax+1),repeat=4)):
        p,q,r,s=list(ells)
        sargs=argsort(ells)
        psort,qsort,rsort,ssort=sorted(list(ells))
        zz1 = k4[3*(psort-lmin):3*(psort+1-lmin),3*(qsort-lmin):3*(qsort+1-lmin),3*(rsort-lmin):3*(rsort+1-lmin),3*(ssort-lmin):3*(ssort+1-lmin)].copy()
        zz2 = zeros_like(zz1)
        for inds in itertools.product(range(0,3),repeat=4):
            tinds=list(inds)
            #sinds=tinds
            sinds2=list()
            for i in sargs:
                sinds2.append(tinds[i])
            #sinds=array(tinds)[sargs]
            zz2[tuple(tinds)] = zz1[tuple(sinds2)]
            k4[3*(p-lmin):3*(p+1-lmin),3*(q-lmin):3*(q+1-lmin),3*(r-lmin):3*(r+1-lmin),3*(s-lmin):3*(s+1-lmin)] = zz2


    return k4

d_symcov12 = sf.d_symcov12
# WARNING: Provisionally adding a b_l[2]^2 to match that added to the coupling matrix...
def _d_symcov12(mat11,mat12,mat22,qty11,qty22,sih):
    lmin,lmax = sih.lmin,sih.lmax 
    mat21=transpose(mat12).copy()
    dcov=zeros((3*(lmax-lmin+1),3*(lmax-lmin+1),3*(lmax-lmin+1)))
    for l in itertools.product(xrange(lmin,lmax+1),repeat=3):
        zz1=zeros((3,3,3))
        m21_ij=matpart(l[0],l[1],mat21,sih)   
        m21_ji=matpart(l[1],l[0],mat21,sih)
        m22_ij=matpart(l[0],l[1],mat22,sih)
        m11_ji=matpart(l[1],l[0],mat11,sih)
        f22_ik=matpart(l[0],l[2],qty22,sih)
        f11_ki=matpart(l[2],l[0],qty11,sih)
        f22_jk=matpart(l[1],l[2],qty22,sih)
        f22_kj=matpart(l[2],l[1],qty22,sih)
        f11_jk=matpart(l[1],l[2],qty11,sih)
        f11_kj=matpart(l[2],l[1],qty11,sih)

        for i in itertools.product([0,1],repeat=6):
            zz1[i[0]+i[1],i[2]+i[3],i[4]+i[5]] += (
                + pairtrace(dot(f22_ik[i[1]::2,i[4]::2], f11_kj[i[5]::2,i[2]::2]), m21_ji[i[3]::2,i[0]::2])
                + pairtrace(dot(m21_ij[i[1]::2,i[2]::2], f22_jk[i[3]::2,i[5]::2]), f11_ki[i[4]::2,i[0]::2])
                + pairtrace(dot(f22_ik[i[1]::2,i[4]::2], f22_kj[i[5]::2,i[3]::2]), m11_ji[i[2]::2,i[0]::2])
                + pairtrace(dot(m22_ij[i[1]::2,i[3]::2], f11_jk[i[2]::2,i[5]::2]), f11_ki[i[4]::2,i[0]::2])
                )*sih.b_l[l[2]]**2

        dcov[3*(l[0]-lmin):3*(l[0]+1-lmin),3*(l[1]-lmin):3*(l[1]+1-lmin),3*(l[2]-lmin):3*(l[2]+1-lmin)]=zz1
    #dcov.tofile("dcov_slow.dat")
    return dcov


def s_romb_part(cls,sih):
    lmin,lmax = sih.lmin,sih.lmax
    grads=zeros(5)
    dcls=(cls-sih.fidcls)
    partdcls=dcls[3*(sih.lminlike-lmin):3*(sih.lmaxlike-lmin+1)].copy()
    grads[0]=dot(transpose(partdcls),sih.fidgrads_part)
    grads[1]=dot(transpose(partdcls),grad_s_part(sih.ci,sih.fidcls+.25*dcls,sih))
    grads[2]=dot(transpose(partdcls),grad_s_part(sih.ci,sih.fidcls+.50*dcls,sih))
    grads[3]=dot(transpose(partdcls),grad_s_part(sih.ci,sih.fidcls+.75*dcls,sih))
    grads[4]=dot(transpose(partdcls),grad_s_part(sih.ci,cls,sih))
    return scipy.integrate.romb(grads,.25,show=False)

