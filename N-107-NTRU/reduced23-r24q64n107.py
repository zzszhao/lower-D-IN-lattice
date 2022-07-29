import random as random321
import sys


from sage.all import *
import time
import datetime

global n,N,Nmax,q,p
n=N=107#19,37,59,211,107,311
Nmax=n+1
q=64
p=3

global x11,x12,x21,x22
x11=x12=x21=x22=0

global jishu
global flag_bukeni


def fuzhi_nomod(a,b,modwho):
   for i in range(1,Nmax):
      a[i]=(b[i]+modwho)%modwho
   a[0]=b[0]

def poly_copy(a,b,Nmax):
    for i in range(Nmax):
        a[i]=b[i]
def Poly_init(r,n):
   for i in range(0,n):
      r[i]=0
def fuzhi(a,b,modwho):
   for i in range(0,Nmax):
      a[i]=(b[i]+modwho)%modwho

def deg(a):
    for i in range(N,-1,-1):
        if(a[i]!=0):
            return i
    return 0

def StarMultiply(a,b,c,n,M):
    for k in range(n-1,-1,-1):
        c[k]=0
        j=k+1
        for i in range(n-1,-1,-1):
            if(j==n):
                j=0
            if( (a[i]!=0) and (b[j]!=0)  ):
                c[k]=(c[k]+(a[i]*b[j]))%M
            if (c[k]<0):
                c[k]=c[k]+M
            j=j+1
def  Poly_Xor(a,b,Nmax):
    for i in range(Nmax):
        a[i]=(a[i]+b[i])%2
def Poly_Exchange(a,b):
    temp=[0 for x in range(0,Nmax)]
    poly_copy(temp,a,Nmax)
    poly_copy(a,b,Nmax)
    poly_copy(b,temp,Nmax)

def Poly_Mult_Num(c,b,a):
    for i in range(0,Nmax):
        b[i]=c*a[i]
def Poly_Minus(a,b,p):
    for i in range(0,Nmax):
        a[i]=(a[i]-b[i])%p


def Inverse_Poly_Fq(a,Fq,q):
    k=v=0
    Poly_init(Fq,Nmax)
    g=[0 for x in range(0,Nmax)]
    g[0]=-1
    g[N]=1
    b=[0 for x in range(0,Nmax)]
    b[0]=1
    c=[0 for x in range(0,Nmax)]
    f=[0 for x in range(0,Nmax)]
    poly_copy(f,a,Nmax)
    temp2=[0 for x in range(0,Nmax)]
    temp=[0 for x in range(0,Nmax)]
    while(1):
        while((f[0]==0) and (deg(f)!=0) ):
            for i in range(1,Nmax):
                f[i-1]=f[i]
                c[N+1-i]=c[N-i]
            f[N]=0
            c[0]=0
            k=k+1
        if(deg(f)==0):
            j=0
            k=k%N
            for i in range(N-1,-1,-1):
                j=i-k
                if(j<0):
                    j=j+N
                Fq[j]=b[i]
            v=2
            while(v<q):
                v=v*2
                StarMultiply(a,Fq,temp,N,v)
                StarMultiply(temp,Fq,temp2,N,v)
                Poly_Mult_Num(2,Fq,Fq)
                Poly_Minus(Fq,temp2,v)
            for i in range(N-1,-1,-1):
                if(Fq[i]<0):
                    Fq[i]=Fq[i]+q
            return 1
        else:
            if(deg(f)<deg(g) ):
                Poly_Exchange(f,g)
                Poly_Exchange(b,c)
            Poly_Xor(f,g,Nmax)
            Poly_Xor(b,c,Nmax)
            





global poly_1,poly_0,poly_N
global poly_x11,poly_x12,poly_x21,poly_x22
global poly_x11_x22,poly_x12_x21
global poly_jian1, poly_ab,poly_q,poly_r

global poly_yanzheng,poly_Fp,poly_Fq


poly_N=[-1]
for x in range(0,n-1):
    poly_N.append(0)
poly_N.append(1)


poly_1=[1]
poly_0=[0 for x in range(0,Nmax)]
for x in range(0,n):
    poly_1.append(0)

poly_x11=[0 for x in range(0,Nmax)]
poly_x12=[0 for x in range(0,Nmax)]
poly_x21=[0 for x in range(0,Nmax)]
poly_x22=[0 for x in range(0,Nmax)]
poly_x11_x22=[0 for x in range(0,Nmax)]
poly_x12_x21=[0 for x in range(0,Nmax)]

poly_jian1=[0 for x in range(0,Nmax)]
poly_ab=[0 for x in range(0,Nmax)]
poly_q=[0 for x in range(0,Nmax)]
poly_r=[0 for x in range(0,Nmax)]


poly_yanzheng=[0 for x in range(0,Nmax)]
poly_Fq=[0 for x in range(0,Nmax)]
poly_Fp=[0 for x in range(0,Nmax)]




def RandPoly( r, n , NumOnes, NumNegOnes):
   for i in range(0,n):
      r[i]=0
   while(NumOnes!=0 or NumNegOnes!=0):
      position=random321.randint(0,n-1)
      if r[position] == 0:
         if NumOnes > 0:
            r[position] = 1
            NumOnes = NumOnes-1
         elif NumNegOnes > 0:
            r[position] = -1
            NumNegOnes = NumNegOnes-1



def maxlen(a,modwho):
   t_maxlen=0
   for i in range(1,Nmax):
      if (a[Nmax-i])%modwho != 0:
         t_maxlen=Nmax-i
         break
      else:
         continue
   return t_maxlen

def poly_mempcy(a,b):
   t=0
   for i in range(0,Nmax):
      if a[i]==b[i]:
         t=t+1
   if t==Nmax:
      return 1
   else :
      return 0

def poly_mempcy_xunhuan(a,b):
   g=len(a)

   for j in range(0,g):
      t=0
      tt=0
      for i in range(0,g):
         if a[(i+j)%g]==b[i]:
            t=t+1
         if a[(i+j)%g]==-b[i]:
            tt=tt+1
      if t==g or tt==g:
         return 1
   return 0

def poly_mempcy_reduced_f(a,b,arrf):
   arrf.reverse()
   g=len(arrf)
   for j in range(g):

       b.pop(arrf[j])

   if poly_mempcy_xunhuan(a,b)==1:
        arrf.reverse()
        for k in range(g):
            b.insert(arrf[k],0)
        return 1
   arrf.reverse()
   for k in range(g):
       b.insert(arrf[k],0)


def poly_mempcy_reduced_a(a,b,arrf,n):
   c=[]
   g=len(arrf)
   for i in range(n-len(arrf)):
      c.append(a[i])
   for k in range(g):
       c.insert(arrf[k],0)
   if poly_mempcy_xunhuan(c,b)==1:
      return 1
   else:
      return 0




def poly_mul_modxn(a,b,ab,modwho):
   for i in range(0,Nmax):
      ab[i]=0
   for i in range(0,Nmax):
      for j in range(0,Nmax):
         if ((i+j)>=n):
            ab[(i+j)%n]+=a[i]*b[j]
         else:
            ab[i+j]+=a[i]*b[j]
   for i in range(0,Nmax):
      ab[i]=(ab[i]+modwho)%modwho

def poly_sub(a,b,a_sub_b,modwho):
   for i in range(0,Nmax):
      a_sub_b[i]=(a[i]-b[i]+modwho)%modwho

def poly_add(a,b,a_add_b,modwho):
   for i in range(0,Nmax):
      a_add_b[i]=(a[i]+b[i])%modwho

def shu_mul_poly(a, b, modwho):
   for i in range(0,Nmax):
      b[i]=(b[i]*a)%modwho



def qinshi(a,m):
   global x11,x12,x21,x22
   x12=a
   x11=1
   x22=m
   x21=0
   qq=0
   r=0
   while x12 != x22:
      if x22>x12:
         qq=(x22-1)//x12
         r=x22-qq*x12
         x21=qq*x11+x21
         x22=r
      elif x12>x22:
         qq=(x12-1)//x22
         r=x12-qq*x22
         x11=qq*x21+x11
         x12=r
   return x11

def jian1_poly(a,b,modwho):
   for i in range(0,Nmax):
      a[i]=b[i]
   a[0]=(a[0]-1+modwho)%modwho


def qin_polymod(a,b,c,d,modwho):
   global poly_0
   global x11
   fuzhi(c,poly_0,modwho)
   
   if(maxlen(a,modwho)<maxlen(b,modwho)):
      fuzhi(d,poly_0,modwho)
      fuzhi(c,a,modwho)
      c[0]=c[0]+1
      return
   elif(maxlen(b,modwho)==0 and b[0]==0):
      print("error b=0")
      print(modwho)
      return
   
   elif(maxlen(b, modwho)==0 and b[0]==modwho ):
      print("irreversible")
      fuzhi_nomod(c,b,modwho)
      return

   elif(maxlen(a, modwho) >= maxlen(b, modwho)):
      poly_gap=0
      poly_gap=maxlen(a, modwho)-maxlen(b,modwho)
      qinshi(b[maxlen(b, modwho)],modwho)
      x11=x11%modwho
      d[poly_gap]=(a[maxlen(a, modwho)]*x11)%modwho
      for i in range(poly_gap,maxlen(a, modwho)+1):
         a[i]=(a[i]-a[maxlen(a,modwho)]*x11*b[i-poly_gap]+a[maxlen(a, modwho)]*x11*modwho)%modwho
      if (maxlen(a,  modwho) < maxlen(b, modwho)):
         fuzhi(c, a, modwho)
         c[0]=(c[0]+1)
         return
      elif (maxlen(a, modwho)==0 and a[0]==0):
         fuzhi(c,a,modwho)
         c[0]=(c[0]+1)
         return
      elif (maxlen(a, modwho) >= maxlen(b,  modwho)):
         qin_polymod(a,b,c,d,modwho)


def qinshi_poly(a,m,modwho):
   global poly_1, poly_0,poly_x11,poly_x21,poly_x12,poly_x22
   global flag_bukeni
   global x11
   flag_bukeni=0

   fuzhi(poly_x11, poly_1, modwho)
   fuzhi(poly_x21, poly_0, modwho)
   fuzhi(poly_x12,a, modwho)
   fuzhi(poly_x22,m, modwho)

   poly_q=[0 for x in range(0,Nmax)]
   poly_r=[0 for x in range(0,Nmax)]

   while( poly_x12[0]!=1 or maxlen(poly_x12,modwho)!=0 ):
      if( maxlen(poly_x22,modwho)==0 and poly_x22[0]==modwho):
         print(" poly_x22  irreversible")
         print(a)
         flag_bukeni=flag_bukeni+1
         break
      elif( maxlen(poly_x12,modwho)==0 and poly_x12[0]==modwho):
         print("poly_x12  irreversible")
         print(a)
         flag_bukeni=flag_bukeni+1
         break
      elif( maxlen(poly_x22,modwho)>maxlen(poly_x12,modwho) or
           ( maxlen(poly_x22,modwho)==maxlen(poly_x12,modwho) and
            poly_x22[maxlen(poly_x22,modwho)] >= poly_x12[maxlen(poly_x12,modwho )])   ):
            jian1_poly(poly_jian1, poly_x22,modwho)
            fuzhi(poly_q, poly_0, modwho)
            qin_polymod(poly_jian1, poly_x12, poly_r, poly_q,modwho)
            poly_mul_modxn(poly_q, poly_x11, poly_ab,modwho)
            poly_add(poly_ab, poly_x21, poly_x21,modwho)
            fuzhi_nomod(poly_x22, poly_r, modwho )
      elif( maxlen(poly_x12,modwho)>maxlen(poly_x22,modwho) or
          ( maxlen(poly_x12,modwho )==maxlen(poly_x22,modwho ) and
           poly_x12[maxlen(poly_x12,modwho )] > poly_x22[maxlen(poly_x22,modwho )])   ):
           jian1_poly(poly_jian1, poly_x12,modwho)
           fuzhi(poly_q, poly_0, modwho)
           qin_polymod(poly_jian1, poly_x22,  poly_r, poly_q,modwho)
           poly_mul_modxn(poly_q, poly_x21, poly_ab,modwho)
           poly_add(poly_ab, poly_x11, poly_x11,modwho)
           fuzhi_nomod(poly_x12, poly_r, modwho )
   if(maxlen(poly_x12,modwho )==0 and poly_x12[0]==0 ):
      print("error")
   elif(maxlen(poly_x12,modwho )==0 and poly_x12[0]!=1):
      qinshi(poly_x12[0], modwho)
      shu_mul_poly(x11, poly_x11,modwho)
   

def panduan_f_inv_pq(r,modp,modq):
   t=0
   global flag_bukeni
   qinshi_poly(r, poly_N, modp)
   if(flag_bukeni!=0):
      return 0
   poly_mul_modxn(poly_x11, r, poly_yanzheng,modp )
   if(poly_mempcy(poly_yanzheng, poly_1)==1 ):
      fuzhi(poly_Fp, poly_x11, modp)
      t=t+1
   if(modq%2==1):
       qinshi_poly(r, poly_N, modq)
       if(flag_bukeni!=0):
           return 0
       poly_mul_modxn(poly_x11, r, poly_yanzheng,modq )
       if(poly_mempcy(poly_yanzheng, poly_1)==1 ):
           fuzhi(poly_Fq, poly_x11,modq)
           t=t+1
   if(modq%2==0):
       Inverse_Poly_Fq(r,Fq,modq)
       poly_mul_modxn(Fq, r, poly_yanzheng,modq )
       if(poly_mempcy(poly_yanzheng, poly_1)==1 ):
           fuzhi(poly_Fq, Fq,modq)
           t=t+1
   if(t==2):
       return 1
   elif(t!=2):
       return 0



def RandPoly_ensure_inv( r, n,  NumOnes,  NumNegOnes,  modp,  modq):
   t=0
   global poly_Fq,poly_Fp,poly_x11,poly_yanzheng,poly_N
   
   poly_x11_1=[0 for x in range(0,Nmax)]
   poly_x11_2=[0 for x in range(0,Nmax)]
   while(t!=2):
      t=0

      RandPoly( r, n,  NumOnes,  NumNegOnes)

      qinshi_poly(r, poly_N, modp)
      fuzhi(poly_x11_1,poly_x11,modp)
      if(modq%2==1):
          qinshi_poly(r, poly_N, modq)
          fuzhi(poly_x11_2,poly_x11,modq)
      elif(modq%2==0):
          Inverse_Poly_Fq(r,Fq,modq)
          fuzhi(poly_x11_2,Fq,modq)

      if(flag_bukeni!=0):
         while(flag_bukeni!=0):
            RandPoly( r, n,  NumOnes,  NumNegOnes)
            qinshi_poly(r, poly_N, modp)
            fuzhi(poly_x11_1,poly_x11,modp)
            if(modq%2==1):
                qinshi_poly(r, poly_N, modq)
                fuzhi(poly_x11_2,poly_x11,modq)
            elif(modq%2==0):
                Inverse_Poly_Fq(r,Fq,modq)
                fuzhi(poly_x11_2,Fq,modq)        
      poly_mul_modxn(poly_x11_1, r, poly_yanzheng,modp )
      if(poly_mempcy(poly_yanzheng, poly_1)==1 ):
         fuzhi(poly_Fp, poly_x11_1, modp)
         t=t+1
      poly_mul_modxn(poly_x11_2, r, poly_yanzheng,modq )

      if(poly_mempcy(poly_yanzheng, poly_1)==1 ):
         fuzhi(poly_Fq, poly_x11_2, modq)
         t=t+1




def CreateKey(aa, bb, cc, dd, ee, ff, gg, hh):
   t=0
   global n,N,q
   global poly_Fp,poly_Fq
   if(panduan_f_inv_pq(dd,p,q)!=1 ):
      if(n==107):
         RandPoly_ensure_inv(dd,N,15,14,cc,bb)#n=107
      elif(n==97):
         RandPoly_ensure_inv(dd,N,14,13,cc,bb)#n=97
      elif(n==83):
         RandPoly_ensure_inv(dd,N,12,11,cc,bb)#n=83
   fuzhi(gg,poly_Fp, cc)
   fuzhi(hh,poly_Fq, bb)
   poly_mul_modxn(hh,ee,ff, q)




ff=[0 for x in range(0,Nmax)]
gg=[0 for x in range(0,Nmax)]
hh=[0 for x in range(0,Nmax)]
Fp=[0 for x in range(0,Nmax)]
Fq=[0 for x in range(0,Nmax)]
flag=[0 for x in range(0,Nmax)]

H= [[0 for i in range(Nmax)] for j in range(Nmax)]
zhuanzhi_H=[[0 for i in range(Nmax)] for j in range(Nmax)]


seconds=0

if __name__ == '__main__':
   success=0
   num=0
   aaa=[0,1,0,0]
   bbb=[0,0,1]
   ccc=[0]
   print(poly_mempcy_reduced_a(bbb,aaa,ccc,4))
   print("*****")
   print(aaa)
   print(bbb)
   print("*****")



   begin = datetime.datetime.now()
   while(num!=100):
      num=num+1
      arr=[0, 3, 7 , 16, 25, 32, 37, 45, 48, 51, 65, 67, 71, 85, 91]
      
      arrf=[0,3,11,16, 25, 32, 37, 45, 48, 51]#
      

      f_norm=0
      if(n==107):
         #n=107
         RandPoly(ff,N,15,14)
         RandPoly(gg,N,12,12)
         f_norm=sqrt(29)
      elif(n==97):
         RandPoly(ff,N,14,13)
         RandPoly(gg,N,12,12)
         f_norm=sqrt(27)
         for i in range(0,1):
            arr.pop()
      elif(n==83):
         RandPoly(ff,N,12,11)
         RandPoly(gg,N,10,10)
         f_norm=sqrt(23)
         for i in range(0,3):
            arr.pop()
      
      
      CreateKey(N,q,p,ff,gg,hh,Fp,Fq)
      
      
      f=[0 for x in range(0,N)]
      g=[0 for x in range(0,N)]
      h=[0 for x in range(0,N)]

      for i in range(0,N):
         f[i]=ff[i]
         g[i]=gg[i]
         h[i]=hh[i]

      position=[]

      j=0
      for i in range(len(g)):
         if(g[i]==0):
            position.insert(j,i)
            j=j+1
         if(j==24):
            break



      for i in range(len(arrf)):
         arrf[i]=position[i]

      position_f=[]
      jj=0
      for i in range(len(f)):
         if(f[i]==0):
            position_f.insert(jj,i)
            jj=jj+1

      position_f_bufen=[]
      for i in range(23):
         position_f_bufen.insert(i,position_f[i])







      for i in range(0,N):
         for j in range(0,N):
            k=(i+j)%n
            H[i][k]=h[j]
            zhuanzhi_H[k][i]=h[j]

      H_reduced=[[0 for i in range(N-len(position_f_bufen))] for j in range(len(position))] 
   
      for m in range(0,len(position)):
         jjj=0
         for i in range(0,N):
            if i in position_f_bufen:
               continue
            H_reduced[m][jjj]=zhuanzhi_H[position[m]][i] 
            jjj=jjj+1

      success_neibu=0
      jishu_neibu=0
      zhuanzhi_H_zhongjian=[[0 for i in range(Nmax)] for j in range(Nmax)] 
      print("*****simultion parallel*****")
      

    
      for m in range(0,len(position)):
         jjj=0
         for i in range(0,N):
            if i in position_f_bufen:
               continue
            H_reduced[m][jjj]=zhuanzhi_H[position[m]][i]  
            jjj=jjj+1
      if(q%2!=0):
         H_reduced_sage=matrix(Zmod(q),H_reduced)
      elif(q%2==0):
         H_reduced_sage=matrix(ZZ,H_reduced)
      Anull = H_reduced_sage.right_kernel_matrix()
      

      LambdaA=block_matrix(2,1,[Anull.lift(),q*identity_matrix(n-len(position_f_bufen))])

      begin_bufen = time.time()
      if(q%2==0):
         #reduced=LambdaA.BKZ(algorithm="NTL", fp='xd',block_size=10, prune=0)
         reduced=LambdaA.BKZ(algorithm="NTL", fp='qd1',block_size=40, prune=15)
      if(q%2==1):
         reduced=LambdaA.BKZ(algorithm="NTL", fp='qd1',block_size=40, prune=15)
      end_bufen = time.time()
      print('one seconds='+str(end_bufen-begin_bufen))
      seconds=seconds+end_bufen-begin_bufen
      min_norm=min(filter(lambda x:x, [r.norm() for r in reduced]))
      
      print(min_norm)
      if (min_norm>f_norm):
         continue
      for r in reduced:
        if r.norm() == f_norm:
           if poly_mempcy_reduced_a(r,f,position_f_bufen,n)==1:
              print(r[0])

              success=success+1
              success_neibu=success_neibu+1
              break
      print("*****simultion end*****")
      

   print('success='+str(success))
   print('all bkz seconds='+str(seconds ))
   print('num='+str(num))
   end = datetime.datetime.now()
   print('all seconds='+str((end-begin).seconds))
   
   filename = 'reduced23-r24n107q64.txt'
   with open(filename, 'w') as file_object:
      file_object.write('all seconds='+str((end-begin).seconds)+"\n")
      file_object.write('all bkz seconds='+str(seconds )+"\n")
      file_object.write("success="+str(success)+"\n")
      file_object.write("num="+str(num)+"\n")
      file_object.write("seconds="+str((end-begin).seconds))
      
   file_object.close()   

