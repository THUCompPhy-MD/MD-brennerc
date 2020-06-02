
       common /integ/ifail,ip,iprint,neigen,iprintstep,
     .            iik,Nntb(natx),Nident(natx,nnx),
     .            NE_yli,Nat,ierr 
       common /lldos/xdos(512),idos(512)      
       common /reall/tbpar(6),p,X02AAF,summa,xkt(10,3),ecl,
     .            HR(nmax,nmax),E(nmax),VR(nmax,nmax),
     .            xnaapurit(natx,nnx,3),scale,rcut,basis(3,3),
     .            E_pot,ed_F(natx,3),ann,
     .            xnopeudet(natx,3),E_kin,dist,vec1(3),vec2(3),
     .            temp,seed,xscale,Emintb,Emax,E_scale, 
     .            fv1(nmax),fv2(nmax), 
     .            V_mkp(3),apu, 
     .            DH(4,4,3),FHF(natx,3),capu2,capu3,capu4, 
     .            wsum(natx),dwsum(natx) 
     .       ,vatb(nmax,nmax) 
c       common / par/ Am,rc,ro,xnc,xm,xmc,do,dc,aa1,aa2 
c     .,f0,f1,f2,f3,f4,s0,s1,s2,s3,t0,t1,t2,t3  
c carbon
       data tbpar/-2.99,3.71,-5.0,4.7,5.5,-1.55/ 

       logical    cont
