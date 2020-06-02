      subroutine PAIRS

      implicit real*8 ( a - h, o - z )

C  This subroutine calculates the forces on the atoms according to
C  the Murty-Atwater-Potential.
C  All variables are declared in common-blocks as include-files.
C  The potential energy is not calculated.
C
C  The atomic pairs are have to be calculated with a linked-cell
C  algorithm and stored in the following way:
C  (i) Two body pairs: nop is the number of 2body pairs,
C       the index of the first atom is in nlista, the second atom
C       in nlistb. If the pair (i,j) is stored, then the pair (j,i) is 
C        not stored. Example: nlista(i)=3,nlistb(i)=7 --> the pair
C        with index i consists of atom number 3 and atom number 7.
C  (ii) Three body pairs: the number of 3body pairs is nod, the
C       index of the first atom is in nlist1,the index of the second in 
C       nlist2 and the index of the third in nlist3. Here, all pairs
C       are stored, i.e., if (i,j,k) is stored then (i,k,j) is also
C       stored.    

C --- Global variables ---
C     nmax                      =   maximum number of atoms
C     nmol                      =   actual number of atoms 
C     nop                       =   number of 2body pairs
C     nod                       =   number of 3body pairs 
C     pos(3,nmax)               =   x,y,z coordinates of the atoms
C     fx(nmax),fy(nmax),fz(nmax)=   x,y,z components of the forces
C     ltype(nmax)               =   integer corresponding to the atom type 
C     xl,yl,zl                  =   box lengths in cm
C     xl2,yl2,zl2               =   2.0/xl, 2.0/yl, 2.0/zl
C     k4(mi,mj,mk)              =   mi,mj,mk are the atom type, k4 labels the kind
C                                   of potential
C     ktype(mi.mj)              =   the same for the parameters which depend only
C                               =   on two types of atoms 

C --- Parameters of the potential, notation
C     according to eqns. (2) and (3)
C     ta  = A
C     tb  = B_0 
C     tl  = lambda_1 
C     tm  = lambda_2
C     tn  = eta
C     tn2 = -delta
C     req = (R_ij)**(e)
C     tc  = c
C     td  = d 
C     hn  = h
C     alf = alpha
C     bet = beta

C --- local variables and arrays, * depends on the number of atoms
C                                 and how many pairs one expects
C     xi(*)      =  Zeta_ij
C     dxi(*,6)   =  derivatives of Zeta_ij with respect
C                   to the coordinates x_i,y_i,z_i,x_j,y_j,z_j
C     bijx(*)    =  the b_ij- term from the Tersoff-notation
C     rrij(*)    =  store some distances 
C     nhelp1(*)  =  first partner and
C     nhelp2(*)  =  second partner for the contributions Phi_ij
C     f3(*)      =  store exponential terms
C     g3(*)      =  store angular terms
C     xn(*,2)    =  number of neighbours weighted with cutoff (=N in eqn.(4)) 

      dimension xi(lp3*nmax),dxi(lp3*nmax,6),bijx(lp3*nmax)
      dimension rrij(lp3*nmax),nhelp1(60*nmax),nhelp2(60*nmax) 
      dimension xn(nmax,2), f3(lp3*nmax),g3(lp3*nmax)


C     initialize arrays

      do 100 i = 1, nmol
        fx(i) = 0 d0
        fy(i) = 0 d0
        fz(i) = 0 d0
        xn(i,1) = 0 d0
        xn(i,2) = 0 d0                   ! is not used here
  100 continue    
      epot=0.0

      do 200,i=1,lp3*nmol
           xi(i)=0.0
           xnij(i)=1.0
           do 220,k=1,6
              dxi(i,k)=0.0
220        continue
200   continue

C     First, count neighbours according to eqn. (4)

      do 300,i=1,nop                     ! use the 2body list
         i1=nlista(i)
         i2=nlistb(i)
         xij=pos(1,i1)-pos(1,i2)
         yij=pos(2,i1)-pos(2,i2)
         zij=pos(3,i1)-pos(3,i2)
C --- periodic boundaries-------
         mx = xij * xl2
         xij = xij - mx * xl
         my = yij * yl2
         yij = yij - my * yl
         mz = zij * zl2
         zij = zij - mz * zl * isurf     ! isurf=0 if no pbc in z-direction
C -------------------------------        ! isurf=1 if pbc in z-direction 
         rij2 = xij * xij + yij * yij + zij * zij
         rij=dsqrt(rij2)
         mi=ltype(i1)                    ! What kind of atom ?
         mj=ltype(i2)
         kp=ktype(mi,mj)                 ! What kind of interaction ?
         xn(i1,1)=xn(i1,1)+fc(rij,kp)    ! Count neighbours, weight with cutoff
300     continue

C     Compute Zeta_ij, every different pair (i,j) gets a label
C     store number of different pairs (ipaar)
C     store angular term c+d(H-cos(theta))**2  (eqn.(3)) in g3
C     store exponential term (eqn.(3)) in f3
      ipaar=0
      iold=0
      jold=0
      do 710 inod = 1, nod               ! use the 3body list
        i = nlist1(inod)                 ! first partner
        j = nlist2(inod)                 ! second partner 
        k = nlist3(inod)                 ! third partner
        mi = ltype(i)
        mj = ltype(j)
        k2ij=ktype(mi,mj)
C ---   pbc -------
        xii = pos(1,i)
        yii = pos(2,i)
        zii = pos(3,i)
        xij = xii - pos(1,j)
        yij = yii - pos(2,j) 
        zij = zii - pos(3,j) 
        mx = xij * xl2
        xij = xij - mx * xl
        my = yij * yl2
        yij = yij - my * yl
        mz = zij * zl2
        zij = zij - mz * zl * isurf
C -------------------
        rij2 = xij * xij + yij * yij + zij * zij
        if ( rij2 .LE. rlc2(k2ij) ) then
           if ((j.ne.jold).or.(i.ne.iold)) then
             ipaar = ipaar+1              ! counts and labels the 3body pairs
             nhelp1(ipaar)=i              ! with different 2body parts
             nhelp2(ipaar)=j              ! (1,2,3) has the same label as
             jold=j                       ! (1,2,4), but (1,3,7) gets a new
             iold=i                       ! label
          endif
          if (k.eq.0) then                ! no third partner
             xi(ipaar)=0.0
             rrij(ipaar)=dsqrt(rij2)      ! store the distance for later use
             f3(inod)=0.0
             g3(inod)=0.0
          else
          mk = ltype(k)
          k2ik=ktype(mi,mk)
          k2jk=ktype(mj,mk)
          kp = k4(mj,mi,mk)
          xik = xii - pos(1,k) 
          yik = yii - pos(2,k) 
          zik = zii - pos(3,k) 
          mx = xik * xl2
          xik = xik - mx * xl
          my = yik * yl2
          yik = yik - my * yl
          mz = zik * zl2
          zik = zik - mz * zl * isurf
          rik2 = xik * xik + yik * yik + zik * zik
          ajik = xij * xik + yij * yik + zij * zik                    
       
          xjk = pos(1,j) - pos(1,k)
          yjk = pos(2,j) - pos(2,k)
          zjk = pos(3,j) - pos(3,k)
          mx = xjk * xl2
          xjk = xjk - mx * xl
          my = yjk * yl2
          yjk = yjk - my * yl
          mz = zjk * zl2
          zjk = zjk - mz * zl * isurf
          rjk2 = xjk * xjk + yjk * yjk + zjk * zjk
          rjk  = dsqrt(rjk2)

             rij = dsqrt(rij2)
             rrij(ipaar)=rij                     ! store distance
             if (rij.lt.0.5d-8) write(*,*) 'BUNG! BUNG!',i,j,mi,mj,rij  
                                                 ! something strange happened
             rik = dsqrt(rik2)
             cosin=ajik/(rij*rik)
             zz=alf(kp)*(rij-rik+req(k2ik)-req(k2ij))**bet(kp)
             f3(inod)=Dexp(zz)
             g3(inod)=tc(kp)+td(kp)*(hn(xn(i,1),kp)-cosin)**2.0
             xi(ipaar)=xi(ipaar)+fc(rik,k2ik)*g3(inod)*f3(inod)
      endif 
      endif
  710 continue                           

C     Compute the bond order term of the potential :
C     F_2 (1+Zeta_ij **eta)**(-delta)
C     It's called 'bijx' according to the usual Tersoff notation (b_ij)
     
      do 720,i=1,ipaar                   
         m1=ltype(nhelp1(i))             
         m2=ltype(nhelp2(i))
         kp=ktype(m1,m2) 
         if (xi(i).ne.0.0) then
         bijx(i)=(1.0d0+xi(i)**tn(kp))**tn2(kp)
         bijx(i)=f2(xn(nhelp1(i),1),kp)*bijx(i) 
         else
          bijx(i)=f2(xn(nhelp1(i),1),kp)
         endif
720   continue

C     Compute the derivatives d Phi_ij / d x_k and so on   
C     only the bond order term contributes

      ip=0
      iold=0
      jold=0 
      do 800,inod=1,nod
         i=nlist1(inod)
         j=nlist2(inod)
         k=nlist3(inod)

         mi = ltype(i)
         mj = ltype(j)
         k2ij=ktype(mi,mj)

         xii = pos(1,i)
         yii = pos(2,i)
         zii = pos(3,i)
         xij = pos(1,j) - xii
         yij = pos(2,j) - yii
         zij = pos(3,j) - zii
         mx = xij * xl2
         xij = xij - mx * xl
         my = yij * yl2
         yij = yij - my * yl
         mz = zij * zl2
         zij = zij - mz * zl * isurf
         rij2 = xij * xij + yij * yij + zij * zij
         if ( rij2 .LE. rlc2(k2ij) ) then
            if ((j.ne.jold).or.(i.ne.iold)) then
               ip=ip+1
               jold=j
               iold=i
            endif
         if (k.ne.0) then 
            mk=ltype(k)
            k2ik=ktype(mi,mk)
            kp = k4(mj,mi,mk)
            xik = pos(1,k) - xii
            yik = pos(2,k) - yii
            zik = pos(3,k) - zii
            mx = xik * xl2
            xik = xik - mx * xl
            my = yik * yl2
            yik = yik - my * yl
            mz = zik * zl2
            zik = zik - mz * zl * isurf
            rik2 = xik * xik + yik * yik + zik * zik
            ajik = xij * xik + yij * yik + zij * zik
            rij = rrij(ip)
            rik = dsqrt(rik2)
            cosin=ajik/(rij*rik)
            xi0=0.0
            if (xi(ip).gt.0.0)
     &      xi0= xi(ip)**(tn(k2ij)-1.0)                        ! if there is no third partner xi(ip) is zero


            hilf=-2.0*td(kp)*(hn(xn(i,1),kp)-cosin)*fc(rik,k2ik)     !=dg3/dcos(theta)   
                                                                     ! hilf means help
            helpx=hilf*(xij/(rij*rik)-cosin*xik/rik2)*f3(inod)       !=hilf*dtheta/dx_k *f3
            helpy=hilf*(yij/(rij*rik)-cosin*yik/rik2)*f3(inod)       !=hilf*dtheta/dy_k *f3
            helpz=hilf*(zij/(rij*rik)-cosin*zik/rik2)*f3(inod)       !=hilf*dtheta/dz_k *f3

            hg=0.5*tn2(k2ij)*(1.0+xi(ip)**tn(k2ij))**(tn2(k2ij)-1.0)*
     &          xi0*tn(k2ij)*
     &          tb(k2ij)*Dexp(-tm(k2ij)*rrij(ip))*f2(xn(i,1),k2ij)             ! exp(-rij) *...* db_ij/dxi_ij   
            ghx=-f3(inod)*alf(kp)*bet(kp)*
     &                (rij-rik+req(k2ik)-req(k2ij))**(bet(kp)-1.0)         ! df3/dx_k 
     &                /rik 
            ghy=ghx*yik
            ghz=ghx*zik
            ghx=ghx*xik
            fx(k)=fx(k)+ fc(rij,k2ij)* hg*(helpx+g3(inod)*ghx)
            fy(k)=fy(k)+ fc(rij,k2ij)* hg*(helpy+g3(inod)*ghy) 
            fz(k)=fz(k)+ fc(rij,k2ij)* hg*(helpz+g3(inod)*ghz)
         endif
         endif
 800  continue                                              ! finished Forces d Phi_ij/dx_k

C     Compute derivatives d Zeta_ij /d x_i, d Zeta_ij /d x_j 
      ip=0
      iold=0
      jold=0
      do 900,inod=1,nod
         i=nlist1(inod)
         j=nlist2(inod)
         k=nlist3(inod)
         mi = ltype(i)
         mj = ltype(j)
         k2ij=ktype(mi,mj)

         xii = pos(1,i)
         yii = pos(2,i)
         zii = pos(3,i)
         xij = xii - pos(1,j) 
         yij = yii - pos(2,j) 
         zij = zii - pos(3,j) 
         mx = xij * xl2
         xij = xij - mx * xl
         my = yij * yl2
         yij = yij - my * yl
         mz = zij * zl2
         zij = zij - mz * zl * isurf
         rij2 = xij * xij + yij * yij + zij * zij
         if ( rij2 .LE. rlc2(k2ij) ) then
            if ((j.ne.jold).or.(i.ne.iold)) then
               ip=ip+1
               jold=j
               iold=i
            endif
         if (k.ne.0) then
             mk = ltype(k)
             k2ik=ktype(mi,mk)
             kp = k4 (mj,mi,mk) 
            xik = xii - pos(1,k) 
            yik = yii - pos(2,k) 
            zik = zii - pos(3,k) 
            mx = xik * xl2
            xik = xik - mx * xl
            my = yik * yl2
            yik = yik - my * yl
            mz = zik * zl2
            zik = zik - mz * zl * isurf
            rik2 = xik * xik + yik * yik + zik * zik
            ajik = xij * xik + yij * yik + zij * zik
            rij = rrij(ip)
            rik = dsqrt(rik2)
            cosin=ajik/(rij*rik)

            hilf=-2.0*td(kp)*(hn(xn(ip,1),kp)-cosin)*f3(inod)
     &            *fc(rik,k2ik)                                    !=dg3/dcos(theta) *f3
            zz= alf(kp)*bet(kp)*
     &            (rij-rik+req(k2ik)-req(k2ij))**(bet(kp)-1.0)

            pjx=-f3(inod)*g3(inod)*zz/rij
            pjy=pjx*yij
            pjz=pjx*zij
            pjx=pjx*xij

            pix=f3(inod)*g3(inod)*zz
            piy=pix*(yij/rij-yik/rik)
            piz=pix*(zij/rij-zik/rik)
            pix=pix*(xij/rij-xik/rik)

            q1= hilf*((xij+xik)/(rij*rik)
     &             +cosin* (-xij/rij2-xik/rik2))             !=dtheta/dx_i
            q2= hilf*((yij+yik)/(rij*rik)
     &              +cosin*(-yij/rij2-yik/rik2))             !=dtheta/dy_i
            q3= hilf*((zij+zik)/(rij*rik)
     &              +cosin*(-zij/rij2-zik/rik2))             !=dtheta/dz_i
            q4= hilf*(-xik/(rij*rik)
     &              +cosin*xij/rij2)                        !=dtheta/dx_j
            q5= hilf*(-yik/(rij*rik)
     &              +cosin*yij/rij2)                        !=dtheta/dy_j
            q6= hilf*(-zik/(rij*rik)
     &              +cosin*zij/rij2)                        !=dtheta/dz_j

            dxi(ip,1)=dxi(ip,1)+q1+pix
            dxi(ip,2)=dxi(ip,2)+q2+piy
            dxi(ip,3)=dxi(ip,3)+q3+piz
            dxi(ip,4)=dxi(ip,4)+q4+pjx
            dxi(ip,5)=dxi(ip,5)+q5+pjy
            dxi(ip,6)=dxi(ip,6)+q6+pjz
        endif
        endif
 900  continue                                              ! finished d Zeta/dx_i 


C     Final loop: sum up all contributions
      do 1000,inop=1,ipaar
         i=nhelp1(inop)
         j=nhelp2(inop)
         m1 = ltype(i)
         m2 = ltype(j)
         kp=ktype(m1,m2)

         xii = pos(1,i)
         yii = pos(2,i)
         zii = pos(3,i)
         xij = xii - pos(1,j)
         yij = yii - pos(2,j) 
         zij = zii - pos(3,j) 
         mx = xij * xl2
         xij = xij - mx * xl
         my = yij * yl2
         yij = yij - my * yl
         mz = zij * zl2
         zij = zij - mz * zl * isurf
         rij=rrij(inop)
         xi0=0.0
         if (xi(inop).gt.0.0) 
     &   xi0= xi(inop)**(tn(kp)-1.0)
         hilfx1  =fc(rij,kp)*0.5*tl(kp)*ta(kp)*f1(xn(i,1),kp)
     &             *Dexp(-tl(kp)*rij)/rij                      ! Contribution of the repulsive part
         hilfy1 = hilfx1*yij
         hilfz1 = hilfx1*zij
         hilfx1 = hilfx1*xij
         hilfx2 = 0.5*tb(kp)*Dexp(-tm(kp)*rij)*fc(rij,kp)*f2(xn(i,1),kp)
         dbij=tn2(kp)*(1.0+xi(inop)**tn(kp))                
     &          **(tn2(kp)-1.0)*
     &          tn(kp)*xi0                                   !dbij/dxi
 
         fix = hilfx1-hilfx2*(tm(kp)*bijx(inop)
     &               *xij/rij-dbij*dxi(inop,1))
         fiy = hilfy1-hilfx2*(tm(kp)*bijx(inop)
     &               *yij/rij-dbij*dxi(inop,2))
         fiz = hilfz1-hilfx2*(tm(kp)*bijx(inop)
     &               *zij/rij-dbij*dxi(inop,3))
         fjx = -hilfx1+hilfx2*(tm(kp)*bijx(inop)
     &               *xij/rij+dbij*dxi(inop,4))
         fjy = -hilfy1+hilfx2*(tm(kp)*bijx(inop)
     &               *yij/rij+dbij*dxi(inop,5))
         fjz = -hilfz1+hilfx2*(tm(kp)*bijx(inop)
     &               *zij/rij+dbij*dxi(inop,6))


         fx(i)=fx(i)+fix
         fy(i)=fy(i)+fiy
         fz(i)=fz(i)+fiz
         fx(j)=fx(j)+fjx
         fy(j)=fy(j)+fjy
         fz(j)=fz(j)+fjz                         
        endif
 
1000  continue                                  

      return
      end

      real*8 function fc(r,kp)
      implicit real*8(a-h,o-z)
C --  This is the cutoff function, it's straightforward
C --  parameter in common-blocks
C --  kp points on the specific potential
      if (r.gt.(rs(kp)+rr(kp)) ) fc=0.0d0
      if (r.lt.(rr(kp)-rs(kp)) ) fc=1.0d0 
      if ( (r.le.(rs(kp)+rr(kp))).and.(r.ge.(rr(kp)-rs(kp))))
     &    fc=0.5d0-9.0/16.0*Dsin(pi*(r-rr(kp))/(2.0*rs(kp)))
     &            -1.0/16.0*Dsin(3.0*pi*(r-rr(kp))/(2.0*rs(kp)))      
      return
      end   


      real*8 function f1(x,kp)
      implicit real*8(a-h,o-z)
C --  This is the Murty-Atwater F_1 function.
C --  kp=6 refers to potential IIa, otherwise the function is 1.0
C --  x corresponds to N in eqn. (4) of the paper
C --  the function splint is from 'Numerical Recipes', before using it
C --  one has to call the subroutine 'spline' once in the main program  
      if (kp.eq.6) then
         f1=splint(f1x,f1y,f1y2,5,x) 
      else
         f1=1.0d0
      endif
      return
      end
   
      real*8 function f2(x,kp)
C --  This is the function F_2, details are the same as in f1
      implicit real*8(a-h,o-z)
      include"cut.cmm"
      if (kp.eq.6) then
          f2=splint(f2x,f2y,f2y2,5,x)
      else
          f2=1.0d0
      endif 
      return
      end  
