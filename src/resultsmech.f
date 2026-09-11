!
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2025 Guido Dhondt
!     
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation(version 2);
!     
!     
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of 
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
!     GNU General Public License for more details.
!     
!     You should have received a copy of the GNU General Public License
!     along with this program; if not, write to the Free Software
!     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!     
      subroutine resultsmech(co,kon,ipkon,lakon,ne,v,
     &     stx,elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,
     &     ielmat,ielorien,norien,orab,ntmat_,t0,t1,ithermal,prestr,
     &     iprestr,eme,iperturb,fn,iout,qa,vold,nmethod,
     &     veold,dtime,time,ttime,plicon,nplicon,plkcon,nplkcon,
     &     xstateini,xstiff,xstate,npmat_,matname,mi,ielas,icmd,
     &     ncmat_,nstate_,stiini,vini,ener,eei,enerini,istep,iinc,
     &     springarea,reltime,calcul_fn,calcul_qa,calcul_cauchy,nener,
     &     ikin,nal,ne0,thicke,emeini,pslavsurf,
     &     pmastsurf,mortar,clearini,nea,neb,ielprop,prop,kscale,
     &     list,ilist,smscale,mscalmethod,enerscal,t0g,t1g,
     &     islavquadel,aut,irowt,jqt,mortartrafoflag,
     &     intscheme,physcon,dam,hasdamage,de12onepass,de12tangent,
     &     ndmat_,ndmcon,dmcon,dambase,damjac,damvisc,damviscini,
     &     damvisceta)
!     
!     calculates stresses and the material tangent at the integration
!     points and the internal forces at the nodes
!     
      implicit none
!     
      character*8 lakon(*),lakonl
      character*80 amat,matname(*)
!     
      integer kon(*),konl(20),nea,neb,mi(*),mint2d,nopes,intscheme,
     &     nelcon(2,*),nrhcon(*),nalcon(2,*),ielmat(mi(3),*),nr,
     &     ielorien(mi(3),*),ntmat_,ipkon(*),ne0,iflag,null,kscale,
     &     istep,iinc,mt,ne,mattyp,ithermal(*),iprestr,i,j,k,m1,m2,jj,
     &     i1,m3,m4,kk,nener,indexe,nope,norien,iperturb(*),iout,
     &     nal,icmd,ihyper,nmethod,kode,imat,mint3d,iorien,ielas,
     &     istiff,ncmat_,nstate_,ikin,ilayer,nlayer,ki,kl,ielprop(*),
     &     nplicon(0:ntmat_,*),nplkcon(0:ntmat_,*),npmat_,calcul_fn,
     &     calcul_cauchy,calcul_qa,nopered,mortar,jfaces,igauss,
     &     istrainfree,nlgeom_undo,list,ilist(*),m,j1,mscalmethod,
     &     irowt(*),jqt(*),jqte(21),irowte(96),icmdcpy,length,id,
     &     islavquadel(*),node1,node2,j2,ii,mortartrafoflag,
     &     hasdamage,de12onepass,de12tangent,ndmat_,ndmcon(2,*),
     &     ivmap(21),jvmap(21),mattypp,ielasp,nlgeom_undop,de12fdskip
!     
      real*8 co(3,*),v(0:mi(2),*),shp(4,20),stiini(6,mi(1),*),
     &     stx(6,mi(1),*),xl(3,20),vl(0:mi(2),20),stre(6),prop(*),
     &     elcon(0:ncmat_,ntmat_,*),rhcon(0:1,ntmat_,*),xs2(3,7),
     &     alcon(0:6,ntmat_,*),vini(0:mi(2),*),thickness,
     &     alzero(*),orab(7,*),stiff(21),rho,fn(0:mi(2),*),
     &     fnl(3,10),skl(3,3),beta(6),q(0:mi(2),20),xl2(3,8),
     &     vkl(0:3,3),t0(*),t1(*),prestr(6,mi(1),*),eme(6,mi(1),*),
     &     ckl(3,3),vold(0:mi(2),*),eloc(6),veold(0:mi(2),*),
     &     springarea(2,*),elconloc(ncmat_),eth(6),xkl(3,3),
     &     voldl(0:mi(2),20),xikl(3,3),ener(2,mi(1),*),
     &     emec(6),eei(6,mi(1),*),enerini(2,mi(1),*),physcon(*),
     &     emec0(6),vel(1:3,20),veoldl(0:mi(2),20),xsj2(3),shp2(7,8),
     &     e,un,al,um,am1,xi,et,ze,tt,senergy,venergy,
     &     xsj,qa(*),vj,t0l,t1l,dtime,weight,pgauss(3),vij,time,ttime,
     &     plicon(0:2*npmat_,ntmat_,*),plkcon(0:2*npmat_,ntmat_,*),
     &     xstiff(27,mi(1),*),xstate(nstate_,mi(1),*),plconloc(802),
     &     vokl(3,3),xstateini(nstate_,mi(1),*),vikl(3,3),
     &     gs(8,4),a,reltime,tlayer(4),dlayer(4),xlayer(mi(3),4),
     &     thicke(mi(3),*),emeini(6,mi(1),*),clearini(3,9,*),
     &     pslavsurf(3,*),pmastsurf(6,*),smscale(*),sum1,sum2,
     &     scal,enerscal,elineng(6),t0g(2,*),t1g(2,*),aut(*),
     &     aute(96),shptil(4,20),xthi(3,3),vthj,
     &     dam(mi(1),*),damageD,damageg,damjac(12,mi(1),*),
     &     damvisc(mi(1),*),damviscini(mi(1),*),damvisceta,
     &     damvbeta,damvphys,damvbase,damvtrial,
     &     dmcon(0:ndmat_,ntmat_,*),dambase(mi(1),*),
     &     strebase(6),strep(6),stiffp(21),emecp(6),betap(6),
     &     damageq(6),damtrial0,damageDtrial,damageDp,damageh,
     &     depviscp,pnewdtp,vjp
      real*8, allocatable :: xstatesav(:)
!
!     residual stiffness floor, settable through CCX_DAMAGE_GMIN
!
      character*32 damgenv
      character*32 fdskipenv
      real*8 damggmin,spmean,damcbulk,damcmu,damcden,damcnumax
      real*8 damtanhi
      integer damgmintan,damcompress
      integer damginit
      save damggmin,damginit,damgmintan,damcompress,damcnumax
      save damtanhi
      data damginit /0/
      data ivmap /1,1,2,1,2,3,1,2,3,4,1,2,3,4,5,1,2,3,4,5,6/
      data jvmap /1,2,2,3,3,3,4,4,4,4,5,5,5,5,5,6,6,6,6,6,6/
!     
      include "gauss.f"
!
      if(damginit.eq.0) then
        damginit=1
        damggmin=1.d-4
        damgmintan=1
        damcompress=0
        call getenv('CCX_DAMAGE_COMPRESSION',damgenv)
        if(damgenv(1:1).eq.'1') damcompress=1
!
!       Preserving the bulk modulus while the shear modulus vanishes drives
!       nu' -> 0.5, and a linear tetrahedron LOCKS there.  The cap bounds how
!       incompressible a fully damaged element may become; it is the price of
!       using C3D4 and is settable so the trade can be measured.
!
        damcnumax=0.45d0
        call getenv('CCX_DAMAGE_COMPRESSION_NUMAX',damgenv)
        if(damgenv(1:1).ne.' ') then
          read(damgenv,*,err=7200,end=7200) damcnumax
          if((damcnumax.lt.0.d0).or.(damcnumax.gt.0.4999d0))
     &         damcnumax=0.45d0
        endif
 7200   continue
        call getenv('CCX_DAMAGE_GMIN_TANGENT',damgenv)
        if(damgenv(1:1).eq.'0') damgmintan=0
        call getenv('CCX_DAMAGE_GMIN',damgenv)
        if(damgenv(1:1).ne.' ') then
          read(damgenv,*,err=7100,end=7100) damggmin
          if((damggmin.lt.1.d-6).or.(damggmin.gt.0.2d0)) damggmin=1.d-4
        endif
 7100   continue
!
!       Upper cut-off of the consistent tangent.  1.999 is D<0.999, the
!       SAME number as the terminal-deletion threshold - two thresholds
!       that mean different things were given one value.  Deletion is
!       batched (64 per pass) and transactional, so an element can sit
!       at D>=0.999 for whole increments carrying g*C_ep alone, and
!       g*C_ep is measured at rho=3443 with the wrong sign on the
!       dominant entry already at D=0.994 (benchmarks/mpfd).
!
!       With the flag the cut-off moves to the RESIDUAL-STIFFNESS FLOOR,
!       which is the boundary that actually matters: on the floor g is
!       pinned at gmin, dg/dD=0 and the rank-1 term genuinely vanishes;
!       below it, it does not.
!
!       THIS CHANGES THE OPERATOR AND THEREFORE THE ANSWER.  Default off
!       = bit-identical.
!
        damtanhi=1.999d0
        call getenv('CCX_DAMAGE_TANGENT_FULL',damgenv)
        if(damgenv(1:1).eq.'1') damtanhi=2.d0-damggmin
      endif
!
      iflag=3
      null=0
!     
      mt=mi(2)+1
      nal=0
      qa(3)=-1.d0
      qa(4)=0.d0
      enerscal=0.d0
      if(((de12tangent.eq.1).or.(de12tangent.eq.2)).and.
     &     (nstate_.gt.0)) then
        allocate(xstatesav(nstate_))
      endif
!
!     CCX_DAMAGE_FD_SKIP=1 - CHANGES THE OPERATOR, despite reading as a
!     diagnostic below.  With it on the tangent is the SECANT one, so it is
!     not safe to sweep into a "diagnostics only" bucket.  Skips the
!     forward-difference
!     construction of dD/d(eps), leaving damageq at zero.  The rank-1 term
!     then vanishes, which is the same matrix CCX_DAMAGE_UNSYM_SCALE=0
!     produces - so running the two together gives an IDENTICAL Newton path
!     and the wall-time difference is exactly the cost of the six extra
!     mechmodel calls per actively softening integration point.  It is a
!     measurement of the prize an analytic derivative could win, not a
!     usable mode: with it on, the tangent is the secant one.
!
      de12fdskip=0
      call getenv('CCX_DAMAGE_FD_SKIP',fdskipenv)
      if(fdskipenv(1:1).eq.'1') de12fdskip=1
!     
      do m=nea,neb
!     
        if(list.eq.1) then
          i=ilist(m)
        else
          i=m
        endif
!     
!     check for strainless reactivated elements
!     
        if(ielmat(1,i).lt.0) then
          istrainfree=1
          ielmat(1,i)=-ielmat(1,i)
        else
          istrainfree=0
        endif
!     
        lakonl=lakon(i)
!     
!     only structural elements (no fluid elements)
!     
        if((ipkon(i).lt.0).or.(lakonl(1:1).eq.'F')) cycle
!     
        if(lakonl(1:7).eq.'DCOUP3D') cycle
!     
!     user elements
!     
        if(lakonl(1:1).eq.'U') then
          call resultsmech_u(co,kon,ipkon,lakon,ne,v,
     &         stx,elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,
     &         ielmat,ielorien,norien,orab,ntmat_,t0,t1,ithermal,prestr,
     &         iprestr,eme,iperturb,fn,iout,qa,vold,nmethod,
     &         veold,dtime,time,ttime,plicon,nplicon,plkcon,nplkcon,
     &         xstateini,xstiff,xstate,npmat_,matname,mi,ielas,icmd,
     &         ncmat_,nstate_,stiini,vini,ener,eei,enerini,istep,iinc,
     &         reltime,calcul_fn,calcul_qa,calcul_cauchy,nener,
     &         ikin,nal,ne0,thicke,emeini,i,ielprop,prop,t0g,t1g)
          cycle
        endif
!     
        if(lakonl(7:8).ne.'LC') then
!     
          imat=ielmat(1,i)
          amat=matname(imat)
          if(norien.gt.0) then
            iorien=max(0,ielorien(1,i))
          else
            iorien=0
          endif
!     
          if(nelcon(1,imat).lt.0) then
            ihyper=1
          else
            ihyper=0
          endif
        elseif(lakonl(4:5).eq.'20') then
!     
!     composite materials    
!     S8R - composite element
!     
          mint2d=4
          nopes=8
!     
!     determining the number of layers
!     
          nlayer=0
          do k=1,mi(3)
            if(ielmat(k,i).ne.0) then
              nlayer=nlayer+1
            endif
          enddo
!     
!     determining the layer thickness and global thickness
!     at the shell integration points
!     
          iflag=1
          indexe=ipkon(i)
          do kk=1,mint2d
            xi=gauss3d2(1,kk)
            et=gauss3d2(2,kk)
            call shape8q(xi,et,xl2,xsj2,xs2,shp2,iflag)
            tlayer(kk)=0.d0
            do k=1,nlayer
              thickness=0.d0
              do j=1,nopes
                thickness=thickness+thicke(k,indexe+j)*shp2(4,j)
              enddo
              tlayer(kk)=tlayer(kk)+thickness
              xlayer(k,kk)=thickness
            enddo
          enddo
          iflag=3
!     
          ilayer=0
          do k=1,4
            dlayer(k)=0.d0
          enddo
!     
!     S6 - composite element
!     
        elseif(lakonl(4:5).eq.'15') then
          mint2d=3
          nopes=6
!     determining the number of layers
!     
          nlayer=0
          do k=1,mi(3)
            if(ielmat(k,i).ne.0) then
              nlayer=nlayer+1
            endif
          enddo
!     
!     determining the layer thickness and global thickness
!     at the shell integration points
!     
          iflag=1
          indexe=ipkon(i)
          do kk=1,mint2d
            xi=gauss3d10(1,kk)
            et=gauss3d10(2,kk)
            call shape6tri(xi,et,xl2,xsj2,xs2,shp2,iflag)
            tlayer(kk)=0.d0
            do k=1,nlayer
              thickness=0.d0
              do j=1,nopes
                thickness=thickness+thicke(k,indexe+j)*shp2(4,j)
              enddo
              tlayer(kk)=tlayer(kk)+thickness
              xlayer(k,kk)=thickness
            enddo
          enddo
          iflag=3
!     
          ilayer=0
          do k=1,3
            dlayer(k)=0.d0
          enddo
!     
        endif
!     
        indexe=ipkon(i)
c     Bernhardi start
        if(lakonl(1:5).eq.'C3D8I') then
          nope=11
        elseif(lakonl(4:5).eq.'20') then
c     Bernhardi end
          nope=20
        elseif(lakonl(4:4).eq.'8') then
          nope=8
        elseif(lakonl(4:5).eq.'10') then
          nope=10
        elseif(lakonl(4:4).eq.'4') then
          nope=4
        elseif(lakonl(4:5).eq.'15') then
          nope=15
        elseif(lakonl(4:4).eq.'6') then
          nope=6
        elseif((lakonl(1:1).eq.'E').and.(lakonl(7:7).ne.'F')) then
!     
!     spring elements, contact spring elements and
!     dashpot elements
!     
          if(lakonl(7:7).eq.'C') then
!     
!     contact spring elements
!     
            if(mortar.eq.1) then
!     
!     face-to-face penalty
!     
              nope=kon(ipkon(i))
            elseif(mortar.eq.0) then
!     
!     node-to-face penalty
!     
              nope=ichar(lakonl(8:8))-47
              konl(nope+1)=kon(indexe+nope+1)
            endif
          else
!     
!     genuine spring elements and dashpot elements
!     
            nope=ichar(lakonl(8:8))-47
          endif
        else
          cycle
        endif
!     
        if(lakonl(4:5).eq.'8R') then
          if(intscheme.eq.0) then
            mint3d=1
          else
            mint3d=8
          endif
        elseif(lakonl(4:7).eq.'20RB') then
          if((lakonl(8:8).eq.'R').or.(lakonl(8:8).eq.'C')) then
            mint3d=50
          else
            call beamintscheme(lakonl,mint3d,ielprop(i),prop,
     &           null,xi,et,ze,weight)
          endif
        elseif((lakonl(4:4).eq.'8').or.
     &         (lakonl(4:6).eq.'20R')) then
          if(lakonl(7:8).eq.'LC') then
            mint3d=8*nlayer
          else
            if((intscheme.eq.0).or.(lakonl(4:8).ne.'20R  ')) then
              mint3d=8
            else
              mint3d=27
            endif
          endif
        elseif(lakonl(4:4).eq.'2') then
          mint3d=27
        elseif(lakonl(4:5).eq.'10') then
          mint3d=4
        elseif(lakonl(4:4).eq.'4') then
          mint3d=1
        elseif(lakonl(4:5).eq.'15') then
          if(lakonl(7:8).eq.'LC') then
            mint3d=6*nlayer
          else
            mint3d=9
          endif
        elseif(lakonl(4:4).eq.'6') then
          mint3d=2
        elseif(lakonl(1:1).eq.'E') then
          mint3d=0
        endif
!     
        do j=1,nope
          konl(j)=kon(indexe+j)
          do k=1,3
            xl(k,j)=co(k,konl(j))
            vl(k,j)=v(k,konl(j))
            voldl(k,j)=vold(k,konl(j))
          enddo
        enddo
!     
!     mortar start
!     
!     calculating the transformation matrix for a quadratic element containing
!     at least one slave node; this matrix transforms the regular 
!     quadratic shape functions into purely positive ones for slave
!     faces.    
!     
        if(mortartrafoflag.gt.0) then
          if(islavquadel(i).gt.0) then
              jqte(1)=1
              ii=1
              do i1=1,nope
                node1=konl(i1)
                length=jqt(node1+1)-jqt(node1)
                do j2=1,nope
                  node2=konl(j2)
                  call nident(irowt(jqt(node1)),node2,length,id)
                  if(id.gt.0) then
                    j1=jqt(node1)+id-1
                    if(irowt(j1).eq.node2) then
                      aute(ii)=aut(j1)
                      irowte(ii)=j2
                      ii=ii+1
                    endif
                  endif
                enddo
                jqte(i1+1)=ii
              enddo
          endif
        endif
!     
!     mortar end
!     
!     q contains the nodal forces per element; initialisation of q
!     
        if((iperturb(1).ge.2).or.((iperturb(1).le.0).and.(iout.lt.1))) 
     &       then
          do m1=1,nope
            do m2=0,mi(2)
              q(m2,m1)=fn(m2,konl(m1))
            enddo
          enddo
        endif
!     
!     calculating the forces for the contact elements
!     
        if(mint3d.eq.0) then
!     
!     "normal" spring and dashpot elements
!     
          kode=nelcon(1,imat)
          if((lakonl(7:7).eq.'A').or.(lakonl(7:7).eq.'1').or.
     &         (lakonl(7:7).eq.'2')) then
            t0l=0.d0
            t1l=0.d0
            if(ithermal(1).eq.1) then
              t0l=(t0(konl(1))+t0(konl(2)))/2.d0
              t1l=(t1(konl(1))+t1(konl(2)))/2.d0
            elseif(ithermal(1).ge.2) then
              t0l=(t0(konl(1))+t0(konl(2)))/2.d0
              t1l=(vold(0,konl(1))+vold(0,konl(2)))/2.d0
            endif
          endif
!     
!     spring elements (including contact springs)
!     
          if(lakonl(2:2).eq.'S') then
            if((lakonl(7:7).eq.'A').or.(lakonl(7:7).eq.'1').or.
     &           (lakonl(7:7).eq.'2').or.((mortar.eq.0).and.
     &           ((nmethod.ne.1).or.(iperturb(1).ge.2).or.
     &           (iout.ne.-1)))) then
!
!     nr is the element number in the assumption that all
!     slave nodes lead to contact, i.e. that in each slave node a
!     spring element was generated in gencontelem_n2f.f
!
              if(nener.eq.1) then
                if(lakonl(7:7).ne.'C') then
                  nr=i
                else
                  nr=ne0+konl(nope+1)
                endif
c                venergy=enerini(2,1,nr)
                venergy=0.d0
              endif
              call springforc_n2f(xl,konl,vl,imat,elcon,nelcon,stiff,
     &             fnl,ncmat_,ntmat_,nope,lakonl,t1l,kode,elconloc,
     &             plicon,nplicon,npmat_,senergy,nener,
     &             stx(1,1,i),mi,springarea(1,konl(nope+1)),nmethod,
     &             ne0,nstate_,xstateini,xstate,reltime,
     &             ielas,venergy,ielorien,orab,norien,i,
     &             smscale,mscalmethod)
              if(nener.eq.1) then
                ener(1,1,nr)=senergy
                ener(2,1,nr)=venergy
              endif
            elseif((mortar.eq.1).and.
     &             ((nmethod.ne.1).or.(iperturb(1).ge.2).or.
     &             (iout.ne.-1))) then
              jfaces=kon(indexe+nope+2)
              igauss=kon(indexe+nope+1)
c              if(nener.eq.1) venergy=enerini(2,1,ne0+igauss)
              if(nener.eq.1) venergy=0.d0
              call springforc_f2f(xl,vl,imat,elcon,nelcon,stiff,
     &             fnl,ncmat_,ntmat_,nope,lakonl,t1l,kode,elconloc,
     &             plicon,nplicon,npmat_,senergy,nener,
     &             stx(1,1,i),mi,springarea(1,igauss),nmethod,
     &             ne0,nstate_,xstateini,xstate,reltime,
     &             ielas,jfaces,igauss,pslavsurf,pmastsurf,
     &             clearini,venergy,kscale,konl,iout,i,
     &             smscale,mscalmethod)
              if(nener.eq.1) then
                ener(1,1,ne0+igauss)=senergy
                ener(2,1,ne0+igauss)=venergy
              endif
            endif
!     
!     next lines are not executed in linstatic.c before the
!     setup of the stiffness matrix (i.e. nmethod=1 and
!     iperturb(1)<1 and iout=-1).
!     
            if((lakonl(7:7).eq.'A').or.
     &           ((nmethod.ne.1).or.(iperturb(1).ge.2).or.(iout.ne.-1)))
     &           then
              do j=1,nope
                do k=1,3
                  fn(k,konl(j))=fn(k,konl(j))+fnl(k,j)
                enddo
              enddo
            endif
          endif
        elseif(ikin.eq.1) then
          do j=1,nope
            do k=1,3
              veoldl(k,j)=veold(k,konl(j))
            enddo
          enddo            
        endif
!     
        do jj=1,mint3d
          if(lakonl(4:5).eq.'8R') then
            if(intscheme.eq.0) then
              xi=gauss3d1(1,jj)
              et=gauss3d1(2,jj)
              ze=gauss3d1(3,jj)
              weight=weight3d1(jj)
            else
!
!     initial acceleration calculation for dynamic calculations:
!     full integration used
!
              xi=gauss3d2(1,jj)
              et=gauss3d2(2,jj)
              ze=gauss3d2(3,jj)
              weight=weight3d2(jj)
            endif
          elseif(lakonl(4:7).eq.'20RB') then
            if((lakonl(8:8).eq.'R').or.(lakonl(8:8).eq.'C')) then
              xi=gauss3d13(1,jj)
              et=gauss3d13(2,jj)
              ze=gauss3d13(3,jj)
              weight=weight3d13(jj)
            else
              call beamintscheme(lakonl,mint3d,ielprop(i),prop,
     &             jj,xi,et,ze,weight)
            endif
          elseif((lakonl(4:4).eq.'8').or.
     &           (lakonl(4:6).eq.'20R'))
     &           then
            if(lakonl(7:8).ne.'LC') then
              if((intscheme.eq.0).or.(lakonl(4:8).ne.'20R  ')) then
                xi=gauss3d2(1,jj)
                et=gauss3d2(2,jj)
                ze=gauss3d2(3,jj)
                weight=weight3d2(jj)
              else
!
!     initial acceleration calculation for dynamic calculations:
!     full integration used
!
                xi=gauss3d3(1,jj)
                et=gauss3d3(2,jj)
                ze=gauss3d3(3,jj)
                weight=weight3d3(jj)
              endif
            else
              kl=mod(jj,8)
              if(kl.eq.0) kl=8
!     
              xi=gauss3d2(1,kl)
              et=gauss3d2(2,kl)
              ze=gauss3d2(3,kl)
              weight=weight3d2(kl)
!     
              ki=mod(jj,4)
              if(ki.eq.0) ki=4
!     
              if(kl.eq.1) then
                ilayer=ilayer+1
                if(ilayer.gt.1) then
                  do k=1,4
                    dlayer(k)=dlayer(k)+xlayer(ilayer-1,k)
                  enddo
                endif
              endif
              ze=2.d0*(dlayer(ki)+(ze+1.d0)/2.d0*xlayer(ilayer,ki))/
     &             tlayer(ki)-1.d0
              weight=weight*xlayer(ilayer,ki)/tlayer(ki)
!     
!     material and orientation
!     
              imat=ielmat(ilayer,i)
              amat=matname(imat)
              if(norien.gt.0) then
                iorien=max(0,ielorien(ilayer,i))
              else
                iorien=0
              endif
!     
              if(nelcon(1,imat).lt.0) then
                ihyper=1
              else
                ihyper=0
              endif
            endif
          elseif(lakonl(4:4).eq.'2') then
            xi=gauss3d3(1,jj)
            et=gauss3d3(2,jj)
            ze=gauss3d3(3,jj)
            weight=weight3d3(jj)
          elseif(lakonl(4:5).eq.'10') then
            xi=gauss3d5(1,jj)
            et=gauss3d5(2,jj)
            ze=gauss3d5(3,jj)
            weight=weight3d5(jj)
          elseif(lakonl(4:4).eq.'4') then
            xi=gauss3d4(1,jj)
            et=gauss3d4(2,jj)
            ze=gauss3d4(3,jj)
            weight=weight3d4(jj)
          elseif(lakonl(4:5).eq.'15') then
            if(lakonl(7:8).ne.'LC') then
              xi=gauss3d8(1,jj)
              et=gauss3d8(2,jj)
              ze=gauss3d8(3,jj)
              weight=weight3d8(jj)
            else
              kl=mod(jj,6)
              if(kl.eq.0) kl=6
!     
              xi=gauss3d10(1,kl)
              et=gauss3d10(2,kl)
              ze=gauss3d10(3,kl)
              weight=weight3d10(kl)
!     
              ki=mod(jj,3)
              if(ki.eq.0) ki=3
!     
              if(kl.eq.1) then
                ilayer=ilayer+1
                if(ilayer.gt.1) then
                  do k=1,3
                    dlayer(k)=dlayer(k)+xlayer(ilayer-1,k)
                  enddo
                endif
              endif
              ze=2.d0*(dlayer(ki)+(ze+1.d0)/2.d0*xlayer(ilayer,ki))/
     &             tlayer(ki)-1.d0
              weight=weight*xlayer(ilayer,ki)/tlayer(ki)
!     
!     material and orientation
!     
              imat=ielmat(ilayer,i)
              amat=matname(imat)
              if(norien.gt.0) then
                iorien=max(0,ielorien(ilayer,i))
              else
                iorien=0
              endif
!     
              if(nelcon(1,imat).lt.0) then
                ihyper=1
              else
                ihyper=0
              endif
            endif
          elseif(lakonl(4:4).eq.'6') then
            xi=gauss3d7(1,jj)
            et=gauss3d7(2,jj)
            ze=gauss3d7(3,jj)
            weight=weight3d7(jj)
          endif
!     
c     Bernhardi start
          if(lakonl(1:5).eq.'C3D8R') then
            call shape8hr(xl,xsj,shp,gs,a)
            icmdcpy=icmd
!     
!           tangent stiffness ds/de is needed in hgforce.f    
!     
            icmd=0
          elseif(lakonl(1:5).eq.'C3D8I') then
            call shape8hu(xi,et,ze,xl,xsj,shp,iflag)
          elseif(nope.eq.20) then
c     Bernhardi end
            if(lakonl(7:7).eq.'A') then
              call shape20h_ax(xi,et,ze,xl,xsj,shp,iflag)
            elseif((lakonl(7:7).eq.'E').or.
     &             (lakonl(7:7).eq.'S')) then
              call shape20h_pl(xi,et,ze,xl,xsj,shp,iflag)
            else
              call shape20h(xi,et,ze,xl,xsj,shp,iflag)
            endif
          elseif(nope.eq.8) then
            call shape8h(xi,et,ze,xl,xsj,shp,iflag)
          elseif(nope.eq.10) then
            call shape10tet(xi,et,ze,xl,xsj,shp,iflag)
          elseif(nope.eq.4) then
            call shape4tet(xi,et,ze,xl,xsj,shp,iflag)
          elseif(nope.eq.15) then
            call shape15w(xi,et,ze,xl,xsj,shp,iflag)
          else
            call shape6w(xi,et,ze,xl,xsj,shp,iflag)
          endif
!     
!     mortar start
!     transforming the shape functions for quadratic elements containing at    
!     lease one slave node into purely positive functions on the slave    
!     faces     
!     
          if(mortartrafoflag.gt.0) then
            if(islavquadel(i).gt.0) then
                do i1=1,nope
                  if(jqte(i1+1)-jqte(i1).gt.0) then
                    shptil(1,i1)=0.0
                    shptil(2,i1)=0.0
                    shptil(3,i1)=0.0
                    shptil(4,i1)=0.0
                  else
                    shptil(1,i1)=shp(1,i1)
                    shptil(2,i1)=shp(2,i1)
                    shptil(3,i1)=shp(3,i1)
                    shptil(4,i1)=shp(4,i1)
                  endif
                  do j1=jqte(i1),jqte(i1+1)-1
                    j2=irowte(j1)
                    shptil(1,i1)=shptil(1,i1)+aute(j1)
     &                   *shp(1,j2)
                    shptil(2,i1)=shptil(2,i1)+aute(j1)
     &                   *shp(2,j2)
                    shptil(3,i1)=shptil(3,i1)+aute(j1)
     &                   *shp(3,j2)
                    shptil(4,i1)=shptil(4,i1)+aute(j1)
     &                   *shp(4,j2)
                  enddo
                enddo
            else
              do i1=1,nope
                shptil(1,i1)=shp(1,i1)
                shptil(2,i1)=shp(2,i1)
                shptil(3,i1)=shp(3,i1)
                shptil(4,i1)=shp(4,i1)
              enddo
            endif
          endif
!
!           mortar end
!     
!     vkl(m2,m3) contains the derivative of the m2-
!     component of the displacement with respect to
!     direction m3
!     
          do m2=1,3
            do m3=1,3
              vkl(m2,m3)=0.d0
            enddo
          enddo
!     
          do m1=1,nope
            do m2=1,3
              do m3=1,3
                vkl(m2,m3)=vkl(m2,m3)+shp(m3,m1)*vl(m2,m1)
              enddo
            enddo
          enddo
!     
!     for frequency analysis or buckling with preload the
!     strains are calculated with respect to the deformed
!     configuration
!     for a linear iteration within a nonlinear increment:
!     the tangent matrix is calculated at strain at the end
!     of the previous increment
!     
          if((iperturb(1).eq.1).or.(iperturb(1).eq.-1)) then
            do m2=1,3
              do m3=1,3
                vokl(m2,m3)=0.d0
              enddo
            enddo
!     
            do m1=1,nope
              do m2=1,3
                do m3=1,3
                  vokl(m2,m3)=vokl(m2,m3)+
     &                 shp(m3,m1)*voldl(m2,m1)
                enddo
              enddo
            enddo
          endif
!     
          kode=nelcon(1,imat)
!     
!     calculating the strain
!     
          call calctotstrain(vkl,vokl,eloc,elineng,iperturb)
!     
!     calculating the deformation gradient (needed to
!     convert the element stiffness matrix from spatial
!     coordinates to material coordinates
!     deformation plasticity)
!     
          if((kode.eq.-50).or.(kode.eq.-53).or.(kode.eq.-54).or.
     &         (kode.le.-100)) then
!     
!     calculating the deformation gradient
!     
c     Bernhardi start
            xkl(1,1)=vkl(1,1)+1.0d0
            xkl(2,2)=vkl(2,2)+1.0d0
            xkl(3,3)=vkl(3,3)+1.0d0
c     Bernhardi end
            xkl(1,2)=vkl(1,2)
            xkl(1,3)=vkl(1,3)
            xkl(2,3)=vkl(2,3)
            xkl(2,1)=vkl(2,1)
            xkl(3,1)=vkl(3,1)
            xkl(3,2)=vkl(3,2)
!     
!     calculating the Jacobian
!     
            vj=xkl(1,1)*(xkl(2,2)*xkl(3,3)-xkl(2,3)*xkl(3,2))
     &           -xkl(1,2)*(xkl(2,1)*xkl(3,3)-xkl(2,3)*xkl(3,1))
     &           +xkl(1,3)*(xkl(2,1)*xkl(3,2)-xkl(2,2)*xkl(3,1))
!     
          endif
!     
!     calculating fields for incremental plasticity
!     
          if(kode.le.-100) then
!     
!     calculating the deformation gradient at the
!     start of the increment
!     
!     calculating the displacement gradient at the
!     start of the increment
!     
            do m2=1,3
              do m3=1,3
                vikl(m2,m3)=0.d0
              enddo
            enddo
!     
            do m1=1,nope
              do m2=1,3
                do m3=1,3
                  vikl(m2,m3)=vikl(m2,m3)
     &                 +shp(m3,m1)*vini(m2,konl(m1))
                enddo
              enddo
            enddo
!     
!     calculating the deformation gradient of the old
!     fields
!     
            xikl(1,1)=vikl(1,1)+1.d0
            xikl(2,2)=vikl(2,2)+1.d0
            xikl(3,3)=vikl(3,3)+1.d0
            xikl(1,2)=vikl(1,2)
            xikl(1,3)=vikl(1,3)
            xikl(2,3)=vikl(2,3)
            xikl(2,1)=vikl(2,1)
            xikl(3,1)=vikl(3,1)
            xikl(3,2)=vikl(3,2)
!     
!     calculating the Jacobian
!     
            vij=xikl(1,1)*(xikl(2,2)*xikl(3,3)
     &           -xikl(2,3)*xikl(3,2))
     &           -xikl(1,2)*(xikl(2,1)*xikl(3,3)
     &           -xikl(2,3)*xikl(3,1))
     &           +xikl(1,3)*(xikl(2,1)*xikl(3,2)
     &           -xikl(2,2)*xikl(3,1))
!     
!     stresses at the start of the increment
!     
            do m1=1,6
              stre(m1)=stiini(m1,jj,i)
            enddo
!     
          endif
!     
!     prestress values
!     
          if(iprestr.eq.1) then
            do kk=1,6
              beta(kk)=-prestr(kk,jj,i)
            enddo
          else
            do kk=1,6
              beta(kk)=0.d0
            enddo
          endif
!     
          if(ithermal(1).ge.1) then
!     
!     calculating the temperature difference in
!     the integration point
!     
            t0l=0.d0
            t1l=0.d0
            if(ithermal(1).eq.1) then
              if((lakonl(4:5).eq.'8 ').or.
     &             (lakonl(4:5).eq.'8I')) then
                do i1=1,8
                  t0l=t0l+t0(konl(i1))/8.d0
                  t1l=t1l+t1(konl(i1))/8.d0
                enddo
              elseif(lakonl(4:6).eq.'20 ') then
                nopered=20
                call lintemp(t0,konl,nopered,jj,t0l)
                call lintemp(t1,konl,nopered,jj,t1l)
              elseif(lakonl(4:6).eq.'10T') then
                call linscal10(t0,konl,t0l,null,shp)
                call linscal10(t1,konl,t1l,null,shp)
              else
                do i1=1,nope
                  t0l=t0l+shp(4,i1)*t0(konl(i1))
                  t1l=t1l+shp(4,i1)*t1(konl(i1))
                enddo
              endif
            elseif(ithermal(1).ge.2) then
              if((lakonl(4:5).eq.'8 ').or.
     &             (lakonl(4:5).eq.'8I')) then
                do i1=1,8
                  t0l=t0l+t0(konl(i1))/8.d0
                  t1l=t1l+vold(0,konl(i1))/8.d0
                enddo
              elseif(lakonl(4:6).eq.'20 ') then
                nopered=20
                call lintemp_th0(t0,konl,nopered,jj,t0l,mi)
                call lintemp_th1(vold,konl,nopered,jj,t1l,mi)
              elseif(lakonl(4:6).eq.'10T') then
                call linscal10(t0,konl,t0l,null,shp)
                call linscal10(vold,konl,t1l,mi(2),shp)
              else
                do i1=1,nope
                  t0l=t0l+shp(4,i1)*t0(konl(i1))
                  t1l=t1l+shp(4,i1)*vold(0,konl(i1))
                enddo
              endif
            endif
            tt=t1l-t0l
          endif
!     
!     calculating the coordinates of the integration point
!     for material orientation purposes (for cylindrical
!     coordinate systems)
!     
          if((iorien.gt.0).or.(kode.le.-100)) then
            do j=1,3
              pgauss(j)=0.d0
              do i1=1,nope
                pgauss(j)=pgauss(j)+shp(4,i1)*co(j,konl(i1))
              enddo
            enddo
          endif
!     
!     material data; for linear elastic materials
!     this includes the calculation of the stiffness
!     matrix
!     
          istiff=0
!     
          call materialdata_me(elcon,nelcon,rhcon,nrhcon,alcon,
     &         nalcon,imat,amat,iorien,pgauss,orab,ntmat_,
     &         stiff,rho,i,ithermal,alzero,mattyp,t0l,t1l,ihyper,
     &         istiff,elconloc,eth,kode,plicon,nplicon,
     &         plkcon,nplkcon,npmat_,plconloc,mi(1),dtime,jj,
     &         xstiff,ncmat_,iperturb)
!     
!     determining the mechanical strain
!     
          if(ithermal(1).ne.0) then
            call calcmechstrain(vkl,vokl,emec,eth,iperturb,nalcon,imat,
     &           xthi,vthj)
          else
            do m1=1,6
              emec(m1)=eloc(m1)
            enddo
          endif
c          write(*,*) 'resultsmech1 ',i,jj,(emec(m1),m1=1,6)
          if(kode.le.-100) then
            do m1=1,6
              emec0(m1)=emeini(m1,jj,i)
            enddo
          endif
!     
!     subtracting the initial strains
!     
          if(iprestr.eq.2) then
            if(istrainfree.eq.0) then
              do m1=1,6
                emec(m1)=emec(m1)-prestr(m1,jj,i)
              enddo
            else
              do m1=1,6
                prestr(m1,jj,i)=emec(m1)
                emeini(m1,jj,i)=emec(m1)
                eme(m1,jj,i)=emec(m1)
                emec(m1)=0.d0
              enddo
            endif
          endif
!     
!     calculating the local stiffness and stress
!     
          nlgeom_undo=0
          call mechmodel(elconloc,stiff,emec,kode,emec0,ithermal,
     &         icmd,beta,stre,xkl,ckl,vj,xikl,vij,
     &         plconloc,xstate,xstateini,ielas,
     &         amat,t1l,dtime,time,ttime,i,jj,nstate_,mi(1),
     &         iorien,pgauss,orab,eloc,mattyp,qa(3),istep,iinc,
     &         ipkon,nmethod,iperturb,qa(4),nlgeom_undo,physcon,
     &         ncmat_,nalcon,imat)
c          write(*,*) 'resultsmech2 ',i,jj,(stre(m1),m1=1,6)
c          write(*,*) 'resultsmech3 ',i,jj,(stiff(m1),m1=1,21)
!
!     modifying the stress and stiffness for a multiplicative
!     decomposition of the deformation gradient in a mechanical and
!     a thermal part
!
          if((ithermal(1).ne.0).and.(iperturb(2).eq.1)) then
            call modifystressstiff(stre,stiff,mattyp,eth,nalcon,imat,
     &     xthi,vthj)
          endif
c          write(*,*) 'resultsmech4 ',i,jj,(stre(m1),m1=1,6)
!
!         BK1 one-pass progressive-damage update.  stre is still the
!         effective constitutive stress here.  damageupdatepoint rebuilds
!         dam(jj,i) from dambase and the trial PEEQ generated by mechmodel.
!
          if(de12onepass.ne.0) then
            call damageupdatepoint(i,jj,ipkon,lakon,kon,co,mi,
     &           ielmat,ne0,ndmat_,ntmat_,ndmcon,dmcon,dam,dambase,
     &           stre,xstate,xstateini,nstate_)
          endif
!
!         BK2 tangent experiment.  Only the scalar derivative dD/d(emec)
!         is evaluated by forward differences; the verified incplas
!         algorithmic tangent remains the effective-stress Jacobian.  The
!         exact outer product is generally nonsymmetric.  xstiff/e_c3d can
!         currently carry only 21 constitutive entries, so this opt-in stage
!         uses its symmetric part.  All point state and damage are restored
!         after every perturbation.
!
          do m1=1,6
            strebase(m1)=stre(m1)
            damageq(m1)=0.d0
          enddo
          damtrial0=0.d0
          if((hasdamage.ne.0).and.(i.le.ne0)) then
            damtrial0=dam(jj,i)
          endif
!
!         UNSYM stage 2.  damjac carries, per integration point, the
!         effective stress and dD/d(strain) that the asymmetric second
!         assembly pass turns into the rank-1 correction
!
!             C = g*C_ep - sigma_eff (x) dD/d(eps)
!
!         It must be cleared whenever the point is not actively
!         softening, because it survives between Newton iterations.
!
          if((de12tangent.eq.2).and.(i.le.ne0)) then
            do m1=1,12
              damjac(m1,jj,i)=0.d0
            enddo
          endif
          if(((de12tangent.eq.1).or.(de12tangent.eq.2)).and.
     &       (de12fdskip.eq.0).and.
     &       (de12onepass.ne.0).and.
     &       (icmd.ne.3).and.(nstate_.gt.0).and.
     &       (mattyp.eq.3).and.(damtrial0.gt.1.d0+1.d-12).and.
     &       (damtrial0.lt.damtanhi).and.
     &       (damtrial0-dambase(jj,i).gt.1.d-14)) then
            damageDtrial=damtrial0-1.d0
            do m1=1,nstate_
              xstatesav(m1)=xstate(m1,jj,i)
            enddo
            do k=1,6
              do m1=1,6
                emecp(m1)=emec(m1)
                betap(m1)=beta(m1)
              enddo
              damageh=1.d-7*dmax1(1.d0,dabs(emec(k)))
              emecp(k)=emecp(k)+damageh
              do m1=1,nstate_
                xstate(m1,jj,i)=xstateini(m1,jj,i)
              enddo
              mattypp=mattyp
              ielasp=ielas
              nlgeom_undop=0
              depviscp=qa(3)
              pnewdtp=qa(4)
              vjp=vj
              call mechmodel(elconloc,stiffp,emecp,kode,emec0,
     &             ithermal,icmd,betap,strep,xkl,ckl,vjp,xikl,vij,
     &             plconloc,xstate,xstateini,ielasp,amat,t1l,dtime,
     &             time,ttime,i,jj,nstate_,mi(1),iorien,pgauss,orab,
     &             eloc,mattypp,depviscp,istep,iinc,ipkon,nmethod,
     &             iperturb,pnewdtp,nlgeom_undop,physcon,ncmat_,
     &             nalcon,imat)
              if((ithermal(1).ne.0).and.(iperturb(2).eq.1)) then
                call modifystressstiff(strep,stiffp,mattypp,eth,
     &               nalcon,imat,xthi,vthj)
              endif
              call damageupdatepoint(i,jj,ipkon,lakon,kon,co,mi,
     &             ielmat,ne0,ndmat_,ntmat_,ndmcon,dmcon,dam,
     &             dambase,strep,xstate,xstateini,nstate_)
              damageDp=dam(jj,i)-1.d0
              if(damageDp.lt.0.d0) damageDp=0.d0
              if(damageDp.gt.1.d0) damageDp=1.d0
              damageq(k)=(damageDp-damageDtrial)/damageh
            enddo
            do m1=1,nstate_
              xstate(m1,jj,i)=xstatesav(m1)
            enddo
            dam(jj,i)=damtrial0
!
!           Store sigma_eff and dD/d(eps) for the asymmetric pass.
!           damageq holds dD/d(emec).  incplas builds the right
!           Cauchy-Green tensor as c(4)=2*emec(4), so emec(4:6) are
!           TENSOR shears while the constitutive matrix that e_c3d
!           consumes is expressed in ENGINEERING shear.  The shear
!           components therefore need the factor 1/2 to live in the
!           same convention as the 21-entry tangent.
!
!           With viscous regularisation the degradation follows Dvis, so
!           the consistent tangent needs d(Dvis)/d(eps) = beta*dD/d(eps).
!
            if((de12tangent.eq.2).and.(i.le.ne0)) then
              damvbeta=1.d0
              if((damvisceta.gt.0.d0).and.(dtime.gt.0.d0)) then
                damvbeta=dtime/(damvisceta+dtime)
              endif
              do m1=1,6
                damjac(m1,jj,i)=strebase(m1)
              enddo
              do m1=1,3
                damjac(6+m1,jj,i)=damvbeta*damageq(m1)
              enddo
              do m1=4,6
                damjac(6+m1,jj,i)=0.5d0*damvbeta*damageq(m1)
              enddo
            endif
          endif
!
!         DE1 scalar stiffness/stress degradation.  dam is a combined
!         state coordinate: omega_D before initiation and 1+D after
!         initiation.  Keeping this operation outside mechmodel/incplas
!         leaves the stock plastic return mapping untouched.
!
!         The tangent used in DE1 is deliberately the symmetric secant-like
!         approximation g(D)*Cep.  The full damage-consistent derivative is
!         deferred to a later patch so the symmetric PARDISO path remains
!         available for this first implementation.
!
!         R4 viscous regularisation.  A brittle phase failing inside a
!         ductile one snaps through at a structural limit point, which no
!         tangent and no convergence criterion can follow: grip
!         displacement control does not control the local event because
!         the surrounding matrix acts as a soft spring in series.
!
!         Relaxing the degradation variable
!
!             Dvis = Dvis_committed + beta*(D - Dvis_committed)
!             beta = dtime/(eta+dtime)
!
!         makes the softening rate dependent, which removes the limit
!         point.  D itself keeps evolving from the committed baseline at
!         the physical rate, so the evolution law and the dissipated
!         fracture energy are untouched; only the stiffness the assembly
!         sees is delayed.  eta=0 reproduces the inviscid branch exactly.
!
!         Dvis is monotone: it can never fall below its committed value,
!         so unloading cannot heal a point.
!
          damvbeta=1.d0
          if((damvisceta.gt.0.d0).and.(dtime.gt.0.d0)) then
            damvbeta=dtime/(damvisceta+dtime)
          endif
          if((hasdamage.ne.0).and.(i.le.ne0).and.
     &         (damvisceta.gt.0.d0)) then
            damvphys=dam(jj,i)-1.d0
            if(damvphys.lt.0.d0) damvphys=0.d0
            if(damvphys.gt.1.d0) damvphys=1.d0
            damvbase=damviscini(jj,i)
            if(damvbase.lt.0.d0) damvbase=0.d0
            if(damvbase.gt.1.d0) damvbase=1.d0
            damvtrial=damvbase+damvbeta*(damvphys-damvbase)
            if(damvtrial.lt.damvbase) damvtrial=damvbase
            if(damvtrial.gt.1.d0) damvtrial=1.d0
            damvisc(jj,i)=damvtrial
          endif
!
          if((hasdamage.ne.0).and.(i.le.ne0)) then
            if(dam(jj,i).gt.1.d0) then
              damageD=dam(jj,i)-1.d0
              if(damageD.lt.0.d0) damageD=0.d0
              if(damageD.gt.1.d0) damageD=1.d0
              if(damvisceta.gt.0.d0) damageD=damvisc(jj,i)
              damageg=1.d0-damageD
!
!             Residual stiffness floor.  A fully damaged element keeps
!             gmin of its tangent until the terminal scan removes it.  At
!             1e-4 a few such elements are harmless, but a notch process
!             zone holds thousands of them at once and the operator
!             becomes badly scaled: that is the regime where the DHC runs
!             stall.  Raising the floor trades a small unreleased load for
!             conditioning; the element is deleted at D>=0.999 either way,
!             so the floor only acts transiently.
!
              if(damageg.lt.damggmin) then
                damageg=damggmin
!
!               On the floor g is constant in D, so d(sigma)/d(eps)
!               is gmin*Cep and the -sigma_eff (x) dD/d(eps) term is
!               identically zero.  damjac is written earlier, before
!               it is known whether the floor will engage, so leaving
!               it populated makes mafilldamas assemble a rank-1
!               correction that does not exist.
!
!               Found with the structural FD probe: at SP1 increment
!               95 one node of the probed element carried a relative
!               Jacobian error of 9.4e-3 while its three element-mates
!               sat at 4e-6, and the error was identical to seven
!               digits over h from 1e-8 to 1e-6, so it was not finite
!               difference truncation.  That node is the one that also
!               touches a floored element.
!
                if((de12tangent.eq.2).and.(i.le.ne0).and.
     &               (damgmintan.eq.1)) then
                  do m1=1,12
                    damjac(m1,jj,i)=0.d0
                  enddo
                endif
              endif
!             Compression is not degraded (CCX_DAMAGE_COMPRESSION=1).
!
!             Abaqus Verification Guide 2.2.21-V states the reference
!             behaviour of ductile damage: the deviatoric and the TENSILE
!             hydrostatic response degrade, the COMPRESSIVE hydrostatic
!             response does not, because a cracked element still carries
!             compression once the crack closes.  cohesive_uc6.f:236 already
!             does exactly this at the facet - "compression retains the full
!             normal penalty" - and the bulk did not, which made the two
!             laws in this tree inconsistent with each other (E-92).
!
!             The tangent is split to match.  For isotropic elasticity
!             C = K*(1(x)1) + 2*mu*Idev, so degrading only the deviator
!             means K' = K and mu' = g*mu, which in (E,nu) form is
!                 E'  = 9*K*g*mu/(3*K+g*mu)
!                 nu' = (3*K-2*g*mu)/(2*(3*K+g*mu))
!             At g=1 this returns (E,nu) identically.  As g->0 it tends to
!             an incompressible solid with no shear stiffness, which is the
!             physically right limit but ill-conditioned in (E,nu) form, so
!             nu' is capped.  Only mattyp=1 is split; anisotropic cards keep
!             the uniform scaling, because the deviatoric projector is not
!             expressible by scaling their constants.
!
              if(damcompress.eq.1) then
                spmean=(stre(1)+stre(2)+stre(3))/3.d0
              else
                spmean=1.d0
              endif
              if((damcompress.eq.1).and.(spmean.lt.0.d0)) then
                do m1=1,3
                  stre(m1)=damageg*(stre(m1)-spmean)+spmean
                enddo
                do m1=4,6
                  stre(m1)=damageg*stre(m1)
                enddo
                if(mattyp.eq.1) then
                  damcbulk=stiff(1)/(3.d0*(1.d0-2.d0*stiff(2)))
                  damcmu=damageg*stiff(1)/(2.d0*(1.d0+stiff(2)))
                  damcden=3.d0*damcbulk+damcmu
                  if(damcden.gt.1.d-30) then
                    stiff(1)=9.d0*damcbulk*damcmu/damcden
                    stiff(2)=(3.d0*damcbulk-2.d0*damcmu)/(2.d0*damcden)
                    if(stiff(2).gt.damcnumax) stiff(2)=damcnumax
                    if(stiff(2).lt.-1.d0) stiff(2)=-1.d0
                  else
                    stiff(1)=damageg*stiff(1)
                  endif
                elseif(mattyp.eq.2) then
                  do m1=1,9
                    stiff(m1)=damageg*stiff(m1)
                  enddo
                else
                  do m1=1,21
                    stiff(m1)=damageg*stiff(m1)
                  enddo
                  if(de12tangent.eq.1) then
                    do m1=1,21
                      stiff(m1)=stiff(m1)-0.5d0*
     &                     (strebase(ivmap(m1))*damageq(jvmap(m1))+
     &                      strebase(jvmap(m1))*damageq(ivmap(m1)))
                    enddo
                  endif
                endif
              else
              do m1=1,6
                stre(m1)=damageg*stre(m1)
              enddo
              if(mattyp.eq.1) then
!               isotropic elastic representation: stiff(1)=E,
!               stiff(2)=nu.  Scale E only; Poisson's ratio is not a
!               stiffness magnitude.
                stiff(1)=damageg*stiff(1)
              elseif(mattyp.eq.2) then
                do m1=1,9
                  stiff(m1)=damageg*stiff(m1)
                enddo
              else
                do m1=1,21
                  stiff(m1)=damageg*stiff(m1)
                enddo
                if(de12tangent.eq.1) then
                  do m1=1,21
                    stiff(m1)=stiff(m1)-0.5d0*
     &                   (strebase(ivmap(m1))*damageq(jvmap(m1))+
     &                    strebase(jvmap(m1))*damageq(ivmap(m1)))
                  enddo
                endif
              endif
              endif
            endif
          endif
!     
          if(((nmethod.ne.4).or.(iperturb(1).ne.0)).and.
     &         (nmethod.ne.5).and.(icmd.ne.3)) then
            do m1=1,21
              xstiff(m1,jj,i)=stiff(m1)
            enddo
          endif
!     
          if((iperturb(1).eq.-1).and.(nlgeom_undo.eq.0)) then
!     
!     if the forced displacements were changed at
!     the start of a nonlinear step, the nodal
!     forces due do this displacements are 
!     calculated in a purely linear way, and
!     the first iteration is purely linear in order
!     to allow the displacements to redistribute
!     in a quasi-static way (only applies to
!     quasi-static analyses (*STATIC))
!     
            eloc(1)=elineng(1)-vokl(1,1)
            eloc(2)=elineng(2)-vokl(2,2)
            eloc(3)=elineng(3)-vokl(3,3)
            eloc(4)=elineng(4)-(vokl(1,2)+vokl(2,1))
            eloc(5)=elineng(5)-(vokl(1,3)+vokl(3,1))
            eloc(6)=elineng(6)-(vokl(2,3)+vokl(3,2))
!     
            if(mattyp.eq.1) then
              e=stiff(1)
              un=stiff(2)
              um=e/(1.d0+un)
              al=un*um/(1.d0-2.d0*un)
              um=um/2.d0
              am1=al*(eloc(1)+eloc(2)+eloc(3))
              stre(1)=am1+2.d0*um*eloc(1)
              stre(2)=am1+2.d0*um*eloc(2)
              stre(3)=am1+2.d0*um*eloc(3)
              stre(4)=um*eloc(4)
              stre(5)=um*eloc(5)
              stre(6)=um*eloc(6)
            elseif(mattyp.eq.2) then
              stre(1)=eloc(1)*stiff(1)+eloc(2)*stiff(2)
     &             +eloc(3)*stiff(4)
              stre(2)=eloc(1)*stiff(2)+eloc(2)*stiff(3)
     &             +eloc(3)*stiff(5)
              stre(3)=eloc(1)*stiff(4)+eloc(2)*stiff(5)
     &             +eloc(3)*stiff(6)
              stre(4)=eloc(4)*stiff(7)
              stre(5)=eloc(5)*stiff(8)
              stre(6)=eloc(6)*stiff(9)
            elseif(mattyp.eq.3) then
              stre(1)=eloc(1)*stiff(1)+eloc(2)*stiff(2)+
     &             eloc(3)*stiff(4)+eloc(4)*stiff(7)+
     &             eloc(5)*stiff(11)+eloc(6)*stiff(16)
              stre(2)=eloc(1)*stiff(2)+eloc(2)*stiff(3)+
     &             eloc(3)*stiff(5)+eloc(4)*stiff(8)+
     &             eloc(5)*stiff(12)+eloc(6)*stiff(17)
              stre(3)=eloc(1)*stiff(4)+eloc(2)*stiff(5)+
     &             eloc(3)*stiff(6)+eloc(4)*stiff(9)+
     &             eloc(5)*stiff(13)+eloc(6)*stiff(18)
              stre(4)=eloc(1)*stiff(7)+eloc(2)*stiff(8)+
     &             eloc(3)*stiff(9)+eloc(4)*stiff(10)+
     &             eloc(5)*stiff(14)+eloc(6)*stiff(19)
              stre(5)=eloc(1)*stiff(11)+eloc(2)*stiff(12)+
     &             eloc(3)*stiff(13)+eloc(4)*stiff(14)+
     &             eloc(5)*stiff(15)+eloc(6)*stiff(20)
              stre(6)=eloc(1)*stiff(16)+eloc(2)*stiff(17)+
     &             eloc(3)*stiff(18)+eloc(4)*stiff(19)+
     &             eloc(5)*stiff(20)+eloc(6)*stiff(21)
            endif
          endif
!     
!     updating the internal energy and mechanical strain
!     for user materials (kode<=-100) the mechanical strain has to
!     be updated at the end of each increment (also if no output
!     is requested), since it is input to the umat routine
!     
c          if((iout.ge.0).or.(iout.eq.-2).or.(kode.le.-100).or.
          if((iout.gt.0).or.(iout.eq.-2).or.(kode.le.-100).or.
!
!              next line was inserted since the extra call
!              of results at the end of an increment is not 
!              executed if no printing information is needed
!              for that increment (due to the FREQUENCY parameter)          
!              07.07.2023
!
     &         ((iout.eq.0).and.(nener.eq.1)).or.
     &         ((nmethod.eq.4).and.
     &         ((iperturb(1).gt.1).and.(nlgeom_undo.eq.0)).and.
     &         (ithermal(1).le.1))) then
!     
            if(nener.eq.1) then
!     sigma.d(epsilon)=(sigma_ini+sigma)*(epsilon-epsilon_ini)/2
              ener(1,jj,i)=enerini(1,jj,i)+
     &             ((emec(1)-emeini(1,jj,i))*
     &             (stre(1)+stiini(1,jj,i))+
     &             (emec(2)-emeini(2,jj,i))*
     &             (stre(2)+stiini(2,jj,i))+
     &             (emec(3)-emeini(3,jj,i))*
     &             (stre(3)+stiini(3,jj,i)))/2.d0+
     &             (emec(4)-emeini(4,jj,i))*(stre(4)+stiini(4,jj,i))+
     &             (emec(5)-emeini(5,jj,i))*(stre(5)+stiini(5,jj,i))+
     &             (emec(6)-emeini(6,jj,i))*(stre(6)+stiini(6,jj,i))
            endif
            eme(1,jj,i)=emec(1)
            eme(2,jj,i)=emec(2)
            eme(3,jj,i)=emec(3)
            eme(4,jj,i)=emec(4)
            eme(5,jj,i)=emec(5)
            eme(6,jj,i)=emec(6)
          endif
!     
          if(iout.gt.0) then
!     
            eei(1,jj,i)=eloc(1)
            eei(2,jj,i)=eloc(2)
            eei(3,jj,i)=eloc(3)
            eei(4,jj,i)=eloc(4)
            eei(5,jj,i)=eloc(5)
            eei(6,jj,i)=eloc(6)
          endif
!     
!     updating the kinetic energy
!     
          if(ikin.eq.1) then
!     
            call materialdata_rho(rhcon,nrhcon,imat,rho,t1l,
     &           ntmat_,ithermal)
            do m1=1,3
              vel(m1,1)=0.d0
              do i1= 1,nope
                vel(m1,1)=vel(m1,1)+shp(4,i1)*veoldl(m1,i1)
              enddo
            enddo
            ener(2,jj,i)=rho*(vel(1,1)*vel(1,1)+
     &           vel(2,1)*vel(2,1)+ vel(3,1)*vel(3,1))/2.d0
          endif
!     
          skl(1,1)=stre(1)
          skl(2,2)=stre(2)
          skl(3,3)=stre(3)
          skl(2,1)=stre(4)
          skl(3,1)=stre(5)
          skl(3,2)=stre(6)
!     
          stx(1,jj,i)=skl(1,1)
          stx(2,jj,i)=skl(2,2)
          stx(3,jj,i)=skl(3,3)
          stx(4,jj,i)=skl(2,1)
          stx(5,jj,i)=skl(3,1)
          stx(6,jj,i)=skl(3,2)
!     
          skl(1,2)=skl(2,1)
          skl(1,3)=skl(3,1)
          skl(2,3)=skl(3,2)
!     
!     calculation of the nodal forces
!     
          if(calcul_fn.eq.1) then
!     
!     calculating fn using skl
!     
!     mortar start
!     
            if(mortartrafoflag.gt.0) then
!
!             using the tilde shape functions
!
              do m1=1,nope
                do m2=1,3
!
!                          linear elastic part
!
                  do m3=1,3
                    fn(m2,konl(m1))=fn(m2,konl(m1))+
     &                   xsj*skl(m2,m3)*shptil(m3,m1)*weight
                  enddo
!
!                          nonlinear geometric part
!
                  if(iperturb(2).eq.1) then
                    do m3=1,3
                      do m4=1,3
                        fn(m2,konl(m1))=fn(m2,konl(m1))+
     &                       xsj*skl(m4,m3)*weight*
     &                       (vkl(m2,m4)*shptil(m3,m1)+
     &                       vkl(m2,m3)*shptil(m4,m1))/2.d0
                      enddo
                    enddo
                  endif
!
                enddo
              enddo
!     
!     mortar end
!     
            else
              do m1=1,nope
                do m2=1,3
!     
!     linear elastic part
!     
                  do m3=1,3
                    fn(m2,konl(m1))=fn(m2,konl(m1))+
     &                   xsj*skl(m2,m3)*shp(m3,m1)*weight
                  enddo
!     
!     nonlinear geometric part
!     
                  if(iperturb(2).eq.1) then
                    do m3=1,3
                      do m4=1,3
                        fn(m2,konl(m1))=fn(m2,konl(m1))+
     &                       xsj*skl(m4,m3)*weight*
     &                       (vkl(m2,m4)*shp(m3,m1)+
     &                       vkl(m2,m3)*shp(m4,m1))/2.d0
                      enddo
                    enddo
                  endif
!     
                enddo
              enddo
            endif
c     Bernhardi start
            if(lakonl(1:5).eq.'C3D8R') then
              call hgforce(fn,stiff,a,gs,vl,mi,konl)
              icmd=icmdcpy
            endif
c     Bernhardi end
          endif
!     
!     calculation of the Cauchy stresses
!     
          if(calcul_cauchy.eq.1) then
!     
!     changing the displacement gradients into
!     deformation gradients
!     
            if((kode.ne.-50).and.(kode.gt.-100)) then
c     Bernhardi start
              xkl(1,1)=vkl(1,1)+1.0d0
              xkl(2,2)=vkl(2,2)+1.0d0
              xkl(3,3)=vkl(3,3)+1.0d0
c     Bernhardi end   
              xkl(1,2)=vkl(1,2)
              xkl(1,3)=vkl(1,3)
              xkl(2,3)=vkl(2,3)
              xkl(2,1)=vkl(2,1)
              xkl(3,1)=vkl(3,1)
              xkl(3,2)=vkl(3,2)
!     
              vj=xkl(1,1)*(xkl(2,2)*xkl(3,3)-xkl(2,3)*xkl(3,2))
     &             -xkl(1,2)*(xkl(2,1)*xkl(3,3)-xkl(2,3)*xkl(3,1))
     &             +xkl(1,3)*(xkl(2,1)*xkl(3,2)-xkl(2,2)*xkl(3,1))
            endif
!     
            do m1=1,3
              do m2=1,m1
                ckl(m1,m2)=0.d0
                do m3=1,3
                  do m4=1,3
                    ckl(m1,m2)=ckl(m1,m2)+
     &                   skl(m3,m4)*xkl(m1,m3)*xkl(m2,m4)
                  enddo
                enddo
                ckl(m1,m2)=ckl(m1,m2)/vj
              enddo
            enddo
!     
            stx(1,jj,i)=ckl(1,1)
            stx(2,jj,i)=ckl(2,2)
            stx(3,jj,i)=ckl(3,3)
            stx(4,jj,i)=ckl(2,1)
            stx(5,jj,i)=ckl(3,1)
            stx(6,jj,i)=ckl(3,2)
          endif
!     
        enddo
!     
!     q contains the contributions to the nodal force in the nodes
!     belonging to the element at stake from other elements (elements
!     already treated). These contributions have to be
!     subtracted to get the contributions attributable to the element
!     at stake only
!     
        if(calcul_qa.eq.1) then
          do m1=1,nope
            do m2=1,3
              qa(1)=qa(1)+dabs(fn(m2,konl(m1))-q(m2,m1))
            enddo
          enddo
          nal=nal+3*nope
        endif
        
!     Calculation of additional kinetic energy if selevtive mass scaling 
!     is used: E_kin_add= v^T *lamda *v *0.5
!     
        if((ikin.eq.1).and.
     &       ((mscalmethod.eq.1).or.(mscalmethod.eq.3))) then
          do m1=1,3
            sum2=0.d0
            do i1=1,nope
              sum1=0.d0
              do j1=1,nope
                if(i1.eq.j1) then
                  scal=smscale(i)*(nope-1)
                else
                  scal=-1*smscale(i)
                endif
                sum1=sum1+(veoldl(m1,j1))*scal
              enddo
              sum2=sum2+sum1*veoldl(m1,i1)
            enddo
            enerscal=enerscal+sum2
          enddo
        endif
!     
      enddo
!     
c          if(j.ne.-1) stop
      if(allocated(xstatesav)) deallocate(xstatesav)
      return
      end
