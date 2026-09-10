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
      subroutine damconnectsets(ipkon,kon,lakon,ne,nk,set,nset,
     &     istartset,iendset,ialset,namea,nameb,iconn,nreach,imode,
     &     ifacdead)
!
!     Load-path test between two named node sets over the surviving
!     elements.  Resolves the set names, then hands the node lists to
!     damconnect (imode=0, two elements conduct through a shared NODE)
!     or damconnectface (imode=1, they must share a whole FACE).
!
!     imode=1 is what CCX_FRACTURE_LINK=FACE selects.  It exists because
!     a specimen can finish in two pieces joined by a single vertex,
!     which the node rule reads as intact (E-75).  It deletes nothing -
!     only the moment the run may stop changes.
!
!     CalculiX stores node sets with a trailing 'N' in set(), so the
!     lookup appends it and compares case-insensitively; the deck may
!     write FACE_X0_NSET while the internal name is FACE_X0_NSETN.
!
!     A negative entry in ialset encodes a generated range a,b,-c meaning
!     a, a+c, ... up to b, which is expanded here.
!
      implicit none
!
      character*8 lakon(*)
      character*81 set(*)
!
!     fixed length, not character*(*): the C side calls through the
!     FORTRAN macro and passes no hidden string lengths, so an assumed
!     length argument reads garbage
!
      character*81 namea,nameb
!
      integer ipkon(*),kon(*),ne,nk,nset,istartset(*),iendset(*),
     &     ialset(*),iconn,nreach,imode,ifacdead(*),ia,ib,na,nb
!
      integer, allocatable :: nodesa(:),nodesb(:)
!
      iconn=1
      nreach=0
!
      call damsetfind(set,nset,namea,ia)
      call damsetfind(set,nset,nameb,ib)
      if((ia.le.0).or.(ib.le.0)) then
        write(*,*) '*WARNING in damconnectsets: node set not found;'
        write(*,*) '         fracture termination is inactive.'
        if(ia.le.0) write(*,*) '         missing: ',namea
        if(ib.le.0) write(*,*) '         missing: ',nameb
        return
      endif
!
      call damsetcount(istartset,iendset,ialset,ia,na)
      call damsetcount(istartset,iendset,ialset,ib,nb)
      if((na.le.0).or.(nb.le.0)) return
!
      allocate(nodesa(na))
      allocate(nodesb(nb))
      call damsetfill(istartset,iendset,ialset,ia,nodesa,na)
      call damsetfill(istartset,iendset,ialset,ib,nodesb,nb)
!
      if(imode.eq.1) then
        call damconnectface(ipkon,kon,lakon,ne,nk,nodesa,na,nodesb,nb,
     &       iconn,nreach,ifacdead)
        deallocate(nodesb,nodesa)
        return
      endif
      call damconnect(ipkon,kon,lakon,ne,nk,nodesa,na,nodesb,nb,
     &     iconn,nreach,ifacdead)
!
      deallocate(nodesb,nodesa)
      return
      end
!
!     ------------------------------------------------------------------
!
      subroutine damsetfind(set,nset,name,id)
      implicit none
      character*81 set(*)
      character*81 name
      character*81 target
      integer nset,id,i,j,n
!
      id=0
      n=len_trim(name)
      if((n.le.0).or.(n.ge.81)) return
      target=' '
      do i=1,n
        target(i:i)=name(i:i)
        if((target(i:i).ge.'a').and.(target(i:i).le.'z')) then
          target(i:i)=char(ichar(target(i:i))-32)
        endif
      enddo
      target(n+1:n+1)='N'
!
      do i=1,nset
        do j=1,81
          if(set(i)(j:j).ne.target(j:j)) goto 10
        enddo
        id=i
        return
 10     continue
      enddo
      return
      end
!
!     ------------------------------------------------------------------
!
      subroutine damsetcount(istartset,iendset,ialset,id,n)
      implicit none
      integer istartset(*),iendset(*),ialset(*),id,n,j,k
!
      n=0
      j=istartset(id)
      do
        if(j.gt.iendset(id)) exit
        if(ialset(j).gt.0) then
          n=n+1
          j=j+1
        else
!
!         generated range: ialset(j-2), ialset(j-1), -ialset(j)
!
          do k=ialset(j-2),ialset(j-1),-ialset(j)
            n=n+1
          enddo
          n=n-2
          j=j+1
        endif
      enddo
      return
      end
!
!     ------------------------------------------------------------------
!
      subroutine damsetfill(istartset,iendset,ialset,id,nodes,n)
      implicit none
      integer istartset(*),iendset(*),ialset(*),id,nodes(*),n,j,k,m
!
      m=0
      j=istartset(id)
      do
        if(j.gt.iendset(id)) exit
        if(ialset(j).gt.0) then
          if(m.lt.n) then
            m=m+1
            nodes(m)=ialset(j)
          endif
          j=j+1
        else
          m=m-2
          do k=ialset(j-2),ialset(j-1),-ialset(j)
            if(m.lt.n) then
              m=m+1
              nodes(m)=k
            endif
          enddo
          j=j+1
        endif
      enddo
      n=m
      return
      end
!
!     ------------------------------------------------------------------
!
      subroutine damsetnodes(set,nset,istartset,iendset,ialset,name,
     &     nodes,n)
!
!     Expose a resolved node set to the C side, so loadpath.c can own the
!     endpoints without a second copy of the set lookup.  Resolution is the
!     existing, in-use damsetfind/damsetcount/damsetfill; nothing new is
!     decided here.
!
!     Called twice: once with n as the only output to size the array (the
!     caller passes a null pointer for nodes, so it must not be touched when
!     n comes in as 0), then again to fill it.
!
      implicit none
!
      character*81 set(*)
      character*81 name
      integer nset,istartset(*),iendset(*),ialset(*),nodes(*),n,
     &     id,na,nwant
!
      nwant=n
      n=0
      call damsetfind(set,nset,name,id)
      if(id.le.0) return
      call damsetcount(istartset,iendset,ialset,id,na)
      if(na.le.0) return
      n=na
      if(nwant.le.0) return
      if(na.gt.nwant) n=nwant
      call damsetfill(istartset,iendset,ialset,id,nodes,n)
      return
      end
