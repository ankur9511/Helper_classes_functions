subroutine mednc_bond_angle(npairsc11, npairsc12, npairsc21, npairsc22, &
        & npairsc12c22, nshellc11, nshellc12, nshellc21, nshellc22, &
        & nshell1c1c2, npairsc1c2, &
        & centers1,ncenters1, &
        & shell,nshell, &
        & centers2,ncenters2, &
        & r,criteria,tot_atom,box,binv,dist2,cosangle_crit, &
        & Hcenters1,Hcenters2,Hshell)
implicit none

!!!!! If angle criteria is present, center must contain all O's that have a bonded H
character*100, intent(in) :: criteria

integer, intent(in) :: tot_atom,ncenters1,ncenters2,nshell
real(4), intent(in), dimension(0:tot_atom-1,0:2) :: r
real(8), intent(in), dimension(0:2) :: box,binv

integer, intent(in), dimension(0:ncenters1-1) :: centers1
integer, intent(in), dimension(0:nshell-1) :: shell
integer, intent(in), dimension(0:ncenters2-1) :: centers2
integer, dimension(0:nshell-1) :: shellc11, shellc12, shellc21, shellc22, notshellc11, notshellc21
real(4), intent(in) :: dist2 ! xl, xh, yl, yh, zl, zh

integer, optional, intent(in), dimension(0:ncenters1-1,0:1) :: Hcenters1
integer, optional, intent(in), dimension(0:ncenters2-1,0:1) :: Hcenters2
integer, optional, intent(in), dimension(0:nshell-1,0:1) :: Hshell

!!!!! If angle criteria is present, H must contain atomID of H's bonded to centers. Max is 2. If only 1 present, atomID is -1
real(4), optional, intent(in) :: cosangle_crit
real(4), dimension(0:2) :: dr
real(4) :: dr2,costheta
integer :: i,j,k,index1,index2,nnshellc11,nnshellc21,flag
integer, intent(inout) :: npairsc11, npairsc12, npairsc21, npairsc22, npairsc12c22, nshellc11, nshellc12, nshellc21, nshellc22, nshell1c1c2, npairsc1c2

!!print *, "Starting initialization"
nshellc11 = 0
nshellc12 = 0
nshellc21 = 0
nshellc22 = 0
nnshellc11 = 0
nnshellc21 = 0

shellc11(:) = -1
shellc12(:) = -1
shellc21(:) = -1
shellc22(:) = -1
notshellc11(:) = -1
notshellc21(:) = -1

npairsc11 = 0
npairsc12 = 0 
npairsc21 = 0 
npairsc22 = 0 
npairsc12c22 = 0
nshell1c1c2 = 0
npairsc1c2 = 0

index1=index(criteria,'dist')
index2=index(criteria,'angle')

!if (index1 * index2 > 0) then
!!print *, "Compute first shell for centers1",ncenters1,nshell
!!!!! Compute first shell for centers1
  do j = 0,ncenters1-1 !! j = Water or any other O in center
        do i = 0,nshell-1    !! i = Water in tetrahedral shell
        flag = -1
        !!print *, "Loop index",i,j,nshellc11,Hcenters1(j,:),Hshell(i,:)
        if (centers1(j)/=shell(i)) then
                call pbcdr(dr,r(centers1(j)-1,:),r(shell(i)-1,:),box,binv)
                dr2 = dr(0)*dr(0) + dr(1)*dr(1) + dr(2)*dr(2)
                costheta = 0.0
!                
!                !!!! Independent "if" logics used because a center atom can donate
!                !multiple h-bonds to same shell atom - and all of them should be calculated as
!                ! h-bonds
!                !print *, "dr2 calc done",i,j,nshellc11,Hcenters1(j,:),Hshell(i,:)
                if (Hcenters1(j,0)/=-1) then
                        costheta = 0.0
                        call cosangle(costheta,r(centers1(j)-1,:),r(Hcenters1(j,0)-1,:),r(shell(i)-1,:),box,binv)

                        if (dr2 <= dist2 .and. costheta >= cosangle_crit) then
                                npairsc11 = npairsc11+1
                                flag = 1
                        end if
                end if
!                !print *, "c11-1 done"
!
                if (Hcenters1(j,1)/=-1) then
                        costheta = 0.0
                        call cosangle(costheta,r(centers1(j)-1,:),r(Hcenters1(j,1)-1,:),r(shell(i)-1,:),box,binv)

                        if (dr2 <= dist2 .and. costheta >= cosangle_crit) then
                                npairsc11 = npairsc11+1
                                flag = 1
                        end if
                end if
!
!                !!!! Independent "if" logics used because a center atom can accept
!                !multiple h-bonds same shell - and all of them should be calculated as
!                ! h-bonds
!                 !print *, "c11-2 done"
!               
                if (Hshell(i,0)/=-1) then
                        costheta = 0.0
                        call cosangle(costheta,r(shell(i)-1,:),r(Hshell(i,0)-1,:),r(centers1(j)-1,:),box,binv)

                        if (dr2 <= dist2 .and. costheta >= cosangle_crit) then
                                npairsc11 = npairsc11+1
                                flag = 1
                        end if
                end if
!                !print *, "c11-3 done"
!
                if (Hshell(i,1)/=-1) then
                        costheta = 0.0
                        call cosangle(costheta,r(shell(i)-1,:),r(Hshell(i,1)-1,:),r(centers1(j)-1,:),box,binv)

                        if (dr2 <= dist2 .and. costheta >= cosangle_crit) then
                                npairsc11 = npairsc11+1
                                flag = 1
                        end if
                end if
!
!                !print *, "shell calc done"
!
                if (flag == 1) then
                        if ( .NOT. ANY(shellc11==i)) then
                        shellc11(nshellc11) = i
                        nshellc11 = nshellc11+1
                        end if
                end if
                if (j == ncenters1-1) then
                        if ( .NOT. ANY(shellc11==i)) then
                        notshellc11(nnshellc11) = i
                        nnshellc11 = nnshellc11+1
                        end if
                end if
                !print *, nshellc11,Hcenters1(j,:),Hshell(i,:)
        end if
        end do
  end do
!print *, nshellc11,nnshellc11,nshellc11+nnshellc11,nshell
!  !print *, "Compute second shell for centers1"
!  !!!! Compute second shell of centers1
  do j = 0,nshellc11-1 !! j = Water or any other O in center
        !print *, "Inside 2nd shell calc",j,nshellc11-1
        do i = 0,nnshellc11-1    !! i = Water in tetrahedral shell
        flag = -1
!        !if (shell(shellc11(j))/=shell(i)) then
                call pbcdr(dr,r(shell(shellc11(j))-1,:),r(shell(notshellc11(i))-1,:),box,binv)
                dr2 = dr(0)*dr(0) + dr(1)*dr(1) + dr(2)*dr(2)
                costheta = 0.0
!                
!                !!!! Independent "if" logics used because a center atom can donate
!                !multiple h-bonds to same shell atom - and all of them should be calculated as
!                ! h-bonds
!                
                if (Hshell(shellc11(j),0)/=-1) then
                        costheta = 0.0
                        call cosangle(costheta,r(shell(shellc11(j))-1,:),r(Hshell(shellc11(j),0)-1,:),r(shell(notshellc11(i))-1,:),box,binv)

                        if (dr2 <= dist2 .and. costheta >= cosangle_crit) then
                                npairsc12 = npairsc12+1
                                flag = 1
                        end if
                end if
!
                if (Hshell(shellc11(j),1)/=-1) then
                        costheta = 0.0
                        call cosangle(costheta,r(shell(shellc11(j))-1,:),r(Hshell(shellc11(j),1)-1,:),r(shell(notshellc11(i))-1,:),box,binv)

                        if (dr2 <= dist2 .and. costheta >= cosangle_crit) then
                                npairsc12 = npairsc12+1
                                flag = 1
                        end if
                end if
!
!                !!!! Independent "if" logics used because a center atom can accept
!                !multiple h-bonds same shell - and all of them should be calculated as
!                ! h-bonds
!                
                if (Hshell(notshellc11(i),0)/=-1) then
                        costheta = 0.0
                        call cosangle(costheta,r(shell(notshellc11(i))-1,:),r(Hshell(notshellc11(i),0)-1,:),r(shell(shellc11(j))-1,:),box,binv)

                        if (dr2 <= dist2 .and. costheta >= cosangle_crit) then
                                npairsc12 = npairsc12+1
                                flag = 1
                        end if
                end if

                if (Hshell(notshellc11(i),1)/=-1) then
                        costheta = 0.0
                        call cosangle(costheta,r(shell(notshellc11(i))-1,:),r(Hshell(notshellc11(i),1)-1,:),r(shell(shellc11(j))-1,:),box,binv)

                        if (dr2 <= dist2 .and. costheta >= cosangle_crit) then
                                npairsc12 = npairsc12+1
                                flag = 1
                        end if
                end if

                if (flag == 1) then
                        if ( .NOT. ANY(shellc12==i)) then
                        shellc12(nshellc12) = i
                        nshellc12 = nshellc12+1
                        end if
                end if
 
!
!                if (flag == 1) then
!                        shellc12(nshellc12) = i
!                        nshellc12 = nshellc12+1
!                end if
!        !end if
        end do
  end do
!
!  !print *, "Compute first shell for centers2"
!  !!!! Compute first shell of centers2
  do j = 0,ncenters2-1 !! j = Water or any other O in center
        do i = 0,nshell-1    !! i = Water in tetrahedral shell
        flag = -1
        if (centers2(j)/=shell(i)) then
                call pbcdr(dr,r(centers2(j)-1,:),r(shell(i)-1,:),box,binv)
                dr2 = dr(0)*dr(0) + dr(1)*dr(1) + dr(2)*dr(2)
                costheta = 0.0
!                
!                !!!! Independent "if" logics used because a center atom can donate
!                !multiple h-bonds to same shell atom - and all of them should be calculated as
!                ! h-bonds
!                
                if (Hcenters2(j,0)/=-1) then
                        costheta = 0.0
                        call cosangle(costheta,r(centers2(j)-1,:),r(Hcenters2(j,0)-1,:),r(shell(i)-1,:),box,binv)
!
                        if (dr2 <= dist2 .and. costheta >= cosangle_crit) then
                                npairsc21 = npairsc21+1
                                flag = 1
                        end if
                end if
!
                if (Hcenters2(j,1)/=-1) then
                        costheta = 0.0
                        call cosangle(costheta,r(centers2(j)-1,:),r(Hcenters2(j,1)-1,:),r(shell(i)-1,:),box,binv)

                        if (dr2 <= dist2 .and. costheta >= cosangle_crit) then
                                npairsc21 = npairsc21+1
                                flag = 1
                        end if
                end if
!
!                !!!! Independent "if" logics used because a center atom can accept
!                !multiple h-bonds same shell - and all of them should be calculated as
!                ! h-bonds
!                
                if (Hshell(i,0)/=-1) then
                        costheta = 0.0
                        call cosangle(costheta,r(shell(i)-1,:),r(Hshell(i,0)-1,:),r(centers2(j)-1,:),box,binv)

                        if (dr2 <= dist2 .and. costheta >= cosangle_crit) then
                                npairsc21 = npairsc21+1
                                flag = 1
                        end if
                end if

                if (Hshell(i,1)/=-1) then
                        costheta = 0.0
                        call cosangle(costheta,r(shell(i)-1,:),r(Hshell(i,1)-1,:),r(centers2(j)-1,:),box,binv)

                        if (dr2 <= dist2 .and. costheta >= cosangle_crit) then
                                npairsc21 = npairsc21+1
                                flag = 1
                        end if
                end if
!
                if (flag == 1) then
                        if ( .NOT. ANY(shellc21==i)) then
                        shellc21(nshellc21) = i
                        nshellc21 = nshellc21+1
                        end if
                end if
                if (j == ncenters2-1) then
                        if ( .NOT. ANY(shellc21==i)) then
                        notshellc21(nnshellc21) = i
                        nnshellc21 = nnshellc21+1
                        end if
                end if

!                if (flag == 1) then
!                        shellc21(nshellc21) = shell(i)
!                        nshellc21 = nshellc21+1
!                else if (flag == -1) then
!                       notshellc21(nnshellc21) = i
!                       nnshellc21 = nnshellc21+1
!                end if
        end if
        end do
  end do
!
!  !print *, "Compute second shell for centers2"
!  !!!! Compute second shell of centers2
  do j = 0,nshellc21-1 !! j = Water or any other O in center
        do i = 0,nnshellc21-1    !! i = Water in tetrahedral shell
        flag = -1
!        !if (shell(shellc11(j))/=shell(i)) then
                call pbcdr(dr,r(shell(shellc21(j))-1,:),r(shell(notshellc21(i))-1,:),box,binv)
                dr2 = dr(0)*dr(0) + dr(1)*dr(1) + dr(2)*dr(2)
                costheta = 0.0
!                
!                !!!! Independent "if" logics used because a center atom can donate
!                !multiple h-bonds to same shell atom - and all of them should be calculated as
!                ! h-bonds
!                
                if (Hshell(shellc21(j),0)/=-1) then
                        costheta = 0.0
                        call cosangle(costheta,r(shell(shellc21(j))-1,:),r(Hshell(shellc21(j),0)-1,:),r(shell(notshellc21(i))-1,:),box,binv)
!
                        if (dr2 <= dist2 .and. costheta >= cosangle_crit) then
                                npairsc22 = npairsc22+1
                                flag = 1
                        end if
                end if
!
                if (Hshell(shellc21(j),1)/=-1) then
                        costheta = 0.0
                        call cosangle(costheta,r(shell(shellc21(j))-1,:),r(Hshell(shellc21(j),1)-1,:),r(shell(notshellc21(i))-1,:),box,binv)
!
                        if (dr2 <= dist2 .and. costheta >= cosangle_crit) then
                                npairsc22 = npairsc22+1
                                flag = 1
                        end if
                end if
!
!                !!!! Independent "if" logics used because a center atom can accept
!                !multiple h-bonds same shell - and all of them should be calculated as
!                ! h-bonds
!                
                if (Hshell(notshellc21(i),0)/=-1) then
                        costheta = 0.0
                        call cosangle(costheta,r(shell(notshellc21(i))-1,:),r(Hshell(notshellc21(i),0)-1,:),r(shell(shellc21(j))-1,:),box,binv)

                        if (dr2 <= dist2 .and. costheta >= cosangle_crit) then
                                npairsc22 = npairsc22+1
                                flag = 1
                        end if
                end if
!
                if (Hshell(notshellc11(i),1)/=-1) then
                        costheta = 0.0
                        call cosangle(costheta,r(shell(notshellc21(i))-1,:),r(Hshell(notshellc21(i),1)-1,:),r(shell(shellc21(j))-1,:),box,binv)

                        if (dr2 <= dist2 .and. costheta >= cosangle_crit) then
                                npairsc22 = npairsc22+1
                                flag = 1
                        end if
                end if
!
                if (flag == 1) then
                        if ( .NOT. ANY(shellc22==i)) then
                        shellc22(nshellc22) = i
                        nshellc22 = nshellc22+1
                        end if
                end if
 
!
!                if (flag == 1) then
!                        shellc22(nshellc22) = i
!                        nshellc22 = nshellc22+1
!                end if
!        !end if
        end do
  end do
!  !print *, "Compute common in second shells for centers1 and centers2"
  do j = 0,nshellc22-1
        do i = 0,nshellc12-1
        if (shell(shellc22(j)) == shell(shellc12(i))) then
                npairsc12c22 = npairsc12c22 + 1
        end if
        end do
  end do
  !npairsc12c22 = int(npairsc12c22/2)

  do j = 0,nshellc21-1
        do i = 0,nshellc11-1
        if (shell(shellc21(j)) == shell(shellc11(i))) then
                nshell1c1c2 = nshell1c1c2 + 1
        end if
        end do
  end do
  !npairsc12c22 = int(npairsc12c22/2)

  do j = 0,ncenters2-1 !! j = Water or any other O in center
        do i = 0,ncenters1-1    !! i = Water in tetrahedral shell
        flag = -1
        if (centers2(j)/=centers1(i)) then
                call pbcdr(dr,r(centers2(j)-1,:),r(centers1(i)-1,:),box,binv)
                dr2 = dr(0)*dr(0) + dr(1)*dr(1) + dr(2)*dr(2)
                costheta = 0.0
!                
!                !!!! Independent "if" logics used because a center atom can donate
!                !multiple h-bonds to same shell atom - and all of them should be calculated as
!                ! h-bonds
!                
                if (Hcenters2(j,0)/=-1) then
                        costheta = 0.0
                        call cosangle(costheta,r(centers2(j)-1,:),r(Hcenters2(j,0)-1,:),r(centers1(i)-1,:),box,binv)
!
                        if (dr2 <= dist2 .and. costheta >= cosangle_crit) then
                                npairsc1c2 = npairsc1c2+1
                                flag = 1
                        end if
                end if
!
                if (Hcenters2(j,1)/=-1) then
                        costheta = 0.0
                        call cosangle(costheta,r(centers2(j)-1,:),r(Hcenters2(j,1)-1,:),r(centers1(i)-1,:),box,binv)

                        if (dr2 <= dist2 .and. costheta >= cosangle_crit) then
                                npairsc1c2 = npairsc1c2+1
                                flag = 1
                        end if
                end if
!
!                !!!! Independent "if" logics used because a center atom can accept
!                !multiple h-bonds same shell - and all of them should be calculated as
!                ! h-bonds
!                
                if (Hcenters1(i,0)/=-1) then
                        costheta = 0.0
                        call cosangle(costheta,r(centers1(i)-1,:),r(Hcenters1(i,0)-1,:),r(centers2(j)-1,:),box,binv)

                        if (dr2 <= dist2 .and. costheta >= cosangle_crit) then
                                npairsc1c2 = npairsc1c2+1
                                flag = 1
                        end if
                end if

                if (Hcenters1(i,1)/=-1) then
                        costheta = 0.0
                        call cosangle(costheta,r(centers1(i)-1,:),r(Hcenters1(i,1)-1,:),r(centers2(j)-1,:),box,binv)

                        if (dr2 <= dist2 .and. costheta >= cosangle_crit) then
                                npairsc1c2 = npairsc1c2+1
                                flag = 1
                        end if
                end if
!
                !if (flag == 1) then
                !        if ( .NOT. ANY(shellc21==i)) then
                !        shellc21(nshellc21) = i
                !        nshellc21 = nshellc21+1
                !        end if
                !end if
                !if (j == ncenters2-1) then
                !        if ( .NOT. ANY(shellc21==i)) then
                !        notshellc21(nnshellc21) = i
                !        nnshellc21 = nnshellc21+1
                !        end if
                !end if

!                if (flag == 1) then
!                        shellc21(nshellc21) = shell(i)
!                        nshellc21 = nshellc21+1
!                else if (flag == -1) then
!                       notshellc21(nnshellc21) = i
!                       nnshellc21 = nnshellc21+1
!                end if
        end if
        end do
  end do

!!else if (index1 > 0) then
!
!!  do j = 0,ncenters-1 !! j = Water or any other O in center 
!!    !print *, "Inside iteration only dist",j
!!    do i = 0,nshell-1    !! i = Water in tetrahedral shell
!!      if (centers(j)/=shell(i)) then
!!        call pbcdr(dr,r(centers(j)-1,:),r(shell(i)-1,:),box,binv)
!!        dr2 = dr(0)*dr(0)+dr(1)*dr(1)+dr(2)*dr(2)
!!        if (dr2 <= dist2) then
!!           npairs = npairs+1      
!!        end if
!!      end if
!!    end do
!!    !print *, "Iteration only dist ",j,"Ended"
!!   end do
! 
!!end if

end subroutine mednc_bond_angle
