program lattice_2dim

use mtmod

implicit none

integer,parameter :: Naa=100    !!アミノ酸の個数
integer,parameter :: Niter = 1000000000
integer,parameter :: Nrep = 100
integer,parameter :: istep_rep = 10
integer,parameter :: istep_save = 100
integer,parameter :: istep_save_xy = 10000

real,parameter :: penalty = 5.0

integer :: iseed = 333  
real,parameter :: T_high = 2.5, T_low=0.025

real, parameter :: prob_single = 0.2
real, parameter :: prob_rigid  = 0.5

integer::i,j,k,istep, irep
integer::i_target

!!アミノ酸が疎水性なら1、親水性なら0として並べる
!integer,dimension(Naa)::aminostate=(/0,0,1,1,1,0,1,1,1,1,1,1,1,1,0,0,0,1,1,1, &
!                                     1,1,1,1,1,1,1,0,1,0,0,0,1,1,1,1,1,1,1,1, &
!                                     1,1,1,1,0,0,0,0,1,1,1,1,1,1,0,1,1,0,1,0/)
integer,parameter,dimension(Naa)::aminostate= &
(/0,0,0,0,0,0,1,0,1,1,0,0,0,0,0,1,1,1,0,1,    &
  1,1,1,1,0,1,1,0,0,0,0,1,1,0,0,1,1,0,1,1,    &
  1,1,1,0,1,1,1,1,1,1,1,1,1,1,0,1,1,0,1,1,    &
  1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,    &
  1,1,1,0,0,1,0,1,1,1,0,0,0,0,0,0,1,0,1,1/)

real::E(Nrep),Emin,Etry,Emin_self
!real::T(Nrep)
integer :: rep2lab(Nrep)
integer :: lab2rep(Nrep)
real :: lab2temp(Nrep)

integer :: xy(2,Naa)
integer :: xy2(2,Naa)
integer :: vec(Naa-1)
integer :: loc(Naa-2,Nrep), loc_try(Naa-2)
integer :: loc_min(Naa-2), loc_min_self(Naa-2)
!logical :: flg_allow
logical :: flg_self
logical :: flg_odd = .true.
integer :: i_start
real    :: delta
integer :: rep_i, rep_j, lab_i, lab_j
real    :: E_i, E_j
real    :: temp_i, temp_j

call sgrnd(iseed)

lab2temp(1) = T_low
do i = 2, Nrep-1
   lab2temp(i) = T_low + (T_high-T_low) / (Nrep-1) * (i-1)
enddo
lab2temp(Nrep) = T_high

do i = 1, Nrep
   rep2lab(i) = i
   lab2rep(i) = i
enddo

open(10,file='ts.data')
open(11,file='rep.data')
open(25,file='min_xy.data')
open(26,file='min_state.data')
open(27,file='self_xy.data')
open(28,file='self_state.data')

!始めの構造のエネルギーの決定-------------------------------------
loc(:,:) = 0
E(:) = 0
loc_min(:) = 0
Emin = 0
loc_min_self(:) = 0
Emin_self = 0


do istep=1,Niter

!   T = Tini-real(istep)*((Tini-Tfin)/real(Niter))
   do irep = 1, Nrep

      loc_try(:) = loc(:,irep)

      !flg_allow = .false.
      !do while(.not. flg_allow)
      do 
         if (grnd() < prob_single) then
            i_target = int(Naa * grnd()) + 1

            if (i_target == 1 .OR. i_target == Naa) then
               call move_terminal(Naa,i_target,loc_try)

            else
               ! rigid body
               if (grnd() < prob_rigid) then
                  call move_rigid(Naa,i_target,loc_try)

               ! corner flip
               else
                  if (loc_try(i_target-1) == 0) then  ! internal line
                     call move_rigid(Naa,i_target,loc_try)
                     !cycle
                  else                                ! not line
                     ! crankshaft
                     if ((i_target > 2     .AND. (loc_try(i_target-2) == loc_try(i_target-1))) .OR.&
                         (i_target < Naa-1 .AND. (loc_try(i_target)   == loc_try(i_target-1)))) then
                        cycle
                     ! internal corner
                     else
                        call move_corner(Naa,i_target,loc_try)
                     endif
                  endif
               endif
            endif

         else
            i_target = int((Naa-3) * grnd()) + 1
            if (loc_try(i_target-1) /= loc_try(i_target)) then
               cycle
            endif
            if (i_target > 2) then
               if (loc_try(i_target-2) == 0) then
                  cycle
               endif
            endif
            if (i_target < Naa-2) then
               if (loc_try(i_target+1) == 0) then
                  cycle
               endif
            endif

            call move_crankshaft(Naa,i_target,loc_try)
   
         endif

!         call check_occupancy(Naa,loc_try,flg_allow)
         !flg_allow = .true.
         exit
      enddo

      !! Energy
      call loc2xy(Naa,loc_try,xy)
      call calc_energy(Naa,aminostate,xy,Etry,flg_self)

      !! Judgement
      if (Etry < E(irep)) then ! accept 
         E(irep) = Etry
         loc(:,irep) = loc_try(:)
      else
         if(grnd() < exp((E(irep)-Etry)/rep2temp(irep))) then !accept
            E(irep) = Etry
            loc(:,irep) = loc_try(:)
         end if
      end if

      if (Etry < Emin) then
         Emin = Etry
         loc_min(:) = loc_try(:)
      endif
      if (flg_self) then
         if (Etry < Emin_self) then
            Emin_self = Etry
            loc_min_self(:) = loc_try(:)
         endif
      endif
   enddo
   
   !! output
   if (mod(istep,istep_save) == 0) then
      write(10,'(i12,a,f5.1,a,f5.1,a)',advance="no") istep,' ', Emin,' ',Emin_self,' '
      write(11,'(i12,a,f5.1,a,f5.1,a)',advance="no") istep,' ', Emin,' ',Emin_self,' '
      do i = 1, Nrep-1
         write(10,'(f5.1,a)',advance="no") E(i) ,' '
         write(11,'(f7.4,a)',advance="no") rep2temp(i) , ' '
      enddo
      write(10,'(f5.1)') E(Nrep)
      write(11,'(f7.4)') rep2temp(Nrep)
      if (mod(istep,istep_save_xy) == 0) then
         call loc2xy(Naa,loc_min,xy)
         call loc2xy(Naa,loc_min_self,xy2)
         write(25,*) ''
         write(25,*) '#',istep
         write(26,*) ''
         write(26,*) '#',istep
         write(27,*) ''
         write(27,*) '#',istep
         write(28,*) ''
         write(28,*) '#',istep
         do i=1,Naa
            write(25,*) xy(1,i),xy(2,i)
            write(27,*) xy2(1,i),xy2(2,i)
            if(aminostate(i)==1)then
               write(26,*) xy(1,i),xy(2,i)
               write(28,*) xy2(1,i),xy2(2,i)
            end if
         end do
      endif
   endif

   !! Replica exchange
   if (mod(istep,istep_rep) == 0) then
      if (flg_odd) then
         i_start = 1
      else
         i_start = 2
      endif

      do lab_i = i_start, (Nrep-1), 2
         lab_j = lab_i + 1
         rep_i = lab2rep(lab_i)
         rep_j = lab2rep(lab_j)
         temp_i = lab2temp(lab_i)
         temp_j = lab2temp(lab_j)
         E_i = E(rep_i)
         E_j = E(rep_j)

         delta = (1/temp_j - 1/temp_i) * (E_i - E_j)
         if (delta > 0) then
            if (grnd() > exp(-delta)) then
               cycle  ! reject
            endif
         endif

         !exchange (accept)
         lab2rep(lab_i) = rep_j
         lab2rep(lab_j) = rep_i
         rep2lab(rep_i) = lab_j
         rep2lab(rep_j) = lab_i
      enddo
            
      flg_odd = .not. flg_odd
   endif

end do

write(*,*)Emin

!!-------------------------------------------------------------------------------
! final output
call loc2xy(Naa,loc_min,xy)
call loc2xy(Naa,loc_min_self,xy2)
write(25,*) ''
write(25,*) '#'
write(26,*) ''
write(26,*) '#'
write(27,*) ''
write(27,*) '#'
write(28,*) ''
write(28,*) '#'
do i=1,Naa
   write(25,*) xy(1,i),xy(2,i)
   write(27,*) xy2(1,i),xy2(2,i)
   if(aminostate(i)==1)then
      write(26,*) xy(1,i),xy(2,i)
      write(28,*) xy2(1,i),xy2(2,i)
   end if
end do

close(25)
close(26)
close(27)
close(28)
close(10)
close(11)

stop
!#########################################################################################
!#########################################################################################

CONTAINS

!!------------------------------------------------------------------------------------------
SUBROUTINE  calc_energy(Naa,aminostate,xy,E,flg)
implicit none
integer,intent(in)::Naa
integer,dimension(Naa),intent(in)::aminostate
integer, intent(in) :: xy(2,Naa)
real,intent(out)::E
logical,intent(out) :: flg  ! self-avoiding or not
integer :: occupy(-Naa:Naa, -Naa:Naa)
integer :: i,j
   occupy(:,:) = 0
   E = 0.0
   flg = .true.
   do i=1,Naa
      !if (xy(1,i) == xy(1,i+1) .AND. xy(2,i) == xy(2,i+1)) then
      !   E=E+penalty
      !endif
      occupy(xy(1,i),xy(2,i)) = occupy(xy(1,i),xy(2,i)) + 1
      do j=i+2,Naa
         if ((abs(xy(1,i)-xy(1,j)) == 1 .AND. abs(xy(2,i)-xy(2,j)) == 0) .OR. &
             (abs(xy(2,i)-xy(2,j)) == 1 .AND. abs(xy(1,i)-xy(1,j)) == 0)) then
            if (aminostate(i) == 1 .AND. aminostate(j) == 1) then 
               E = E - 1.0
            endif
         end if
      end do
   end do 

   do i=-Naa,Naa
      do j=-Naa,Naa
         if (occupy(i,j) > 1) then
            E = E + penalty * ((occupy(i,j)-1) ** 2)
            flg = .false.
         endif
      enddo
   enddo
   return
END SUBROUTINE calc_energy


!========================================================================
!========================================================================


subroutine loc2vec(Naa,loc,vec)
   implicit none
   integer, intent(in) :: Naa
   integer, intent(in) :: loc(Naa-2)
   integer, intent(out) :: vec(Naa-1)
   integer :: i

   vec(1) = 3
   do i = 1, Naa-2
      if (loc(i) == -1) then
         vec(i+1) = vec(i) + 3
         if (vec(i+1) == 12) then
            vec(i+1) = 0
         endif
      elseif (loc(i) == 0) then
         vec(i+1) = vec(i)
      elseif (loc(i) == +1) then
         vec(i+1) = vec(i) - 3
         if (vec(i+1) == -3) then
            vec(i+1) = 9
         endif
      else
         write(*,*) 'Error'
         stop
      endif
   enddo
endsubroutine loc2vec

subroutine vec2xy(Naa,vec,xy)
   implicit none
   integer, intent(in) :: Naa
   integer, intent(in) :: vec(Naa-1)
   integer, intent(out) :: xy(2,Naa)
   integer :: i
   
   xy(:,1) = (/0,0/)
   do i = 1, Naa-1
      if (vec(i) == 3) then
         xy(1,i+1) = xy(1,i) + 1
         xy(2,i+1) = xy(2,i)
      elseif (vec(i) == 0) then
         xy(1,i+1) = xy(1,i) 
         xy(2,i+1) = xy(2,i) + 1
      elseif (vec(i) == 9) then
         xy(1,i+1) = xy(1,i) - 1
         xy(2,i+1) = xy(2,i) 
      elseif (vec(i) == 6) then
         xy(1,i+1) = xy(1,i) 
         xy(2,i+1) = xy(2,i) - 1
      else
         write(*,*) 'Error'
         stop
      endif
   enddo
endsubroutine vec2xy

subroutine loc2xy(Naa,loc,xy)
   implicit none
   integer, intent(in) :: Naa
   integer, intent(in) :: loc(Naa-2)
   integer, intent(out) :: xy(2,Naa)
   integer :: vec, vec_pre

   xy(1,1) = 0
   xy(2,1) = 0
   xy(1,2) = 1
   xy(2,2) = 0
   vec_pre = 3
   do i = 2, Naa-1
      if (loc(i-1) == -1) then
         vec = vec_pre + 3
         if (vec == 12) then
            vec = 0
         endif
      elseif (loc(i-1) == 0) then
         vec = vec_pre
      elseif (loc(i-1) == +1) then
         vec = vec_pre - 3
         if (vec == -3) then
            vec = 9
         endif
      else
         write(*,*) 'Error 1'
         stop
      endif

      if (vec == 3) then
         xy(1,i+1) = xy(1,i) + 1
         xy(2,i+1) = xy(2,i)
      elseif (vec == 0) then
         xy(1,i+1) = xy(1,i) 
         xy(2,i+1) = xy(2,i) + 1
      elseif (vec == 9) then
         xy(1,i+1) = xy(1,i) - 1
         xy(2,i+1) = xy(2,i) 
      elseif (vec == 6) then
         xy(1,i+1) = xy(1,i) 
         xy(2,i+1) = xy(2,i) - 1
      else
         write(*,*) 'Error 2', i, vec
         stop
      endif
      vec_pre = vec
   enddo
endsubroutine loc2xy

subroutine move_terminal(Naa,i_target,loc_try)
   implicit none
   integer, intent(in) :: Naa, i_target
   integer, intent(inout) :: loc_try(Naa-2)
   integer :: loc_target

   if (i_target == 1) then
      loc_target = 1
   elseif (i_target == Naa) then
      loc_target = Naa-2
   else
      write(*,*) 'Error'
      stop
   endif

   if (grnd() < 0.5) then
      loc_try(loc_target) = loc_try(loc_target) + 1
   else
      loc_try(loc_target) = loc_try(loc_target) - 1
   endif

   if (loc_try(loc_target) == 2) then
      loc_try(loc_target) = -1
   elseif (loc_try(loc_target) == -2) then
      loc_try(loc_target) = +1
   endif
endsubroutine move_terminal

subroutine move_rigid(Naa,i_target,loc_try)
   implicit none
   ! 2 <= i_target <= Naa-1
   integer, intent(in) :: Naa, i_target
   integer, intent(inout) :: loc_try(Naa-2)
   integer :: loc_target 

   loc_target = i_target - 1

   if (grnd() < 0.5) then
      loc_try(loc_target) = loc_try(loc_target) + 1
   else
      loc_try(loc_target) = loc_try(loc_target) - 1
   endif

   if (loc_try(loc_target) == -2) then
      loc_try(loc_target) =  1
   else if (loc_try(loc_target) == 2) then
      loc_try(loc_target) = -1
   endif
endsubroutine move_rigid

subroutine move_corner(Naa,i_target,loc_try)
   implicit none
   ! 2 <= i_target <= Naa-1
   integer, intent(in) :: Naa, i_target
   integer, intent(inout) :: loc_try(Naa-2)
   integer :: loc_target 

   loc_target = i_target - 1

   if (loc_target-1 > 0) then
      loc_try(loc_target-1) = loc_try(loc_target-1) + loc_try(loc_target)
   endif

   if (loc_target+1 < Naa-1) then
      loc_try(loc_target+1) = loc_try(loc_target+1) + loc_try(loc_target)
   endif

   loc_try(loc_target) = loc_try(loc_target) * (-1)
endsubroutine move_corner

subroutine move_crankshaft(Naa,i_target,loc_try)
   implicit none
   integer, intent(in) :: Naa, i_target
   integer, intent(inout) :: loc_try(Naa-2)
   integer :: loc_target

   loc_target = i_target - 1
   
   if (loc_target > 1) then
      loc_try(loc_target-1) = loc_try(loc_target-1) * (-1)
   endif
   loc_try(loc_target)   = loc_try(loc_target)   * (-1)
   loc_try(loc_target+1) = loc_try(loc_target+1) * (-1)
   if (loc_target < Naa-3) then
      loc_try(loc_target+2) = loc_try(loc_target+2) * (-1)
   endif
endsubroutine move_crankshaft

subroutine check_occupancy(Naa,loc,flg)
   implicit none
   integer, intent(in) :: Naa, loc(Naa-2)
   logical, intent(out) :: flg

   integer :: i
   integer :: vec(Naa-1), xy(2,Naa)
   
   call loc2vec(Naa,loc,vec)
   call vec2xy(Naa,vec,xy)

   do i = 1, Naa
      do j = i+1, Naa
         if (xy(1,i) == xy(1,j) .AND. xy(2,i) == xy(2,j)) then
            flg = .false.
            return
         endif
      enddo
   enddo

   flg = .true.
   return
endsubroutine check_occupancy

real function rep2temp(ii)
   implicit none
   integer, intent(in) :: ii  
   rep2temp = lab2temp(rep2lab(ii))
endfunction rep2temp

END PROGRAM lattice_2dim
