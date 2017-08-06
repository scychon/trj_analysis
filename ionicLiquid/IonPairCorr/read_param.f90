MODULE routines

  !****************************************
  ! this module contains miscellaneous subroutines that are used in the program
  !****************************************

CONTAINS


  !***********************************************************
  ! this subroutine gets the number of atoms for each molecules, as well as the number
  ! of datapoints in the file.  The number of data points is recognized as the number of lines that
  ! start with an integer, divided by two




  !***********************************************
  ! this subroutine reads in force field parameters, as well as any 
  ! hard constraints
  !
  !**********************************************
  subroutine getparameters(param_file)
    use variables
    CHARACTER(*),INTENT(in)::param_file
    INTEGER::i,j,k,atoms1,atoms2,inputstatus,atom,site1,site2,ind=0,flag
    CHARACTER(50)::line
    character(5)::junk1
    character(len(atomtypemon1(1)))::test1
    REAL*8::temp,junk,exch,e1,e2,temp_pen(4)
    real*8,dimension(:),allocatable::exp_molecule

    atoms1=size(chargeatoms1)
    atoms2=SIZE(chargeatoms2)

    OPEN(unit=7,file=param_file,status='old')


    !********************************* hard constraints******************************
    if ( hard_constraints .eq. "yes" ) then
       read(7,*) atom 
       allocate(hard_cons_type(atom),hard_cons_param(atom,4)) 
       do i=1,atom
          read(7,*) hard_cons_type(i),(hard_cons_param(i,j),j=1,4)
       enddo

       ! fill in hard constraints for both molecules
       do i=1,atoms2
          do j=1,size(hard_cons_type)
             if ( atomtypemon2(i) .eq. hard_cons_type(j) ) then
                Exchatom2(i) = hard_cons_param(j,1)
                Elecatom2(i) = hard_cons_param(j,2)
                Inducatom2(i) = hard_cons_param(j,3)
                dhf2(i) = hard_cons_param(j,4)
             endif
          enddo
       enddo
       do i=1,atoms1
          do j=1,size(hard_cons_type)
             if ( atomtypemon1(i) .eq. hard_cons_type(j) ) then
                Exchatom1(i) = hard_cons_param(j,1)
                Elecatom1(i) = hard_cons_param(j,2)
                Inducatom1(i) = hard_cons_param(j,3)
                dhf1(i) = hard_cons_param(j,4)
             endif
          enddo
       enddo

    endif


!!!!!!!!!!!!!!!!!!!read exchange exponents, create cross terms
    allocate(exp_molecule(atoms1))
    READ(7,'(A)') line
    DO i=1,atoms1
       READ(7,*) junk1,exp_molecule(i)
    ENDDO
    write(*,*) ""
    write(*,*) "**************************************************************"
    write(*,*) "explicitly creating cross-term exponents using combination rule"
    write(*,*) "exp_tot = (exp1 + exp2 ) * exp1*exp2 / (exp1^2 + exp2^2 )"
    write(*,*) "**************************************************************"  
    write(*,*) "" 
    do i=1,atoms1
       e1 = exp_molecule(i)
       do j=1,atoms2
          e2=exp_molecule(j)
          exponents(i,j) = (e1+e2) * (e1*e2)/(e1**2+e2**2)
       enddo
    enddo


!!!!!!!!!!!!!!!!!!!!! read dispersion coefficients, C6-C12, create cross terms
    READ(7,'(A)') line
    DO i=1,atoms1
       READ(7,*) junk1, (Cn_cross(i,i,j),j=1,4)
    ENDDO
    write(*,*) ""
    write(*,*) "**************************************************************"
    write(*,*) "explicitly creating cross-term dispersion coeffs using combination rule"
    write(*,*) "Cn(AB) = sqrt(Cn(A) * Cn(B))"
    write(*,*) "**************************************************************"  
    write(*,*) "" 
    do i=1,atoms1
       do j=1,atoms2
          do k=1,4
             Cn_cross(i,j,k) = dsqrt( Cn_cross(i,i,k) * Cn_cross(j,j,k) )
          enddo
       enddo
    enddo




!!!!!!!!!!!!!!!!!!!!!!read charges
    READ(7,'(A)') line
    DO i=1,atoms1
       READ(7,*) atom, chargeatoms1(i)
    ENDDO
    chargeatoms2=chargeatoms1


!!!!!!!!!!!!!!!!!!!!!!! read drude charges and springcon
    READ(7,'(A)') line
    DO i=1,atoms1
       READ(7,*) atom, shellcharge1(i)
    ENDDO
    READ(7,'(A)') line
    READ(7,*) springcon1      

    close(7)



  END SUBROUTINE getparameters



!**********************************************************************

       

END MODULE routines
