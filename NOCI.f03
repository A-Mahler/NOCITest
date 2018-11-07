program NOCI

  use mqc_gaussian
  use mqc_algebra
  use iso_fortran_env

  implicit none
  real(kind=real64),parameter::ZERO_THRESH = 1.0E-2
  type(mqc_gaussian_unformatted_matrix_file)::temp_file
  type(mqc_wavefunction)::common_wave
  type(mqc_molecule_data)::common_mol
  type(mqc_scf_integral)::temp_int
  type(mqc_vector)::S_vec,S_inv,E_NOCI
  type(mqc_matrix)::core_ham,overlap,Hmat,Nmat,MO_I,MO_J,MIJ,rho,Gmat,MIJ_inv,eigenvecs
  type(mqc_matrix)::Umat,Vmat,tmoI,tmoJ,Smat
  type(mqc_r4tensor)::ERI
  type(mqc_scalar)::Vnn,NIJ,core_con,eri_con,half,zero,one
  integer(kind=real64)::nElectrons,nAlpha,nBeta,nBasis
  character(len=80),dimension(:),allocatable::fileList
  character(len=:),allocatable::fileName
  integer::a,i,j,k,l,m,n,io_stat_number,unitno,numFile = 0,printLevel=3
  type(mqc_matrix),dimension(:),allocatable::mo_list
  type(mqc_matrix)::S_temp,mij_temp,rho_temp,rotate,t1,t2,t3,r_Jmatrix
  type(mqc_scalar)::temp1,temp2,temp3,temp4,cos_t,sin_t
!
  call mqc_get_command_argument(1,fileName)

  open(newunit=unitno,file=fileName,status='old',iostat=io_stat_number)
  if(io_stat_number/=0) then
    call mqc_error('Error opening file',6)
  end if

  read(unit=unitno, fmt='(i2)', iostat=io_stat_number) numFile
  if(io_stat_number/=0) then
    call mqc_error('Error reading file number',6)
  end if

  allocate(fileList(numFile))

  do i = 1, numFile
    read(unit=unitno, fmt='(a)', iostat=io_stat_number) fileList(i)
    if((io_stat_number<0).and.(i<=numFile)) then
      if(allocated(fileList)) deallocate(fileList)
      if(allocated(fileName)) deallocate(fileName)
      call mqc_error('File EOF reached early',6)
    end if
  end do
!
!   If no errors thrown in opening and reading file,
!   initialize H matrix and N matrix to N_solutions x N_solutions
!   and allocate memory for mo list
!
  call Hmat%init(numFile,numFile)
  call Nmat%init(numFile,numFile)
  allocate(mo_list(numFile))

!
!   Release filename allocation and close input stream
!
  if(allocated(fileName)) deallocate(fileName)
  close(unit=unitno)
!
!   Load MO Coefficients 
!
  do i = 1, numFile
    call temp_file%getESTObj('mo coefficients',est_integral=temp_int,filename=fileList(i))
    mo_list(i) = temp_int%getBlock('full')
    if(printLevel.gt.2) then
      call mo_list(i)%print(6,'MO COEFFICIENTS')
    end if
  end do
!
!   Assign common matrices and scalars:
!   overlap, core_ham, ERI, Vnn, nelectrons,
!   nbasis, nalpha, nbeta
!
!   As all files should have the same molecule and
!   same basis set, the common elements are drawn
!   from the first file given in the file list
!
  call temp_file%load(fileList(1))
  call temp_file%getESTObj('wavefunction',common_wave)
  call temp_file%getMolData(common_mol)
  call temp_file%getArray('REGULAR 2E INTEGRALS',r4TensorOut=ERI)
  nbasis = common_wave%nbasis%ival()
  nelectrons = common_wave%nelectrons%ival()
  nalpha = common_wave%nalpha%ival()
  nbeta = common_wave%nbeta%ival()
  Vnn = mqc_get_nuclear_repulsion(molecule_info=common_mol)
  core_ham = common_wave%core_hamiltonian%getBlock('full')
  overlap = common_wave%overlap_matrix%getBlock('full')

  if(printLevel.ge.2) then
    call overlap%print(6,'overlap')
    call core_ham%print(6,'core ham')
  end if

!
!   Initialize matrices necessary for HIJ and NIJ construction
!   and scalar value
!
  call MIJ%init(nelectrons,nelectrons)
  call rho%init(nbasis*2,nbasis*2)
  half = 0.5
  zero = 0.0
  one = 1.0

! call t1%identity(4,4)
! temp1 = cos(pi/8)
! temp2 = sin(pi/8)
! temp3 = -1 * sin(pi/8)
! rotate = temp1 * t1
! call t2%identity(2,2)
! call t3%identity(2,2)
! t2 = temp2 * t2
! t3 = temp3 * t3
! call t2%print(6,'t2')
! call t3%print(6,'t3')
! call rotate%mput(t2,[3,4],[1,2])
! call rotate%mput(t3,[1,2],[3,4])

! call rotate%print(6,'rotate')

! mo_list(1) = matmul(rotate,mo_list(1))

!
!   Begin filling HIJ and NIJ matricies
!
  do i = 1, numFile
    do j = 1, numFile
      call MO_I%init(nbasis*2,nelectrons)
      call MO_J%init(nbasis*2,nelectrons)
      call rho%init(nbasis*2,nbasis*2)
      call Gmat%init(nbasis*2,nbasis*2)

!
!   Fill MO_I and MO_J with occupied orbitals; spin blocking
!   is preserved
!
      call MO_I%mput(mo_list(i)%mat([1,nBasis*2],[1,nAlpha]),[1,nBasis*2],[1,nAlpha])
      call MO_J%mput(mo_list(j)%mat([1,nBasis*2],[1,nAlpha]),[1,nBasis*2],[1,nAlpha])
      
      if(nBeta.gt.0) then
        call MO_I%mput(mo_list(i)%mat([1,nBasis*2],[nBasis+1,nBasis+nBeta]), &
          [1,nBasis*2],[nAlpha+1,nAlpha+nBeta])
        call MO_J%mput(mo_list(j)%mat([1,nBasis*2],[nBasis+1,nBasis+nBeta]), &
          [1,nBasis*2],[nAlpha+1,nAlpha+nBeta])
      end if

!
!   Build MIJ and check for NIJ threshold
!
      MIJ = matmul(dagger(MO_I),matmul(overlap,MO_J))
      NIJ = MIJ%det()

      if(printLevel.ge.2) then
        call MO_I%print(6,'MO_I')
        call MO_J%print(6,'MO_J')
        call NIJ%print(6,'NIJ')
      end if
!
!   Check for orthogonality condition
!
      if((NIJ%absval()).lt.ZERO_THRESH) then
        call MIJ%svd(EVals=S_vec,EUVecs=Umat,EVVecs=Vmat)
        
        tmoI = matmul(dagger(Umat),dagger(MO_I))
        tmoJ = matmul(MO_J,dagger(Vmat))
        
        if(printLevel.ge.3) then
          call Umat%print(6,'Umat')
          call Vmat%print(6,'Vmat')
          call S_vec%print(6,'Svec')
          call tmoI%print(6,'tmoI')
          call tmoJ%print(6,'tmoJ')
        endif

call MO_J%mput(mo_list(j)%mat([1,nBasis*2],[1,nAlpha]),[1,nBasis*2],[1,nAlpha])
142        call S_vec%diag(S_temp)
        call S_inv%init(nelectrons)

        call orth_inv(S_vec,S_inv)

        call S_inv%diag(Smat)
        
        if(printLevel.ge.3) call Smat%print(6,'Smat')

        rho = matmul(tmoJ,matmul(Smat,tmoI))
      else
        MIJ_inv = MIJ%inv()
        rho = matmul(MO_J,matmul(MIJ_inv,dagger(MO_I)))
      end if

      if(printLevel.ge.2) then
        call MIJ_inv%print(6,'MIJ_inv')
        call MIJ%print(6,'MIJ')
        call rho%print(6,'rho')
      end if
!
!   Rezero and build G matrix (ERI + Rho contraction)
!
        do k = 1, nBasis
          do l = 1, nBasis
            do m = 1, nBasis
              do n = 1, nBasis
                ! AA Block
                call Gmat%put(Gmat%at(k,l)+ERI%at(k,l,m,n)*rho%at(m,n),k,l)
                call Gmat%put(Gmat%at(k,l)+ERI%at(k,l,m,n)*rho%at(m+nBasis,n+nBasis),k,l)
                call Gmat%put(Gmat%at(k,l)-ERI%at(k,n,m,l)*rho%at(m,n),k,l)
                ! BB Block
                call Gmat%put(Gmat%at(k+nBasis,l+nBasis)+ERI%at(k,l,m,n)*rho%at(m+nBasis,n+nBasis), &
                  k+nBasis,l+nBasis)
                call Gmat%put(Gmat%at(k+nBasis,l+nBasis)+ERI%at(k,l,m,n)*rho%at(m,n), &
                  k+nBasis,l+nBasis)
                call Gmat%put(Gmat%at(k+nBasis,l+nBasis)-ERI%at(k,n,m,l)*rho%at(m+nBasis,n+nBasis), &
                  k+nBasis,l+nBasis)
                ! AB Block
                call Gmat%put(Gmat%at(k+nBasis,l)-ERI%at(k,n,m,l) * &
                  rho%at(m+nBasis,n), k+nBasis,l)
                ! BA BLock
                call Gmat%put(Gmat%at(k,l+nBasis)-ERI%at(k,n,m,l) * &
                  rho%at(m,n+nBasis), k,l+nBasis)
              end do
            end do
          end do
        end do
!
!   Perform contractions
!

        core_con = Contraction(rho,core_ham)
        eri_con = Contraction(rho,Gmat)

        if(printLevel.ge.3) then
          call Gmat%print(6,'Gmat')
          call core_con%print(6,'core contraction')
          call eri_con%print(6,'eri contraction')
        end if
!
!   Assign HIJ and NIJ 
!
        if((NIJ%absval()).lt.ZERO_THRESH) then
          call Hmat%put((core_con + (half*eri_con)),I,J)
        else
          call Hmat%put(NIJ*((core_con + (half*eri_con))),I,J)
        endif
        call Nmat%put(NIJ,I,J)
    end do
  end do

  call E_NOCI%init(numFile)

  call mqc_matrix_generalized_eigensystem(Hmat,Nmat,E_NOCI,eigenvecs)

  do i = 1, numFile
    call E_NOCI%put(E_NOCI%at(i) + Vnn, i)
  end do
  call E_NOCI%print(6,'NOCI Energies')
  if(printLevel.ge.1) then
    call eigenvecs%print(6,'Eigenvectors')
    call Hmat%print(6,'HIJ')
    call Nmat%print(6,'NIJ')
    call Vnn%print(6,'Vnn')
  end if

!
!   Release allocatables
!
  if(allocated(fileList)) deallocate(fileList)
  if(allocated(mo_list)) deallocate(mo_list)

end program NOCI

subroutine orth_inv(Svals,Sinv)

  use iso_fortran_env
  use mqc_algebra

  implicit none
  real(kind=real64),parameter::Z_THRESH = 1.0E-3
  type(mqc_vector)::Svals,Sinv
  type(mqc_scalar)::one,zero,pNIJ,hold,s_temp
  real(kind=real64)::temp
  integer(kind=int64)::i,nvals

  one = 1.0
  zero = 0.0
  nvals = Svals%size()

  do i = 1, nvals
    hold = Svals%at(i)
    temp = hold%rval()

    if(temp.le.Z_THRESH) then
      s_temp = 1.0
    else if((1-temp).le.Z_THRESH) then
      s_temp = 0.0
    else
      s_temp = 0.0
    end if

    call Sinv%put(s_temp,i)
  end do
  
  call Svals%print(6,'Svals')
  call Sinv%print(6,'Sinv')

end subroutine orth_inv
