program NOCI

  use mqc_gaussian
  use mqc_algebra
  use iso_fortran_env

  implicit none
  real(kind=real64),parameter::NIJ_THRESH = 1.0E-10
  type(mqc_gaussian_unformatted_matrix_file)::temp_file
  type(mqc_wavefunction)::common_wave
  type(mqc_molecule_data)::common_mol
  type(mqc_scf_integral)::temp_int
  type(mqc_vector)::E_NOCI
  type(mqc_matrix)::core_ham,overlap,Hmat,Nmat,MO_I,MO_J,MIJ,rho,Gmat
  type(mqc_r4tensor)::ERI
  type(mqc_scalar)::Vnn,NIJ,core_con,eri_con,half
  integer(kind=real64)::nElectrons,nAlpha,nBeta,nBasis
  character(len=80),dimension(:),allocatable::fileList
  character(len=:),allocatable::fileName
  integer::i,j,k,l,m,n,argnum,io_stat_number,unitno,numFile = 0
  type(mqc_matrix),dimension(:),allocatable::mo_list

  call mqc_get_command_argument(1,fileName)
  
  open(newunit=unitno,file=fileName,status='old',iostat=io_stat_number)
  if(io_stat_number/=0) then
    call mqc_error('Error opening file',6)
  end if

  read(unit=unitno, fmt='(i1)', iostat=io_stat_number) numFile
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
  Vnn = mqc_get_nuclear_repulsion(6,common_mol)
  core_ham = common_wave%core_hamiltonian%getBlock('full')
  overlap = common_wave%overlap_matrix%getBlock('full')
!
!   Initialize matrices necessary for HIJ and NIJ construction
!   and scalar value
!
  call MIJ%init(nelectrons,nelectrons,cmplx(0.0,0.0))
  call rho%init(nbasis*2,nbasis*2,cmplx(0.0,0.0))
  half = 0.5
!
!   Begin filling HIJ and NIJ matricies
!
  do i = 1, numFile
    do j = 1, numFile
      call MO_I%init(nbasis*2,nelectrons,cmplx(0.0,0.0))
      call MO_J%init(nbasis*2,nelectrons,cmplx(0.0,0.0))
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

      if((i>2).or.(j>2)) then
        call MO_I%print(6, 'MO_I')
        call MO_J%print(6, 'MO_J')
      end if
!
!   Build MIJ and check for NIJ threshold
!
      MIJ = matmul(matmul(dagger(MO_I),overlap),MO_J)
      if(i.eq.j) then
        call MIJ%print(6,'MIJ')
      endif
      NIJ = MIJ%det()

      if((NIJ%rval()).gt.NIJ_THRESH) then
        rho = matmul(matmul(MO_J,MIJ%inv()),dagger(MO_I))
!
!   Rezero and build G matrix (ERI + Rho contraction)
!
        call Gmat%init(nBasis*2,nBasis*2,cmplx(0.0,0.0))
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
                call Gmat%put(Gmat%at(k+nBasis,l)-ERI%at(k,n,m,l)*rho%at(m+nBasis,n), &
                  k+nBasis,l)
                ! BA BLock
                call Gmat%put(Gmat%at(k,l+nBasis)-ERI%at(k,n,m,l)*rho%at(m,n+nBasis), &
                  k,l+nBasis)
              end do
            end do
          end do
        end do
!
!   Perform contractions
!
        core_con = Contraction(rho,core_ham)
        eri_con = Contraction(rho,Gmat)
!
!   Assign HIJ and NIJ 
!   
        call Hmat%put(NIJ*((core_con + (half*eri_con)) + Vnn),I,J)
        call Nmat%put(NIJ,I,J)
      end if
    end do
  end do

  call Hmat%print(6,'HIJ')
  call Nmat%print(6,'NIJ')

  call E_NOCI%init(numFile)
  call mqc_matrix_generalized_eigensystem(Hmat,Nmat,E_NOCI)

  call E_NOCI%print(6,'NOCI Energies')
!
!   Release allocatables
!
  if(allocated(fileList)) deallocate(fileList)
  if(allocated(mo_list)) deallocate(mo_list) 

end program NOCI
