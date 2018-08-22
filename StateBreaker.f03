program StateBreaker

  use mqc_gaussian
  use mqc_algebra
  use iso_fortran_env

  implicit none
  real(kind=real64),parameter::offset1 = 0.5 * (pi/180)
  real(kind=real64),parameter::offset2 = 0.5 * (pi/180)
  integer::i,j,k,l,unitno,numFile,io_stat_number,iOut=6
  character(len=80),dimension(2)::fileList
  character(len=:),allocatable::fileName
  type(mqc_gaussian_unformatted_matrix_file)::temp_file
  type(mqc_scf_integral)::temp_int
  type(mqc_r4tensor)::ERI
  type(mqc_matrix)::core_ham,overlap
  type(mqc_matrix),dimension(2)::mo_list,occ_mo
  type(mqc_matrix),dimension(4)::psi
  type(mqc_vector)::E_NOCI
  type(mqc_scalar)::theta1,theta2,cos_theta1,cos_theta2,sin_theta1,sin_theta2,e1,e2,Vnn
  type(mqc_wavefunction)::common_wave
  type(mqc_molecule_data)::common_mol

  call mqc_get_command_argument(1,fileName)

  open(newunit=unitno, file=fileName, status='old', iostat=io_stat_number)
  if(io_stat_number/=0) call mqc_error('Error opening file',iOut)

  do i = 1, 2
    read(unit=unitno, fmt='(a)', iostat=io_stat_number) fileList(i)
    if((io_stat_number<0).and.(i<=numFile)) then
      if (allocated(fileName)) deallocate(fileName)
      call mqc_error('EOF reached early',iOut)
    end if
  end do

  if(allocated(fileName)) deallocate(fileName)
  close(unit=unitno)

! Load MO coefficients into mo_list array
  do i = 1, 2
    call temp_file%getESTObj('mo coefficients',est_integral=temp_int,filename=fileList(i))
    mo_list(i) = temp_int%getBlock('full')
  end do

! Extract common matrix elements

  call temp_file%load(fileList(1))
  call temp_file%getArray('REGULAR 2E INTEGRALS',r4TensorOut=ERI)
  call temp_file%getESTObj('wavefunction',common_wave)
  call temp_file%getMolData(common_mol)
  core_ham = common_wave%core_hamiltonian%getBlock('full')
  overlap = common_wave%overlap_matrix%getBlock('full')
  Vnn = mqc_get_nuclear_repulsion(molecule_info=common_mol)

! Initialize occ_mo matricies
  do i = 1, 2
    call occ_mo(i)%init(4,2)
  end do

! Extract the occupied orbitals from MO coefficients
  do i = 1, 2
    call occ_mo(i)%mput(mo_list(i)%mat([1,4],[1,1]),[1,4],[1,1])
    call occ_mo(i)%mput(mo_list(i)%mat([1,4],[3,3]),[1,4],[2,2])
  end do

! Begin primary algorithm
!
!
  do i = 1, 4
    call psi(i)%init(4,1)
  end do

  do i = 1, 90
    do j = 1, 90
      theta1 = 0.0 + (i*offset1)
      theta2 = 0.0 + (i*offset2)
      cos_theta1 = cos(theta1%rval())
      sin_theta1 = sin(theta1%rval())
      cos_theta2 = cos(theta2%rval())
      sin_theta2 = sin(theta2%rval())


! Rotate ground state orbitals
      call psi(1)%mput((cos_theta1*mo_list(1)%mat([1,4],[1,1]))+((sin_theta1)*mo_list(1)%mat([1,4],[2,2])),[1,4],[1,1])
      call psi(2)%mput((cos_theta1*mo_list(1)%mat([1,4],[3,3]))-((sin_theta1)*mo_list(1)%mat([1,4],[4,4])),[1,4],[1,1])
      call psi(3)%mput((cos_theta2*mo_list(2)%mat([1,4],[1,1]))+((sin_theta2)*mo_list(2)%mat([1,4],[2,2])),[1,4],[1,1])
      call psi(4)%mput((cos_theta2*mo_list(2)%mat([1,4],[3,3]))-((sin_theta2)*mo_list(2)%mat([1,4],[4,4])),[1,4],[1,1])
!
      call occ_mo(1)%mput(psi(1)%mat([1,4],[1,1]),[1,4],[1,1])
      call occ_mo(1)%mput(psi(2)%mat([1,4],[1,1]),[1,4],[2,2])
      call occ_mo(2)%mput(psi(3)%mat([1,4],[1,1]),[1,4],[1,1])
      call occ_mo(2)%mput(psi(4)%mat([1,4],[1,1]),[1,4],[2,2])
      E_NOCI = NOCI(occ_mo,2,2,overlap,core_ham,ERI)

      e1 = E_NOCI%at(1)
      e2 = E_NOCI%at(2)

      e1 = e1 + Vnn
      e2 = e2 + Vnn

      print *, theta1%rval(), theta2%rval(), e1%rval(), e2%rval()
    end do
  end do
 
contains

function NOCI(states,nStates,nBasis,overlap,core_ham,ERI)

  implicit none
  integer::nStates,nBasis,i,j,k,l,m,n
  type(mqc_scalar)::NIJ,core_con,eri_con,half
  type(mqc_matrix),dimension(:)::states
  type(mqc_matrix)::overlap,core_ham,MIJ,rho,Gmat,Hmat,Nmat,Dmat
  type(mqc_r4tensor)::ERI
  type(mqc_vector)::NOCI

! Initialize Matricies
  half = 0.5
  call NOCI%init(nStates,1.0)
  call Hmat%init(nStates,nStates)
  call Nmat%init(nStates,nStates)


! Build HIJ and NIJ matricies
  do i = 1, nStates
    do j = 1, nStates

!     Zero or Re-zero G Matrix

      call Gmat%init(nBasis*2,nBasis*2)

!     Build MIJ and rho matricies
      MIJ = matmul(matmul(dagger(states(i)),overlap),states(j))
      NIJ = MIJ%det()
      rho = matmul(matmul(states(j),MIJ%inv()),dagger(states(i)))
      
      if(NIJ%absval().gt.(1E-10)) then

!     Perform G-matrix contraction
        do k = 1, nBasis
          do l = 1, nBasis
            do m = 1, nBasis
              do n = 1, nBasis
                ! AA Block
                call Gmat%put(Gmat%at(k,l)+ERI%at(k,l,m,n)*rho%at(m,n),k,l)
                call Gmat%put(Gmat%at(k,l)+ERI%at(k,l,m,n)*rho%at(m+nBasis,n+nBasis),k,l)
                call Gmat%put(Gmat%at(k,l)-ERI%at(k,n,m,l)*rho%at(m,n),k,l)
                !BB Block
                call Gmat%put(Gmat%at(k+nBasis,l+nBasis)+ERI%at(k,l,m,n)*rho%at(m+nBasis,n+nBasis), &
                  k+nBasis,l+nBasis)
                call Gmat%put(Gmat%at(k+nBasis,l+nBasis)+ERI%at(k,l,m,n)*rho%at(m,n), &
                  k+nBasis,l+nBasis)
                call Gmat%put(Gmat%at(k+nBasis,l+nBasis)-ERI%at(k,n,m,l)*rho%at(m+nBasis,n+nBasis), &
                  k+nBasis,l+nBasis)
                !AB Block
                call Gmat%put(Gmat%at(k+nBasis,l)-ERI%at(k,n,m,l)*rho%at(m+nBasis,n), &
                  k+nBasis,l)
                !BA Block
                call Gmat%put(Gmat%at(k,l+nBasis)-ERI%at(k,n,m,l)*rho%at(m,n+nBasis), &
                  k,l+nBasis)
              end do
            end do
          end do
        end do

!   Perform final contractions

        core_con = Contraction(rho,core_ham)
        eri_con = Contraction(rho,Gmat)

!   Assign elements of HIJ and NIJ
        call Hmat%put(NIJ*((core_con + (half*eri_con))),i,j)
        call Nmat%put(NIJ,i,j)
      end if
    end do
  end do

! Solve Eigensystem

  call mqc_matrix_generalized_eigensystem(Hmat,Nmat,NOCI,Dmat)

end function NOCI

end program StateBreaker
