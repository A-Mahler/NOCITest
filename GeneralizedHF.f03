program GeneralizedHF

  use mqc_gaussian
  use mqc_algebra
  use iso_fortran_env

  implicit none
  type(mqc_gaussian_unformatted_matrix_file)::matFile
  type(mqc_wavefunction)::wavefunc
  type(mqc_molecule_data)::molInfo
  type(mqc_matrix)::MO_I,MO_J,rho,MIJ,MIJ_inv,MOcoeff,overlap,core_ham,Gmat,p_Fock,Idem,ps
  type(mqc_r4tensor)::ERI
  type(mqc_scalar)::half,core_con,eri_con,Vnn,HIJ,ps_trace,tAA,tBB,tAB,tBA
  character(len=:),allocatable::fileName
  integer::i,j,k,l,nelectrons,nalpha,nbeta,nbasis

  call mqc_get_command_argument(1,fileName)
  call matFile%load(fileName)
  call matFile%getESTObj('wavefunction',wavefunc)
  call matFile%getArray('REGULAR 2E INTEGRALS',r4TensorOut=ERI)
  call matFile%getMolData(molInfo)
  call overlap%initialize(nBasis*2,nBasis*2,(0.0,0.0))
  nelectrons = wavefunc%nElectrons%ival()
  nalpha = wavefunc%nAlpha%ival()
  nbeta = wavefunc%nBeta%ival()
  nbasis = wavefunc%nbasis%ival()
  overlap = wavefunc%overlap_matrix%getBlock('full')
  MOcoeff = wavefunc%mo_coefficients%getBlock('full')
  core_ham = wavefunc%core_hamiltonian%getBlock('full')
  Vnn = mqc_get_nuclear_repulsion(6,molInfo)
  half = 0.5

  call overlap%print(6,'Overlap')
  call MOcoeff%print(6,'MO Coefficients')

!  if(IAND(NELECTRONS,1) .eq. 0) then
!    nalpha = (nelectrons/2)
!  else
!    nalpha = (nelectrons/2) + 1
!  end if
!  nbeta = nelectrons - nalpha
!  initialize occupied MO coefficient matrices (2*nBasis,nAlpha+nBeta)
  
  call core_ham%print(6,'Core Hamiltonian')
  call MO_I%initialize(nBasis*2,nElectrons,(0.0,0.0))
  call MO_J%initialize(nBasis*2,nElectrons,(0.0,0.0))

!  propogate MO_I matrices for alpha/beta electrons
  print *, 'alpha : beta: ', nalpha, ' : ', nbeta
  call MO_I%mput(MOcoeff%mat([1,nBasis*2],[1,nalpha]),[1,nBasis*2],[1,nalpha])
!  call MO_I%mput(MOcoeff%mat([1,nBasis*2],[nBasis+1,nBasis+nbeta]),[1,nBasis*2],[nalpha+1,nalpha+nbeta])
  if(nbeta.gt.0) then
    call MO_I%mput(MOcoeff%mat([1,nBasis*2],[nBasis+1,nBasis+nbeta]),&
      [1,nBasis*2],[nalpha+1,nalpha+nbeta])
  end if

  call MO_I%print(6,'MO_I')
  MO_J = MO_I

  call MIJ%initialize(nElectrons,nElectrons,(0.0,0.0))
  call MIJ_inv%initialize(nElectrons,nElectrons,(0.0,0.0))

  MIJ = matmul(dagger(MO_I),matmul(overlap,MO_J))
  call MIJ%print(6,'MIJ')
  MIJ_inv = MIJ
  rho = matmul(MO_J,matmul(MIJ_inv,dagger(MO_I)))
  call rho%print(6,'rho')

  call Idem%initialize(nBasis*2,nBasis*2,(0.0,0.0))
  call ps%initialize(nBasis*2,nBasis*2,(0.0,0.0))
  Idem = matmul(matmul(rho,overlap),matmul(rho,overlap)) - matmul(rho,overlap)
  ps = matmul(rho,overlap)
  ps_trace = ps%trace()
  call Idem%print(6,'Idempotence Check')
  call ps_trace%print(6, 'PS Trace Check')

  core_con = Contraction(rho,core_ham)
  call core_con%print(6,'Core Contraction',.true.)
! intitialize G matrix and contract ERI w/ rho

  call Gmat%initialize(nBasis*2,nBasis*2,(0.0,0.0))
  call p_Fock%initialize(nBasis*2,nBasis*2,(0.0,0.0))

  DO I = 1, NBASIS
    DO J = 1, NBASIS
      tAA=0;tBB=0;tAB=0;tBA=0
      DO K = 1, NBASIS
        DO L = 1, NBASIS
          !alpha alpha
          tAA = tAA + ERI%at(I,J,K,L)*rho%at(K,L)
          tAA = tAA + ERI%at(I,J,K,L)*rho%at(K+nBasis,L+nBasis)
          tAA = tAA - ERI%at(I,L,K,J)*rho%at(K,L)
          call Gmat%put(tAA,I,J)

          !beta beta
          tBB = tBB + ERI%at(I,J,K,L)*rho%at(K+nBasis,L+nBasis)
          tBB = tBB + ERI%at(I,J,K,L)*rho%at(K,L)
          tBB = tBB - ERI%at(I,L,K,J)*rho%at(K+nBasis,L+nBasis)
          call Gmat%put(tBB,I+nBasis,J+nBasis)

          !alpha beta
          tAB = tAB - ERI%at(I,L,K,J)*rho%at(K+nBasis,L)
          call Gmat%put(tAB,I+nBasis,J)

          !beta alpha
          tBA = tBA - ERI%at(I,L,K,J)*rho%at(K,L+nBasis)
          call Gmat%put(tBA,I,J+nBasis)
        END DO
      END DO
    END DO
  END DO
 
  call Gmat%print(6,'G Matrix')
  p_Fock = core_ham + Gmat
  call p_Fock%print(6,'P Fock')
  eri_con = Contraction(Gmat,rho)
  call eri_con%print(6,'ERI contraction',.true.)

  HIJ = core_con + (half*eri_con) + Vnn
  call HIJ%print(6,'HIJ',.true.)

end program GeneralizedHF
