program GeneralizedHF

  use mqc_gaussian
  use mqc_algebra
  use iso_fortran_env

  implicit none
  type(mqc_gaussian_unformatted_matrix_file)::matFile
  type(mqc_wavefunction)::wavefunc
  type(mqc_molecule_data)::molInfo
  type(mqc_matrix)::MO_I,MO_J,rho,MIJ,MIJ_inv,MOcoeff,overlap,core_ham,Gmat,Fock,p_Fock,dens,&
    dens_array,core_array,MO_J_dagger,rho_alt,Idem,ps,real_G
  type(mqc_r4tensor)::ERI
  type(mqc_scalar)::half,two,core_con,eri_con,tAA,tBB,tAB,tBA,Vnn,HIJ,alt_con,array_con,ps_trace, &
    eri_alt_con,HIJ_alt
  character(len=:),allocatable::fileName
  integer::i,j,k,l,nelectrons,nalpha,nbeta,nbasis

  call setvbuf3f(6, 2, 0)
  call mqc_get_command_argument(1,fileName)
  call matFile%load(fileName)
  call matFile%getESTObj('wavefunction',wavefunc)
  call matFile%getArray('REGULAR 2E INTEGRALS',r4TensorOut=ERI)
  call matFile%getMolData(molInfo)
  call overlap%initialize(nBasis*2,nBasis*2,(0.0,0.0))
  nelectrons = wavefunc%nElectrons%ival
  nalpha = wavefunc%nAlpha%ival
  nbeta = wavefunc%nBeta%ival
  nbasis = wavefunc%nbasis%ival
  overlap = wavefunc%overlap_matrix%getBlock('full')
  dens = wavefunc%density_matrix%getBlock('full')
  call matFile%getArray('ALPHA DENSITY MATRIX',matrixOut=dens_array)
  call matFile%getARray('CORE HAMILTONIAN ALPHA',matrixOut=core_array)
  call overlap%print(6,'Overlap')
  
  MOcoeff = wavefunc%mo_coefficients%getBlock('full')
  core_ham = wavefunc%core_hamiltonian%getBlock('full')
  Fock = wavefunc%fock_matrix%getBlock('full')
  Vnn = mqc_get_nuclear_repulsion(6,molInfo)
  half = 0.5
  two = 2.0

  nalpha = 8
  nbeta = 8
!  call overlap%print(6,'Overlap')
  call MOcoeff%print(6,'MO Coefficients')
  
!  initialize occupied MO coefficient matrices (2*nBasis,nAlpha+nBeta)
  
  call core_ham%print(6,'Core Hamiltonian')
  call MO_I%initialize(nBasis*2,nElectrons)
  call MO_J%initialize(nBasis*2,nElectrons)

!  propogate MO_I matrices for alpha/beta electrons

  DO I = 1, NALPHA
    call MO_I%mput(MOcoeff%mat([1,nBasis*2],[I,I]),[1,nBasis*2],[I,I])
  END DO

  DO I = 1, NBETA
    call MO_I%mput(MOcoeff%mat([1,nBasis*2],[nBasis+I,nBasis+I]),[1,nBasis*2],[nAlpha+I,nAlpha+I])
  END DO

  call MO_I%print(6,'MO_I')
  MO_J = MO_I
  MO_J_dagger = dagger(MO_J)
  call MO_J_dagger%print(6,'MO_J Dagger')

  MIJ = matmul(dagger(MO_I),matmul(overlap,MO_J))
  call MIJ%print(6,'MIJ')

  call ERI%print(6,'ERI output')
  MIJ_inv = MIJ%inv()
  rho = matmul(MO_J,matmul(MIJ_inv,dagger(MO_I)))
  rho_alt = matmul(MO_J,dagger(MO_I))
  call rho%print(6,'rho')
  call rho_alt%print(6,'rho_alt')
  call dens%print(6,'density')
  call dens_array%print(6,'density array')

  Idem = matmul(matmul(rho,overlap),matmul(rho,overlap)) - matmul(rho,overlap)
  ps = matmul(rho,overlap)
  ps_trace = ps%trace
  call Idem%print(6,'Idempotence Check')
  call ps_trace%print(6, 'PS Trace Check')
  
  Idem = matmul(matmul(dens,overlap),matmul(dens,overlap)) - matmul(dens,overlap)
  call Idem%print(6,'Second Idempotence Check')
  ps = matmul(dens,overlap)
  ps_trace = ps%trace
  call ps_trace%print(6,'Second PS Trace Check')

  core_con = Contraction(rho,core_ham)
  array_con = Contraction(dens_array,core_array)
  call core_con%print(6,'Core Contraction',.true.)
  call array_con%print(6,'Array Contraction',.true.)
! intitialize G matrix and contract ERI w/ rho

  call Gmat%initialize(nBasis*2,nBasis*2)
  call p_Fock%initialize(nBasis*2,nBasis*2,(0.0,0.0))

  print *, 'Before contraction loop'
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
 !         tAB = tAB + ERI%at(I,J,K,L)*rho%at(K+nBasis,L)
          tAB = tAB - ERI%at(I,L,K,J)*rho%at(K+nBasis,L)
          call Gmat%put(tAB,I+nBasis,J)

          !beta alpha
 !         tBA = tBA + ERI%at(I,J,K,L)*rho%at(K,L+nBasis)
          tBA = tBA - ERI%at(I,L,K,J)*rho%at(K,L+nBasis)
          call Gmat%put(tBA,I,J+nBasis)
        END DO
      END DO
    END DO
  END DO
  
  real_G = Fock - core_ham
  call Gmat%print(6,'G Matrix')
  call real_G%print(6,'Real G Matrix')
  p_Fock = core_ham + Gmat
  call p_Fock%print(6,'P Fock')
  call Fock%print(6,'Fock')
  eri_con = Contraction(Gmat,rho)
  eri_alt_con = Contraction((Fock - core_ham),dens)
  call eri_con%print(6,'ERI contraction',.true.)
  call eri_alt_con%print(6,'ERI alt contraction',.true.)

  HIJ = core_con + (half*eri_con) + Vnn
  HIJ_alt = array_con + (half*eri_alt_con) + Vnn 
  call HIJ%print(6,'HIJ',.true.)
  call HIJ_alt%print(6,'HIJ alt',.true.)

end program GeneralizedHF
