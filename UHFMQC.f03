program uhfMQC

  use mqc_gaussian
  use mqc_algebra
  use iso_fortran_env

  implicit none
  character(len=:),allocatable::filename
  type(mqc_gaussian_unformatted_matrix_file)::matFile
  type(mqc_matrix)::MO_IA,MO_IB,MO_JA,MO_JB,core_ham,MIJ_alpha,MIJ_beta,&
    rho_alpha,rho_beta,overlap,rho_total,R2_ERI,G_alpha,G_beta,p_Fock
  type(mqc_molecule_data)::mol_info
  type(mqc_scalar)::NIJ,HIJ,core_con,eri_con,Vnn,temp_a,temp_b,two,half
  type(mqc_R4Tensor)::ERI
  integer::I,J,K,L,nElectrons,nBasis,nAlpha,nBeta,iOut=6

  !Variable assignments
  two = 2.0;half=0.5
  call mqc_get_command_argument(1,filename)
  call matFile%load(filename)
  call matFile%getMolData(mol_info)
  nElectrons = matFile%nelectrons
  nBasis = matFile%nbasis
  nBeta = nElectrons/2
  nAlpha = nElectrons - nBeta
  Vnn = mqc_get_nuclear_repulsion(iOut,mol_info)
  call matFile%getArray('OVERLAP',matrixOut=overlap)
  call matFile%getArray('ALPHA MO COEFFICIENTS',matrixOut=MO_IA)
  call matFile%getArray('BETA MO COEFFICIENTS',matrixOut=MO_IB)
  call matFile%getArray('CORE HAMILTONIAN ALPHA',matrixOut=core_ham)
  call matFile%getArray('REGULAR 2E INTEGRALS',r4TensorOut=ERI)
  MO_JA = MO_IA
  MO_JB = MO_IB
  MO_IA = MO_IA%mat([1,nBasis],[1,nAlpha])
  MO_IB = MO_IB%mat([1,nBasis],[1,nBeta])
  MO_JA = MO_JA%mat([1,nBasis],[1,nAlpha])
  MO_JB = MO_JB%mat([1,nBasis],[1,nBeta])
  
  MIJ_alpha = matmul(transpose(MO_IA),matmul(overlap,MO_JA))
  MIJ_beta = matmul(transpose(MO_IB),matmul(overlap,MO_JB))
  rho_alpha = matmul(matmul(MO_JA,MIJ_alpha),transpose(MO_IA))
  rho_beta = matmul(matmul(MO_JB,MIJ_beta),transpose(MO_IB))
  rho_total = rho_alpha + rho_beta
  core_con = half*(contraction(rho_total,core_ham) + contraction(rho_alpha,core_ham) &
    + contraction(rho_beta,core_ham))
  
  call G_alpha%initialize(nBasis,nBasis,0.0)
  call G_beta%initialize(nBasis,nBasis,0.0)
  DO I = 1, nBasis
    DO J = 1, nBasis
      temp_a = 0.0
      temp_b = 0.0
      DO K = 1, nBasis
        DO L = 1, nBasis
          temp_a = temp_a + (ERI%at(I,J,K,L) * rho_total%at(K,L))
          temp_a = temp_a - (ERI%at(I,K,J,L) * rho_alpha%at(K,L))
          temp_b = temp_b + (ERI%at(I,J,K,L) * rho_total%at(K,L))
          temp_b = temp_b - (ERI%at(I,K,J,L) * rho_beta%at(K,L))
          call G_alpha%put(temp_a,I,J)
          call G_beta%put(temp_b,I,J)
        END DO
      END DO
    END DO
  END DO

  call p_Fock%initialize(nBasis,nBasis,0.0)
  DO I = 1, nBasis
    DO J = 1, nBasis
      call p_Fock%put((G_alpha%at(I,J) + core_ham%at(I,J)),I,J)
    END DO
  END DO

  eri_con = contraction(G_alpha,rho_alpha) + contraction(G_beta,rho_beta)

  HIJ = core_con + (half*eri_con) + Vnn

  !Diagnostics printing
  print *,'NElectrons: ', nElectrons
  print *,'nAlpha: ', nAlpha, ' nBeta: ', nBeta
  call overlap%print(iOut,'Overlap')
  call MO_IA%print(iOut,'Alpha MO Coefficients')
  call MO_IB%print(iOut,'Beta MO Coefficients')
  call MIJ_alpha%print(iOut,'MIJ Alpha')
  call MIJ_beta%print(iOut,'MIJ Beta')
  call rho_alpha%print(iOut,'Alpha Density')
  call rho_beta%print(iOut,'Beta Density')
  call rho_total%print(iOut,'Total Density')
  call core_ham%print(iOut,'Core Hamiltonian')
  call core_con%print(iOut,'Core Contraction',.true.,.true.)
  call G_alpha%print(iOut,'G Alpha')
  call p_Fock%print(iOut,'Psuedo Fock Alpha')
  call eri_con%print(iOut,'ERI Contraction',.true.,.true.)
  call Vnn%print(iOut,'VNN',.true.,.true.)
  call HIJ%print(iOut,'HIJ',.true.,.true.)
end program uhfMQC
