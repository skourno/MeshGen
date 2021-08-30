program MainProg
   Use System ,      Only: ReadInput, PrintSystem, MyGrid, NOP
   Use GridGen,      Only: InitializeMesh, CalculateMetrics, MSIP_Interface, InputSorenson
   Use MeshDis,      Only: DisplaceMesh
   Implicit None
   Integer           ::  INode, i, j, iter, itermax, IMode
   Real*8            ::  Residual
   Logical           ::  LDis

   print*, ''
   call ReadInput
   print*, 'READ INPUT FILE'

   call InitializeMesh( 1.15 , 1.1 )
   print*, 'INITIALIZED MESH'
   print*, '-------------------------------------'
   print*, 'INPUT WORKING MODE'
   print*, '1 - Laplace | 2 - Thomas-Middlecoff | 3 - Sorenson '
   read(*,*) IMode
   print*, '-------------------------------------'

   if (IMode == 3) call InputSorenson

   itermax = 100
   do iter = 1, itermax
       print*, 'Initiate MSIP Iteration:', iter

       do i = 2, MyGrid%NNodes(2) - 1
           do j = 2, MyGrid%NNodes(1) - 1
               INode  = NOP(i,j)
               call   CalculateMetrics( INode , 0 )    ! 0 flags the node as interior
           enddo
       enddo

       print*, 'Calculated Node Metrics'

       Call MSIP_Interface ( IMode, 0.5d0, Residual )
       print*, 'Residual:', Residual
       print*, '---------------------------------'

       ! Define exit criteria here
       if (Residual <= 1e-4) EXIT
   enddo

   call PrintSystem ( 23 )

   print*, 'GENERATION SUCCESS'
   print*, 'YOUR GRID HAS BEEN SAVED @ MYGRID.txt'

   print*, '-------------------------------------'
   print*, 'DO YOU WANT TO USE MESH DISPLACEMENT MOD ? '
   print*, '            T - yes , F - no               '
   read*,   LDis

   if ( LDis .EQV. .TRUE. ) then
      call DisplaceMesh
      print*, ' DISPLACEMENT SUCCESFUL '
      print*, ' NEW NODE COORDINATES SAVED @ DISGRID.TXT '
   endif

end



   !*************************************
   ! DEBUG HELP
   !*************************************
   ! Remember to include Node variable if used in main


   !PRINT METRICS OF INTERIOR NODES
   !-------------------------------
   !do i = 2, MyGrid%NNodes(2) - 1
   !    do j = 2, MyGrid%NNodes(1) - 1
   !        INode  = NOP(i,j)
   !        write(3,*) metrics(INode)%Jac
   !    enddo
   !enddo
