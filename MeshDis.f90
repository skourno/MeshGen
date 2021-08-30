Module MeshDis

    Type Map
       Integer   ::  n, s, e, w     ! # of north-south-east-west neighbor
    End Type Map

    Type(Map), Allocatable  ::  kMap(:)    ! k "map" in the mesh
    Real*8   , Allocatable  ::  k(:)       ! k values - numbering is correlated
                                           ! to the appropriate nodes, through kMap

contains

    !**************************************************************************
    ! Used as a master subroutine to perform Mesh Displacment. Call this in
    ! main program after mesh generation to move it based on a displacement
    ! input file
    !**************************************************************************
    Subroutine DisplaceMesh
       Use System,     Only: ReadBoundaryDisplacement, PrintSystem, MyGrid, Node
       Implicit None
       Integer         :: nk, i, INode

       ! Read Input Displacement File
       ! -----------------------------
       call ReadBoundaryDisplacement
       print*, '-------------------------------------'
       print*, ' UPDATED BOUNDARY NODE COORDINATES ACCORDING TO INPUT FILE '

      ! Read West Boundary Node Coordinates Clockwise
      INode = 1
      do i = 1, MyGrid%NNodes(2)
         write(44,*) Node(INode)%R(:)
         INode            = INode + 1
      enddo

      ! Read North Boundary Node Coordinates Clockwise
      INode = 2*MyGrid%NNodes(2)
      do i = 2, MyGrid%NNodes(1)
         write(44,*) Node(INode)%R(:)
         INode     = INode + MyGrid%NNodes(2)
      enddo

      ! Read East Boundary Node Coordinates Clockwise
      INode = MyGrid%TotalNodes - 1
      do i = 2, MyGrid%NNodes(2)
         write(44,*) Node(INode)%R(:)
         INode     = INode - 1
      enddo

      ! Read South Boundary Node Coordinates Clockwise
      INode = MyGrid%TotalNodes - 2*MyGrid%NNodes(2) + 1
      do i = 2, MyGrid%NNodes(1) - 1
         write(44,*) Node(INode)%R(:)
         INode     = INode - MyGrid%NNodes(2)
      enddo


       ! Allocate module arrays
       ! -----------------------
       Allocate ( kMap (MyGrid%TotalNodes) )

       nk  =  ( MyGrid%NNodes(1) - 2 )*( MyGrid%NNodes(2) - 1 ) + &
              ( MyGrid%NNodes(2) - 2 )*( MyGrid%NNodes(1) - 1 )

       Allocate ( k ( nk ) )


       ! Calculate k values
       ! ------------------
       call Calculate_k_values ( nk )

       ! Solve and calculate new node coordinates
       ! ----------------------------------------
       call MSIP_k_Interface ( 0.5d0 )

       ! Output Displaced Grid
       ! ---------------------
       call PrintSystem ( 24 )

    End Subroutine DisplaceMesh


    !**************************************************************************
    ! Used to calculate and "map" k spring constant for every internal Node.
    ! Can be used to initialize or update k values. Input # of k values within
    ! the mesh
    !**************************************************************************
    Subroutine Calculate_k_values ( Nk )
       Use System,     Only: MyGrid, Node, NOP
       Implicit None
       Integer         :: INode, JNode, Ik, Nk, i, j
       Real*8          :: dist
       Intent(IN)      :: Nk

       ! Initialize k
       kMap(:)%n = 0
       kMap(:)%s = 0
       kMap(:)%w = 0
       kMap(:)%e = 0
       k(:)      = 0
       Ik        = 1

       ! Loop over all interior Nodes and calculate north and east k values
       do i = 2, MyGrid%NNodes(2) - 1
           do j = 2, MyGrid%NNodes(1) - 1
               INode  = NOP(i,j)

               ! North
               JNode         =  Node(INode)%North
               dist          =  DISTANCE ( Node(INode)%R(:) , Node(JNode)%R(:) )
               k(Ik)         =  1 / dist**(0.1d0)
               kMap(INode)%n =  Ik
               Ik            =  Ik + 1

               ! East
               JNode         =  Node(INode)%East
               dist          =  DISTANCE ( Node(INode)%R(:) , Node(JNode)%R(:) )
               k(Ik)         =  1 / dist**(0.1d0)
               kMap(INode)%e =  Ik
               Ik            =  Ik + 1
           enddo
       enddo

       ! The only k values that need to be calculated are the west and south boundary

       ! Start by calculating West Boundary
       do INode = 2, MyGrid%NNodes(2) - 1
           JNode         =  Node(INode)%East
           dist          =  DISTANCE ( Node(INode)%R(:) , Node(JNode)%R(:) )
           k(Ik)         =  1 / dist**(0.1d0)
           kMap(INode)%e =  Ik
           Ik            =  Ik + 1
       enddo

       ! Calculate South Boundary
       INode = MyGrid%TotalNodes - 2*MyGrid%NNodes(2) + 1
       do i = 2, MyGrid%NNodes(1) - 1
           JNode         =  Node(INode)%North
           dist          =  DISTANCE ( Node(INode)%R(:) , Node(JNode)%R(:) )
           k(Ik)         =  1 / dist**(0.1d0)
           kMap(INode)%n =  Ik
           Ik            =  Ik + 1

           INode         =  INode - MyGrid%NNodes(2)
       enddo

       ! Check that k array is filled and each has been calculated once
       ! if ( PRODUCT( k(:) ) == 0 ) STOP ' ERROR WITH K ARRAY - LOCATED A ZERO VALUE WITHIN K '
       if ( Nk.NE.( Ik - 1 ) )     STOP ' ERROR WITH # OF K VALUES - CALCULATED WRONG # OF K VALUES '


       ! All Nodes have been assigned a northern and an eastern k value. Based on this
       ! the remaining map is completed.

       ! Complete East Boundary Mapping
       INode = MyGrid%TotalNodes - 1
       do i = 2, MyGrid%NNodes(2) - 1
           JNode         =  Node(INode)%West
           kMap(INode)%w =  kMap(JNode)%e
           INode         =  INode - 1
       enddo

       ! Complete North Boundary Mapping
       INode = 2*MyGrid%NNodes(2)
       do i = 2, MyGrid%NNodes(1) - 1
           JNode         =  Node(INode)%South
           kMap(INode)%s =  kMap(JNode)%n
           INode         =  INode + MyGrid%NNodes(2)
       enddo

       ! Loop over all interior
       do i = 2, MyGrid%NNodes(2) - 1
           do j = 2, MyGrid%NNodes(1) - 1
               INode  =  NOP(i,j)

               ! Assign South
               JNode         =  Node(INode)%South
               kMap(INode)%s =  kMap(JNode)%n

               ! Assign West
               JNode         =  Node(INode)%West
               kMap(INode)%w =  kMap(JNode)%e

           enddo
       enddo

    End Subroutine Calculate_k_values


    !**************************************************************************
    ! Used for caclulating the new coordinates for the displaced system using
    ! MSIP.
    !**************************************************************************
    Subroutine MSIP_k_Interface ( psi )
       Use System,   Only: MyGrid
       Use GridGen,  Only: MSIP9
       Implicit None
       Integer       :: iter, itermax, FirstIter_Flag
       Real*8        :: psi, Residual
       Intent(IN)    :: psi

       itermax  =  200

       print*, '-------------------------------------'
       print*, ' INITIATING MSIP ITERATIONS '
       print*, '-------------------------------------'
       FirstIter_Flag  =  1

       do iter = 1, itermax
           call Create_kStencil ( FirstIter_Flag )
           call MSIP9 ( psi , MyGrid%TotalNodes , Residual )

           print*, 'Iteration : ', iter
           print*, 'Residual  : ', Residual
           print*, '---------------------------------'

           if ( Residual <= 1e-4 ) EXIT
           FirstIter_Flag = 0
       enddo

    EndSubroutine MSIP_k_Interface



    !**************************************************************************
    ! Used to create the stencil arrays Qr. Results are imported to module GridGen
    ! ready for use by MSIP9 subroutine. IFlag = 1 (iter = 1), else IFlag = 0
    !**************************************************************************
    Subroutine Create_kStencil ( IFlag )
       Use System,   Only: MyGrid, Node, NOP
       Use GridGen,  Only: Qr
       Implicit None
       Integer       :: IFlag, INode, i, j
       Integer       :: n, s, e, w
       Intent(IN)    :: IFlag

       ! for iter = 1 Qr should be allocated
       if (IFlag == 1) Allocate ( Qr(MyGrid%TotalNodes) )

       ! Initialize Qr
       do INode = 1,  MyGrid%TotalNodes
           Qr(INode)%n      =  0.d0
           Qr(INode)%s      =  0.d0
           Qr(INode)%w      =  0.d0
           Qr(INode)%e      =  0.d0
           Qr(INode)%ne     =  0.d0
           Qr(INode)%nw     =  0.d0
           Qr(INode)%se     =  0.d0
           Qr(INode)%sw     =  0.d0
           Qr(INode)%p      =  0.d0

           Qr(INode)%rhs(:) =  0.d0
       enddo


       ! Loop over all interior nodes
       do i = 2, MyGrid%NNodes(2) - 1
           do j = 2, MyGrid%NNodes(1) - 1
               INode  =  NOP(i,j)

               n      =  kMap(INode)%n
               s      =  kMap(INode)%s
               e      =  kMap(INode)%e
               w      =  kMap(INode)%w

               if (n*s*e*w == 0) STOP '!!! INTERNAL NODE MAPPING ERROR !!!'

               Qr(INode)%p  =  k(n) + k(s) + k(w) + k(e)
               Qr(INode)%n  =  k(n)
               Qr(INode)%s  =  k(s)
               Qr(INode)%w  =  k(w)
               Qr(INode)%e  =  k(e)
               Qr(INode)%ne =  0.d0
               Qr(INode)%nw =  0.d0
               Qr(INode)%se =  0.d0
               Qr(INode)%sw =  0.d0
           enddo
       enddo

       ! West Boundary Conditions
       INode = 1
       do i = 1, MyGrid%NNodes(2)
           Qr(INode)%n      =  0.d0
           Qr(INode)%s      =  0.d0
           Qr(INode)%w      =  0.d0
           Qr(INode)%e      =  0.d0
           Qr(INode)%ne     =  0.d0
           Qr(INode)%nw     =  0.d0
           Qr(INode)%se     =  0.d0
           Qr(INode)%sw     =  0.d0
           Qr(INode)%p      =  1.d0

           Qr(INode)%rhs(:) =  Node(INode)%R(:)
           INode            =  INode + 1
       enddo


       ! North Boundary Conditions
       INode = 2*MyGrid%NNodes(2)
       do i = 2, MyGrid%NNodes(1)
           Qr(INode)%n      =  0.d0
           Qr(INode)%s      =  0.d0
           Qr(INode)%w      =  0.d0
           Qr(INode)%e      =  0.d0
           Qr(INode)%ne     =  0.d0
           Qr(INode)%nw     =  0.d0
           Qr(INode)%se     =  0.d0
           Qr(INode)%sw     =  0.d0
           Qr(INode)%p      =  1.d0

           Qr(INode)%rhs(:) =  Node(INode)%R(:)
           INode     = INode + MyGrid%NNodes(2)
       enddo

       ! East Boundary Conditions
       INode = MyGrid%TotalNodes - 1
       do i = 2, MyGrid%NNodes(2)
           Qr(INode)%n      =  0.d0
           Qr(INode)%s      =  0.d0
           Qr(INode)%w      =  0.d0
           Qr(INode)%e      =  0.d0
           Qr(INode)%ne     =  0.d0
           Qr(INode)%nw     =  0.d0
           Qr(INode)%se     =  0.d0
           Qr(INode)%sw     =  0.d0
           Qr(INode)%p      =  1.d0

           Qr(INode)%rhs(:) =  Node(INode)%R(:)
           INode            =  INode - 1
       enddo

      ! South Boundary Conditions
      INode = MyGrid%TotalNodes - 2*MyGrid%NNodes(2) + 1
      do i = 2, MyGrid%NNodes(1) - 1
           Qr(INode)%n      =  0.d0
           Qr(INode)%s      =  0.d0
           Qr(INode)%w      =  0.d0
           Qr(INode)%e      =  0.d0
           Qr(INode)%ne     =  0.d0
           Qr(INode)%nw     =  0.d0
           Qr(INode)%se     =  0.d0
           Qr(INode)%sw     =  0.d0
           Qr(INode)%p      =  1.d0

           Qr(INode)%rhs(:) =  Node(INode)%R(:)
           INode            = INode - MyGrid%NNodes(2)
      enddo


    EndSubroutine Create_kStencil





    !**************************************************************************
    ! Given the vectors of 2 points in 2D space, returns the distance between them
    !**************************************************************************
    Real*8 Function DISTANCE ( R1, R2 )
       Implicit None
       Real*8          ::  R1(2), R2(2), a, b

       a        =  R1(1) - R2(1)
       b        =  R1(2) - R2(2)

       Distance =  SQRT( a*a + b*b )

    End function



End Module MeshDis
